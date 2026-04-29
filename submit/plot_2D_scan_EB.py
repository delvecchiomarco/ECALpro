import ROOT
import math
import argparse
import numpy as np
import time
import csv


ROOT.gROOT.SetBatch(True)  # no graphical output
ROOT.EnableImplicitMT()  # multithreading, internal parallelization
ROOT.TGaxis.SetMaxDigits(4)

nbins = 120
xmin ,xmax = 0.05, 0.25  #histograms range
mean_min, mean_max = 0.115, 0.140  # mean range
sigma_min, sigma_max = 0.015, 0.025  # sigma range
fit_min, fit_max = 0.09, 0.17  # fit window
nsigma = 3.0  # for significance computation
pt_min_vals = np.arange(0.5, 0.75 + 1e-6, 0.01)
S4S9_low_min_vals = np.arange(0.75, 0.90 + 1e-6, 0.01)
S4S9_high_min_vals = np.arange(0.75, 0.90 + 1e-6, 0.01)

REGIONS = {
  "EB_low": {
    "eta_expr": "abs_etaPi0_cor < 1.0",
    "ptG1_min": 0.70,  # before was 0.65
    "ptG2_min": 0.70,  # before was 0.65
    "ptPi0_min": 1.50,  # before was 2.0
    "S4S9_1_min": 0.65,  # before was 0.88
    "S4S9_2_min": 0.65,  # before was 0.88
    "HLTIso_max": 0.5,
  },
  "EB_high": {
    "eta_expr": "abs_etaPi0_cor > 1.0 && abs_etaPi0_cor < 1.479",
    "ptG1_min": 0.70,  # before was 0.65
    "ptG2_min": 0.68,  # before was 0.65
    "ptPi0_min": 1.50,  # before was 1.75
    "S4S9_1_min": 0.65,  # before was 0.92
    "S4S9_2_min": 0.65,  # before was 0.92
    "HLTIso_max": 0.5,
  }
}


def make_mask(cfg):
    return (
        f"({cfg['eta_expr']})"
        f" && (ptG1_cor > {cfg['ptG1_min']})"
        f" && (ptG2_cor > {cfg['ptG2_min']})"
        f" && (ptPi0_cor > {cfg['ptPi0_min']})"
        f" && (HLTIsoPi0 < {cfg['HLTIso_max']})"
        f" && (S4S9_1 > {cfg['S4S9_1_min']})"
        f" && (S4S9_2 > {cfg['S4S9_2_min']})"
    )

def style_hist(h, x_title="m_{#pi^{0}} [GeV]", y_title=None,
               marker_style=20, marker_size=0.5, line_width=1):
    h.SetStats(0)
    h.GetXaxis().SetTitle(x_title)
    if y_title is None:
        y_title = f"Entries / {h.GetBinWidth(1)*1000:.1f} MeV"
    h.GetYaxis().SetTitle(y_title)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetMaxDigits(4)
    h.SetMarkerStyle(marker_style)
    h.SetMarkerSize(marker_size)
    h.SetLineWidth(line_width)

def fill_hist(df, colname, hname="", htitle="", nbins=120, xmin=0.05, xmax=0.25):
    hptr = df.Histo1D((hname, htitle, nbins, xmin, xmax), colname)
    h = hptr.GetValue()
    h.SetDirectory(0)
    h.Sumw2()
    style_hist(h)
    return h

def save_hist(h, outdir, fout, plot_name=""):
    c = ROOT.TCanvas(f"c_{plot_name}", "", 1200, 1000)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.1)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.1)
    h.Draw("E1")
    c.SaveAs(outdir + plot_name + ".pdf")
    c.SaveAs(outdir + plot_name + ".png")
    fout.cd()
    h.Write(h.GetName(), ROOT.TObject.kOverwrite)

def initial_guesses(hist, xlo, xhi):
    ax = hist.GetXaxis()
    b1 = ax.FindBin(xlo)
    b2 = ax.FindBin(xhi)
    max_y = -1.0
    max_x = 0.135
    for b in range(b1, b2 + 1):
        y = hist.GetBinContent(b)
        x = ax.GetBinCenter(b)
        if y > max_y:
            max_y = y
            max_x = x
    return max_y if max_y > 0 else 1.0, max_x, 0.004

def pol2_fit(hist, func):
    amp0, mean0, sig0 = initial_guesses(hist, mean_min, mean_max)
    # fit: gaus + pol2
    func.SetParameters(amp0, mean0, sig0, 1.0, 0.0, 0.0)
    func.SetParLimits(0, 0.0, 1e9)  # gaussian amplitude >= 0
    func.SetParLimits(1, mean_min, mean_max)  # mean limits [0.115, 0.140] GeV
    func.SetParLimits(2, sigma_min, sigma_max)  # sigma limits [15, 25] MeV
    # "R": range, "I": integral bins, "S": FitResult
    res_pol2 = hist.Fit(func, "R I S 0")
    return res_pol2

def compute_chi2_ndf(res_pol2_status, func):
    status = int(res_pol2_status.Status())
    if status != 0:
        print(f"[WARN] Fit status = {status}")
    chi2 = func.GetChisquare()
    ndf  = func.GetNDF()
    return chi2, ndf

def extract_fit_results(name, func, n_sigma, binw):
    mean  = func.GetParameter(1)
    sigma = func.GetParameter(2)
    x1 = mean - n_sigma * sigma
    x2 = mean + n_sigma * sigma
    fs = ROOT.TF1(f"fs_pol2_{name}", "gaus(0)", x1, x2)  # signal
    fb = ROOT.TF1(f"fb_pol2_{name}", "pol2(0)", x1, x2)  # background
    for ip in range(3):
        fs.SetParameter(ip, func.GetParameter(ip))
        fb.SetParameter(ip, func.GetParameter(3 + ip))
    # statistical significance
    Nsig = fs.Integral(x1, x2) / binw
    Nbkg = fb.Integral(x1, x2) / binw
    Z = Nsig / math.sqrt(Nsig + Nbkg) if (Nsig > 0 and (Nsig + Nbkg) > 0) else 1e-6
    return mean, sigma, x1, x2, Nsig, Nbkg, Z

def save_hist_postfit(h, func, outdir, fout, plot_name=""):
    c = ROOT.TCanvas(f"c_{plot_name}", "", 1200, 1000)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.1)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.1)
    h.Draw("E1")
    # vertical lines for fit range
    l1 = ROOT.TLine(fit_min, 0, fit_min, h.GetMaximum()*1.05)
    l2 = ROOT.TLine(fit_max, 0, fit_max, h.GetMaximum()*1.05)
    l1.SetLineStyle(2)
    l2.SetLineStyle(2)
    l1.Draw("same")
    l2.Draw("same")
    func.SetLineWidth(1)
    func.Draw("same")
    # fit results
    pave = ROOT.TPaveText(0.28, 0.20, 0.58, 0.50, "NDC")
    pave.SetFillColor(0)
    pave.SetBorderSize(0)
    pave.SetTextFont(42)
    pave.SetTextAlign(12)
    pave.AddText("Fit: gaus + pol2")
    pave.AddText(f"A = {func.GetParameter(0):.2e} #pm {func.GetParError(0):.2e}")
    pave.AddText(f"#mu = {func.GetParameter(1)*1000:.2f} #pm {func.GetParError(1)*1000:.2f} MeV")
    pave.AddText(f"#sigma = {func.GetParameter(2)*1000:.2f} #pm {func.GetParError(2)*1000:.2f} MeV")
    pave.AddText(f"pol2[0] = {func.GetParameter(3):.2e} #pm {func.GetParError(3):.2e}")
    pave.AddText(f"pol2[1] = {func.GetParameter(4):.2e} #pm {func.GetParError(4):.2e}")
    pave.AddText(f"pol2[2] = {func.GetParameter(5):.2e} #pm {func.GetParError(5):.2e}")
    pave.AddText(f"#chi^{{2}} / ndf = {func.GetChisquare():.1f} / {func.GetNDF()}")
    pave.Draw()
    c.SaveAs(outdir + plot_name + ".pdf")
    c.SaveAs(outdir + plot_name + ".png")
    fout.cd()
    h.Write(h.GetName(), ROOT.TObject.kOverwrite)

def print_fit_results(chi2, ndf, mean, sigma, x1, x2, Nsig, Nbkg, Z, n_sigma):
    print("\n=== Fit results (range {:.3f}-{:.3f} GeV) ===".format(fit_min, fit_max))
    print("Model: gaus + pol2")
    print(f"  chi2/ndf = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}")
    print(f"  mean = {mean:.6f} GeV  sigma = {sigma*1000:.2f} MeV")
    print(f"  window = [{x1:.6f}, {x2:.6f}] (+/- {n_sigma:.0f} sigma)")
    print(f"  Nsig = {Nsig:.0f}  Nbkg = {Nbkg:.0f}  Z = {Z:.2f}")

def compute_Z_for_ptcuts(name, df_base, base_cfg, pt1, pt2, tag, counter, outdir, fout):
    cfg = dict(base_cfg)
    cfg["ptG1_min"] = float(pt1)
    cfg["ptG2_min"] = float(pt2)
    mask = make_mask(cfg)
    df_tmp = df_base.Define(f"mPi0_scan", f"mPi0_cor[{mask}]")
    h = fill_hist(df_tmp, f"mPi0_scan", f"h_scan_{name}_{tag}", htitle="")
    if h.GetEntries() < 50:
        return 1e-6, 1e-6, 1e-6
    f_pol2_scan = ROOT.TF1(f"f_pol2_{name}_scan_{tag}", "gaus(0)+pol2(3)", fit_min, fit_max)
    res_pol2_scan = pol2_fit(h, f_pol2_scan)
    status = int(res_pol2_scan.Status())
    if status != 0:
        return 1e-6, 1e-6, 1e-6
    if counter % 4 == 0:
        save_hist_postfit(h, f_pol2_scan, outdir, fout, f"mpi0_{name}_postfit_{tag}")
    mean, sigma, x1, x2, Nsig, Nbkg, Z = extract_fit_results(f"{name}_scan_{tag}", f_pol2_scan, nsigma, h.GetBinWidth(1))
    return Nsig, Nbkg, Z

def hist_scan_options(hist, x_title, y_title, z_title):
    hist.SetStats(0)
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetZaxis().SetTitle(z_title)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetZaxis().SetTitleOffset(1.2)
    hist.SetContour(256)
    hist.GetZaxis().SetNdivisions(510)
    hist.GetZaxis().SetMaxDigits(3)

def scan2D_ptG1G2(name, df_base, outdir, fout):
    x0 = pt_min_vals[0] - 0.005
    x1 = pt_min_vals[-1] + 0.005
    nb = len(pt_min_vals)
    hsig = ROOT.TH2D(f"hSig_{name}", "", nb, x0, x1, nb, x0, x1)
    hbkg = ROOT.TH2D(f"hBkg_{name}", "", nb, x0, x1, nb, x0, x1)
    hZ = ROOT.TH2D(f"hZ_{name}", "", nb, x0, x1, nb, x0, x1)
    Zmax = -1.0
    best_pt1 = None
    best_pt2 = None
    # 2D scan
    counter = 0
    for i, pt1 in enumerate(pt_min_vals, start=1):
        print(f"\n>>> Scanning ptG1_min = {pt1:.2f} GeV")
        for j, pt2 in enumerate(pt_min_vals, start=1):
            if pt2 > pt1:
                break
            print(f">>> Scanning ptG2_min = {pt2:.2f} GeV")
            tag = f"pt1_{pt1:.2f}_pt2_{pt2:.2f}".replace(".", "p")
            Nsig, Nbkg, Z = compute_Z_for_ptcuts(name, df_base, REGIONS[name], pt1, pt2, tag, counter, outdir, fout)
            hsig.SetBinContent(i, j, Nsig)
            hbkg.SetBinContent(i, j, Nbkg)
            hZ.SetBinContent(i, j, Z)
            if Z > Zmax:
                Zmax = Z
                best_pt1 = float(pt1)
                best_pt2 = float(pt2)
            counter += 1
    c = ROOT.TCanvas(f"cZ_{name}", "", 1200, 1000)
    c.SetRightMargin(0.16)
    c.SetLeftMargin(0.1)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.1)
    # signal
    hist_scan_options(hsig, "pt_{#gamma_{1}} [GeV]", "pt_{#gamma_{2}} [GeV]", "N_{sig}")
    hsig.Draw("COLZ")
    c.SaveAs(outdir + f"Nsig_scan_{name}_ptG1G2.pdf")
    c.SaveAs(outdir + f"Nsig_scan_{name}_ptG1G2.png")
    fout.cd()
    hsig.Write(hsig.GetName(), ROOT.TObject.kOverwrite)
    # background
    hist_scan_options(hbkg, "pt_{#gamma_{1}} [GeV]", "pt_{#gamma_{2}} [GeV]", "N_{bkg}")
    hbkg.Draw("COLZ")
    c.SaveAs(outdir + f"Nbkg_scan_{name}_ptG1G2.pdf")
    c.SaveAs(outdir + f"Nbkg_scan_{name}_ptG1G2.png")
    fout.cd()
    hbkg.Write(hbkg.GetName(), ROOT.TObject.kOverwrite)
    # significance
    # hist_scan_options(hZ, "pt_{#gamma_{1}} [GeV]", "pt_{#gamma_{2}} [GeV]", "Z = #frac{N_{sig}}{#sqrt{N_{sig} + N_{bkg}}}")
    hist_scan_options(hZ, "pt_{#gamma_{1}} [GeV]", "pt_{#gamma_{2}} [GeV]", "Z")
    hZ.Draw("COLZ")
    c.SaveAs(outdir + f"Z_scan_{name}_ptG1G2.pdf")
    c.SaveAs(outdir + f"Z_scan_{name}_ptG1G2.png")
    fout.cd()
    hZ.Write(hZ.GetName(), ROOT.TObject.kOverwrite)
    print(f"\n>>> {name} scan done")
    print(f"    best ptG1_min = {best_pt1:.2f} GeV")
    print(f"    best ptG2_min = {best_pt2:.2f} GeV")
    print(f"    Z_max         = {Zmax:.2f}")
    return best_pt1, best_pt2, Zmax

def compute_Z_for_s4s9cuts(name, df_base, base_cfg, s4s9_1, s4s9_2, tag, counter, outdir, fout):
    cfg = dict(base_cfg)
    cfg["S4S9_1_min"] = float(s4s9_1)
    cfg["S4S9_2_min"] = float(s4s9_2)
    mask = make_mask(cfg)
    df_tmp = df_base.Define(f"mPi0_scan", f"mPi0_cor[{mask}]")
    h = fill_hist(df_tmp, f"mPi0_scan", f"h_scan_{name}_{tag}", htitle="")
    if h.GetEntries() < 50:
        return 1e-6, 1e-6, 1e-6
    f_pol2_scan = ROOT.TF1(f"f_pol2_{name}_scan_{tag}", "gaus(0)+pol2(3)", fit_min, fit_max)
    res_pol2_scan = pol2_fit(h, f_pol2_scan)
    status = int(res_pol2_scan.Status())
    if status != 0:
        return 1e-6, 1e-6, 1e-6
    if counter % 4 == 0:
        save_hist_postfit(h, f_pol2_scan, outdir, fout, f"mpi0_{name}_postfit_{tag}")
    mean, sigma, x1, x2, Nsig, Nbkg, Z = extract_fit_results(f"{name}_scan_{tag}", f_pol2_scan, nsigma, h.GetBinWidth(1))
    return Nsig, Nbkg, Z

def scan2D_S4S9(name, df_base, outdir, fout):
    if name == "EB_low":
        x0 = S4S9_low_min_vals[0] - 0.005
        x1 = S4S9_low_min_vals[-1] + 0.005
        nb = len(S4S9_low_min_vals)
    else:
        x0 = S4S9_high_min_vals[0] - 0.005
        x1 = S4S9_high_min_vals[-1] + 0.005
        nb = len(S4S9_high_min_vals)
    hsig = ROOT.TH2D(f"hSig_{name}", "", nb, x0, x1, nb, x0, x1)
    hbkg = ROOT.TH2D(f"hBkg_{name}", "", nb, x0, x1, nb, x0, x1)
    hZ = ROOT.TH2D(f"hZ_{name}", "", nb, x0, x1, nb, x0, x1)
    Zmax = -1.0
    best_S4S9_1 = None
    best_S4S9_2 = None
    # 2D scan
    if name == "EB_low":
        scan_vals = S4S9_low_min_vals
    else:
        scan_vals = S4S9_high_min_vals
    counter = 0
    for i, s4s9_1 in enumerate(scan_vals, start=1):
        print(f"\n>>> Scanning S4S9_1_min = {s4s9_1:.3f}")
        for j, s4s9_2 in enumerate(scan_vals, start=1):
            print(f">>> Scanning S4S9_2_min = {s4s9_2:.3f}")
            tag = f"S4S9_1_{s4s9_1:.3f}_S4S9_2_{s4s9_2:.3f}".replace(".", "p")
            Nsig, Nbkg, Z = compute_Z_for_s4s9cuts(name, df_base, REGIONS[name], s4s9_1, s4s9_2, tag, counter, outdir, fout)
            hsig.SetBinContent(i, j, Nsig)
            hbkg.SetBinContent(i, j, Nbkg)
            hZ.SetBinContent(i, j, Z)
            if Z > Zmax:
                Zmax = Z
                best_s4s9_1 = float(s4s9_1)
                best_s4s9_2 = float(s4s9_2)
            counter += 1
    c = ROOT.TCanvas(f"cZ_{name}", "", 1200, 1000)
    c.SetRightMargin(0.16)
    c.SetLeftMargin(0.1)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.1)
    # signal
    hist_scan_options(hsig, "S4S9_{#gamma_{1}}", "S4S9_{#gamma_{2}}", "N_{sig}")
    hsig.Draw("COLZ")
    c.SaveAs(outdir + f"Nsig_scan_{name}_S4S9.pdf")
    c.SaveAs(outdir + f"Nsig_scan_{name}_S4S9.png")
    fout.cd()
    hsig.Write(hsig.GetName(), ROOT.TObject.kOverwrite)
    # background
    hist_scan_options(hbkg, "S4S9_{#gamma_{1}}", "S4S9_{#gamma_{2}}", "N_{bkg}")
    hbkg.Draw("COLZ")
    c.SaveAs(outdir + f"Nbkg_scan_{name}_S4S9.pdf")
    c.SaveAs(outdir + f"Nbkg_scan_{name}_S4S9.png")
    fout.cd()
    hbkg.Write(hbkg.GetName(), ROOT.TObject.kOverwrite)
    # significance
    hist_scan_options(hZ, "S4S9_{#gamma_{1}}", "S4S9_{#gamma_{2}}", "Z")
    hZ.Draw("COLZ")
    c.SaveAs(outdir + f"Z_scan_{name}_S4S9.pdf")
    c.SaveAs(outdir + f"Z_scan_{name}_S4S9.png")
    fout.cd()
    hZ.Write(hZ.GetName(), ROOT.TObject.kOverwrite)
    print(f"\n>>> {name} scan done")
    print(f"    best S4S9_1_min = {best_s4s9_1:.3f}")
    print(f"    best S4S9_2_min = {best_s4s9_2:.3f}")
    print(f"    Z_max         = {Zmax:.2f}")
    return best_s4s9_1, best_s4s9_2, Zmax


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="root file to process")
    parser.add_argument("--output-dir", "-o", required=True, help="output directory plots")
    args = parser.parse_args()
    fname = args.input
    outdir = args.output_dir.rstrip("/") + "/"
    # histo from RDataFrame
    df = ROOT.RDataFrame("Tree_Optim", fname)
    df = (df
        .Define("ptG1_cor", "enG1_cor / cosh(etaG1_cor)")
        .Define("ptG2_cor", "enG2_cor / cosh(etaG2_cor)")
        .Define("abs_etaPi0_cor", "abs(etaPi0_cor)")
    )

    fout = ROOT.TFile(f"{outdir}EB_histograms.root", "RECREATE")
    csv_path = outdir + "EB_best_results.csv"
    csv_file = open(csv_path, "w", newline="")
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(["region", "best_ptG1_min", "best_ptG2_min",
        "mean", "sigma", "x1", "x2",
        "Nsig", "Nbkg", "Z"])
    for name, cfg in REGIONS.items():
        print(f"\n========== REGION {name} (default cuts) ==========")
        mask = make_mask(cfg)
        df_reg = df.Define(f"mPi0_{name}", f"mPi0_cor[{mask}]")
        h = fill_hist(df_reg, f"mPi0_{name}", hname=f"h_mpi0_{name}_prefit_default", htitle="")
        save_hist(h, outdir, fout, plot_name=f"mpi0_{name}_prefit_default")
        print(f"\n========== Fit for {name} (default cuts) ==========")
        f_pol2 = ROOT.TF1(f"f_pol2_{name}", "gaus(0)+pol2(3)", fit_min, fit_max)
        res_pol2 = pol2_fit(h, f_pol2)
        chi2_pol2, ndf_pol2 = compute_chi2_ndf(res_pol2, f_pol2)
        mean_pol2, sigma_pol2, x1_pol2, x2_pol2, Nsig_pol2, Nbkg_pol2, Z_pol2 = extract_fit_results(name, f_pol2, nsigma, h.GetBinWidth(1))
        print_fit_results(chi2_pol2, ndf_pol2, mean_pol2, sigma_pol2, x1_pol2, x2_pol2, Nsig_pol2, Nbkg_pol2, Z_pol2, nsigma)
        save_hist_postfit(h, f_pol2, outdir, fout, plot_name=f"mpi0_{name}_postfit_default")
        print(f"\n========== 2D scan (ptG1, ptG2) for {name} ==========")
        print(f"Scanning ptG1_min in [{pt_min_vals[0]:.2f}, {pt_min_vals[-1]:.2f}] GeV")
        print(f"Scanning ptG2_min in [{pt_min_vals[0]:.2f}, {pt_min_vals[-1]:.2f}] GeV")
        best_pt1, best_pt2, best_Z = scan2D_ptG1G2(name, df, outdir, fout)
        cfg_best = dict(REGIONS[name])
        cfg_best["ptG1_min"] = best_pt1
        cfg_best["ptG2_min"] = best_pt2
        print(f"\n========== {name} BEST POINT postfit ==========")
        print(f"Using ptG1_min={best_pt1:.2f}, ptG2_min={best_pt2:.2f} (Zmax={best_Z:.2f})")
        # print(f"\n========== 2D scan (S4S9_1, S4S9_2) for {name} ==========")
        # best_s4s9_1, best_s4s9_2, best_Z = scan2D_S4S9(name, df, outdir, fout)
        # cfg_best = dict(REGIONS[name])
        # cfg_best["S4S9_1_min"] = best_s4s9_1
        # cfg_best["S4S9_2_min"] = best_s4s9_2
        # print(f"\n========== {name} BEST POINT postfit ==========")
        # print(f"Using s4s9_1_min={best_s4s9_1:.2f}, s4s9_2_min={best_s4s9_2:.2f} (Zmax={best_Z:.2f})")
        mask_best = make_mask(cfg_best)
        df_best = df.Define(f"mPi0_{name}_best", f"mPi0_cor[{mask_best}]")
        h_best = fill_hist(df_best, f"mPi0_{name}_best", hname=f"h_mpi0_{name}_postfit_best", htitle="")
        f_best = ROOT.TF1(f"f_pol2_{name}_best", "gaus(0)+pol2(3)", fit_min, fit_max)
        res_best = pol2_fit(h_best, f_best)
        chi2_best, ndf_best = compute_chi2_ndf(res_best, f_best)
        mean_best, sigma_best, x1_best, x2_best, Nsig_best, Nbkg_best, Z_best = extract_fit_results(f"{name}_best", f_best, nsigma, h_best.GetBinWidth(1))
        print_fit_results(chi2_best, ndf_best, mean_best, sigma_best, x1_best, x2_best, Nsig_best, Nbkg_best, Z_best, nsigma)
        save_hist_postfit(h_best, f_best, outdir, fout, plot_name=f"mpi0_{name}_postfit_best")
        csv_writer.writerow([name, f"{best_pt1:.3f}", f"{best_pt2:.3f}",
            f"{mean_best:.6f}", f"{sigma_best:.6f}", f"{x1_best:.6f}", f"{x2_best:.6f}", 
            f"{Nsig_best:.1f}", f"{Nbkg_best:.1f}", f"{Z_best:.3f}"])
        # csv_writer.writerow([name, f"{best_s4s9_1:.3f}", f"{best_s4s9_2:.3f}",
        #     f"{mean_best:.6f}", f"{sigma_best:.6f}", f"{x1_best:.6f}", f"{x2_best:.6f}", 
        #     f"{Nsig_best:.1f}", f"{Nbkg_best:.1f}", f"{Z_best:.3f}"])
    fout.Close()
    csv_file.close()

if __name__ == "__main__":
    main()
