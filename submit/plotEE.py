import ROOT

# --- open file and get tree ---
f = ROOT.TFile.Open("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run3/marco/AlCaP0_2025G_TestCalib_prova/AlCaP0_2025G_TestCalib/iter_0/AlCaP0_2025G_TestCalib_EcalNtp_8.root")
tree = f.Get("Tree_Optim")

# --- create RDataFrame ---
df = ROOT.RDataFrame(tree)

# --- helper functions ---
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> computePt(
    const ROOT::VecOps::RVec<float>& en,
    const ROOT::VecOps::RVec<float>& eta
){
    ROOT::VecOps::RVec<float> out;
    for(size_t i=0;i<en.size();++i){
        out.push_back(en[i]/cosh(eta[i]));
    }
    return out;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> computePtDiG(
    const ROOT::VecOps::RVec<float>& ptG1,
    const ROOT::VecOps::RVec<float>& ptG2,
    const ROOT::VecOps::RVec<float>& phiG1,
    const ROOT::VecOps::RVec<float>& phiG2
){
    ROOT::VecOps::RVec<float> out;
    for(size_t i=0;i<ptG1.size();++i){
        float px = ptG1[i]*cos(phiG1[i]) + ptG2[i]*cos(phiG2[i]);
        float py = ptG1[i]*sin(phiG1[i]) + ptG2[i]*sin(phiG2[i]);
        out.push_back(sqrt(px*px + py*py));
    }
    return out;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> absVec(const ROOT::VecOps::RVec<float>& vec){
    ROOT::VecOps::RVec<float> out;
    for(auto x: vec) out.push_back(fabs(x));
    return out;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<bool> maskEE(
    const ROOT::VecOps::RVec<float>& absEtaPi0,
    const ROOT::VecOps::RVec<float>& ptG1,
    const ROOT::VecOps::RVec<float>& ptG2,
    const ROOT::VecOps::RVec<float>& ptDiG,
    const ROOT::VecOps::RVec<float>& S4S9_1,
    const ROOT::VecOps::RVec<float>& S4S9_2,
    float etaMin, float etaMax,
    float ptGmin1, float ptGmin2,
    float ptDiGmin,
    float s4s9Min1, float s4s9Min2
){
    ROOT::VecOps::RVec<bool> mask;
    for(size_t i=0;i<absEtaPi0.size();++i){
        bool pass = 
            absEtaPi0[i] > etaMin &&
            absEtaPi0[i] < etaMax &&
            ptG1[i] > ptGmin1 &&
            ptG2[i] > ptGmin2 &&
            ptDiG[i] > ptDiGmin &&
            S4S9_1[i] > s4s9Min1 &&
            S4S9_2[i] > s4s9Min2;
        mask.push_back(pass);
    }
    return mask;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<bool> maskEB(
    const ROOT::VecOps::RVec<float>& absEtaPi0,
    const ROOT::VecOps::RVec<float>& ptG1,
    const ROOT::VecOps::RVec<float>& ptG2,
    const ROOT::VecOps::RVec<float>& ptDiG,
    const ROOT::VecOps::RVec<float>& S4S9_1,
    const ROOT::VecOps::RVec<float>& S4S9_2
){
    ROOT::VecOps::RVec<bool> mask;
    for(size_t i=0;i<absEtaPi0.size();++i){
        bool pass = 
            absEtaPi0[i] < 1.479 &&
            ptG1[i] > 0.5 &&
            ptG2[i] > 0.5 &&
            ptDiG[i] > 1.5 &&
            S4S9_1[i] > 0.65 &&
            S4S9_2[i] > 0.65;
        mask.push_back(pass);
    }
    return mask;
}
""")

# --- define photon pt and diphoton pt ---
df = df.Define("ptG1", "computePt(enG1_cor, etaG1_cor)")
df = df.Define("ptG2", "computePt(enG2_cor, etaG2_cor)")
df = df.Define("ptDiG", "computePtDiG(ptG1, ptG2, phiG1_cor, phiG2_cor)")
df = df.Define("absEtaPi0", "absVec(etaPi0_cor)")

# --- define masks ---
# df = df.Define("EE1_mask", "maskEE(absEtaPi0, ptG1, ptG2, ptDiG, S4S9_1, S4S9_2, 1.479, 1.8, 1.1, 1.1, 3.75, 0.85, 0.85)")
df = df.Define("EE1_mask", "maskEE(absEtaPi0, ptG1, ptG2, ptDiG, S4S9_1, S4S9_2, 1.479, 1.8, 1.1, 1.1, 0, 0, 0)")
# df = df.Define("EE2_mask", "maskEE(absEtaPi0, ptG1, ptG2, ptDiG, S4S9_1, S4S9_2, 1.8, 2.0, 0.95, 0.95, 2.0, 0.92, 0.92)")
df = df.Define("EE2_mask", "maskEE(absEtaPi0, ptG1, ptG2, ptDiG, S4S9_1, S4S9_2, 1.8, 2.0, 0, 0, 0, 0, 0)")
df = df.Define("EE3_mask", "maskEE(absEtaPi0, ptG1, ptG2, ptDiG, S4S9_1, S4S9_2, 2.0, 100.0, 0.95, 0.95, 2.0, 0.92, 0.92)")
df = df.Define("EB_mask", "maskEB(absEtaPi0, ptG1, ptG2, ptDiG, S4S9_1, S4S9_2)")

# --- filtered mPi0 ---
def filter_selected(col_name, mask_name):
    return f"""
    ROOT::VecOps::RVec<float> out;
    for (size_t i=0; i<{col_name}.size(); ++i){{
        if ({mask_name}[i]) out.push_back({col_name}[i]);
    }}
    return out;
    """

df = df.Define("mPi0_EE1_sel", filter_selected("mPi0_cor", "EE1_mask"))
df = df.Define("mPi0_EE2_sel", filter_selected("mPi0_cor", "EE2_mask"))
df = df.Define("mPi0_EE3_sel", filter_selected("mPi0_cor", "EE3_mask"))
df = df.Define("mPi0_EB_sel",  filter_selected("mPi0_cor", "EB_mask"))

# --- histograms ---
hist_EE1 = df.Histo1D(("mPi0_EE1", "mPi0 EE1; m [GeV]; Entries", 100, 0, 0.3), "mPi0_EE1_sel")
hist_EE2 = df.Histo1D(("mPi0_EE2", "mPi0 EE2; m [GeV]; Entries", 100, 0, 0.3), "mPi0_EE2_sel")
hist_EE3 = df.Histo1D(("mPi0_EE3", "mPi0 EE3; m [GeV]; Entries", 100, 0, 0.3), "mPi0_EE3_sel")
hist_EB  = df.Histo1D(("mPi0_EB", "mPi0 EB; m [GeV]; Entries", 100, 0, 0.3), "mPi0_EB_sel")

# --- create separate canvases ---
def draw_hist(hist, color, title, filename):
    c = ROOT.TCanvas(title, title, 1200, 1000)
    h_ptr = hist.GetValue()
    h_ptr.SetLineColor(color)
    h_ptr.Draw()
    c.SaveAs(filename)

draw_hist(hist_EE1, ROOT.kRed, "EE1", "/eos/user/d/delvecch/www/pi0_calib/mPi0_EE1.png")
draw_hist(hist_EE2, ROOT.kBlue, "EE2", "/eos/user/d/delvecch/www/pi0_calib/mPi0_EE2.png")
draw_hist(hist_EE3, ROOT.kGreen, "EE3", "/eos/user/d/delvecch/www/pi0_calib/mPi0_EE3.png")
draw_hist(hist_EB,  ROOT.kMagenta, "EB", "/eos/user/d/delvecch/www/pi0_calib/mPi0_EB.png")
