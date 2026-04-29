#!/usr/bin/env python3
import ROOT

ROOT.gROOT.SetBatch(True)  # niente GUI, utile su lxplus
ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetMaxDigits(3)

path = "/eos/user/d/delvecch/www/pi0_calib/EE_ptG_0p00_1p10/"
# path = "/eos/user/d/delvecch/www/pi0_calib/EEptG1p1DOWN/"
f = ROOT.TFile.Open(f"{path}EE_histograms.root")
if not f or f.IsZombie():
    raise RuntimeError(f"Impossibile aprire: {infile}")

name = "EE_low"
hist_type = "hZ"
h = f.Get(f"{hist_type}_{name}")
if not h:
    raise RuntimeError(f"Istogramma 'h{hist_type}_{name}' non trovato nel file")

# IMPORTANTISSIMO: stacca l'istogramma dal file, così puoi chiudere il file senza perdere l'oggetto
h.SetDirectory(0)
f.Close()

# --- ZOOM Z: imposta qui i limiti che vuoi ---
zmin = 30
zmax = 110
h.SetMinimum(zmin)
h.SetMaximum(zmax)

c = ROOT.TCanvas("c", "", 1200, 1000)
c.SetRightMargin(0.14)  # spazio per la palette
h.Draw("COLZ")
c.SaveAs(f"{path}{hist_type}_{name}_zoom.png")
c.SaveAs(f"{path}{hist_type}_{name}_zoom.pdf")
