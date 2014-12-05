#! /usr/bin/env python
import ROOT as r

def set_graph_style(graph, color, width, marker):
    graph.SetMarkerStyle(marker)
    graph.SetMarkerColor(color)
    graph.SetMarkerSize(0.)
    graph.SetFillColor(0)
    graph.SetLineWidth(width)
    graph.SetLineColor(color)

if __name__ == '__main__':

    canvas = r.TCanvas('canvas', 'canvas', 500, 500)
    legend = r.TLegend(0.6, 0.6, 0.85, 0.85)
    legend.SetFillColor(0)
    legend.SetFillStyle(3001)
    legend.SetLineWidth(0)
    legend.SetLineColor(0)
    legend.SetTextSize(0.03)

    ### Get SF files ###
    triggerFile = r.TFile('data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root', 'OPEN')
    idFile      = r.TFile('data/MuonEfficiencies_Run2012ReReco_53X.root', 'OPEN')
    isoFile     = r.TFile('data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root', 'OPEN');

    ### Trigger SFs ###
    h2_MuTriggerSF1 = triggerFile.Get("DATA_over_MC_Mu17Mu8_OR_Mu17TkMu8_Tight_Mu1_10To20_&_Mu2_20ToInfty_with_STAT_uncrt")
    h2_MuTriggerSF2 = triggerFile.Get("DATA_over_MC_Mu17Mu8_OR_Mu17TkMu8_Tight_Mu1_20ToInfty_&_Mu2_20ToInfty_with_STAT_uncrt")

    h2_MuTriggerSF1.Draw('colz text')
    h2_MuTriggerSF1.SetTitle('Mu17Mu8 OR Mu17TkMu8 scale factors;|#eta_{1}|;|#eta_{2}|')
    canvas.Print('plots/muon_efficiency_official/Mu17Mu8_OR_Mu17TkMu8_Mu1_10To20_&_Mu2_20ToInfty.pdf')
    h2_MuTriggerSF2.Draw('colz text')
    h2_MuTriggerSF2.SetTitle('Mu17Mu8 OR Mu17TkMu8 scale factors;|#eta_{1}|;|#eta_{2}|')
    canvas.Print('plots/muon_efficiency_official/Mu17Mu8_OR_Mu17TkMu8_Mu1_20ToInfty_&_Mu2_20ToInfty.pdf')

    muSF2012_ID_0   = idFile.Get("DATA_over_MC_Tight_pt_abseta<0.9");
    muSF2012_ID_1   = idFile.Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
    muSF2012_ID_2   = idFile.Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");
    muSF2012_ID_3   = idFile.Get("DATA_over_MC_Tight_pt_abseta2.1-2.4");

    muSF2012_ID_0.SetTitle('muon ID scale factors;p_{T,#mu};#varepsilon')
    muSF2012_ID_0.SetMinimum(0.8)
    muSF2012_ID_0.SetMaximum(1.2)

    set_graph_style(muSF2012_ID_0, r.kBlue, 2, 22)
    set_graph_style(muSF2012_ID_1, r.kRed, 2, 22)
    set_graph_style(muSF2012_ID_2, r.kCyan, 2, 22)
    set_graph_style(muSF2012_ID_3, r.kMagenta, 2, 22)

    legend.AddEntry(muSF2012_ID_0, '|#eta| < 0.9')
    legend.AddEntry(muSF2012_ID_1, '0.9 < |#eta| < 1.2')
    legend.AddEntry(muSF2012_ID_2, '1.2 < |#eta| < 2.1')
    legend.AddEntry(muSF2012_ID_3, '2.1 < |#eta| < 2.4')

    muSF2012_ID_0.Draw('AL')
    muSF2012_ID_1.Draw('L SAME')
    muSF2012_ID_2.Draw('L SAME')
    muSF2012_ID_3.Draw('L SAME')
    legend.Draw()
    canvas.Print('plots/muon_efficiency_official/ID_Tight_pt.pdf')

    muSF2012_ISO_0  = isoFile.Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9");
    muSF2012_ISO_1  = isoFile.Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2");
    muSF2012_ISO_2  = isoFile.Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1");
    muSF2012_ISO_3  = isoFile.Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4");

    muSF2012_ISO_0.SetTitle('muon ISO scale factors;p_{T,#mu};#varepsilon')
    muSF2012_ISO_0.SetMinimum(0.9)
    muSF2012_ISO_0.SetMaximum(1.15)

    set_graph_style(muSF2012_ISO_0, r.kBlue, 2, 22)
    set_graph_style(muSF2012_ISO_1, r.kRed, 2, 22)
    set_graph_style(muSF2012_ISO_2, r.kCyan, 2, 22)
    set_graph_style(muSF2012_ISO_3, r.kMagenta, 2, 22)

    legend.Clear()
    legend.AddEntry(muSF2012_ISO_0, '|#eta| < 0.9')
    legend.AddEntry(muSF2012_ISO_1, '0.9 < |#eta| < 1.2')
    legend.AddEntry(muSF2012_ISO_2, '1.2 < |#eta| < 2.1')
    legend.AddEntry(muSF2012_ISO_3, '2.1 < |#eta| < 2.4')

    muSF2012_ISO_0.Draw('AL')
    muSF2012_ISO_1.Draw('L SAME')
    muSF2012_ISO_2.Draw('L SAME')
    muSF2012_ISO_3.Draw('L SAME')
    legend.Draw()
    canvas.Print('plots/muon_efficiency_official/ISO_Tight_pt.pdf')

