#! /usr/bin/env python
import ROOT as r


if __name__ == '__main__':
    inFile_nominal = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20141021_164352.root')

### MUONS ###

    muFile_up   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20141021_164514.root')
    muFile_down = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20141021_164732.root')
    elFile_up   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20141022_113102.root')
    elFile_down = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20141022_113216.root')

    categories      = ['ss_mumu', 'ss_emu', 'ss_ee', '3l_mumumu', '3l_emumu', '3l_eemu', '3l_eee']
    signalSamples   = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal
    bgSamples       = ['WZJets3LNu', 'ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttZ', 'ttW', 'ttG'] # Irreducible backgrounds

    ### Do background MC ###
    for sample in bgSamples + signalSamples:
            
        for category in categories:
            yieldsNominal   = inFile_nominal.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
            muScaleUp   = muFile_up.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
            muScaleDown = muFile_down.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
            elScaleUp   = elFile_up.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
            elScaleDown = elFile_down.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))

            if yieldsNominal.GetBinContent(6) == 0: continue

            systematic = (
                        muScaleUp.GetBinContent(6)/yieldsNominal.GetBinContent(6),
                        muScaleDown.GetBinContent(6)/yieldsNominal.GetBinContent(6),
                        elScaleUp.GetBinContent(6)/yieldsNominal.GetBinContent(6),
                        elScaleDown.GetBinContent(6)/yieldsNominal.GetBinContent(6)
                        )


            print '{0} :: {1} :: {2[0]:.2f} -- {2[1]:.2f} :: {2[2]:.2f} -- {2[3]:.2f}'.format(sample, category, systematic)
        
