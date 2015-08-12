#! /usr/bin/env python
import ROOT as r

if __name__ == '__main__':
    inFile_nominal = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20150222_234009.root')

### MUONS ###




    muFile_up   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20150223_122316.root')
    muFile_down = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20150223_122342.root')

    categories      = ['ss_inclusive', 'ss_mumu', 'ss_emu', 'ss_ee', '3l_inclusive', '3l_mumumu', '3l_emumu', '3l_eemu', '3l_eee']
    signalSamples   = ['FCNC_M125_t'] # signal
    bgSamples       = ['WZJets3LNu'] # Irreducible backgrounds

    ### Do background MC ###
    for sample in bgSamples + signalSamples:
           
        for category in categories:
            yieldsNominal   = inFile_nominal.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
            muScaleUp       = muFile_up.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
            muScaleDown     = muFile_down.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))

            if yieldsNominal.GetBinContent(9) == 0: continue

            systematic = (
                        1-muScaleUp.GetBinContent(9)/yieldsNominal.GetBinContent(9),
                        1-muScaleDown.GetBinContent(9)/yieldsNominal.GetBinContent(9),
                        )


            print '{0} :: {1} :: {2[0]:.3f} -- {2[1]:.3f} '.format(sample, category, systematic)
        
