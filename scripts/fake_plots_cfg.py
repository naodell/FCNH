#! /usr/bin/env python
import subprocess
import shlex
import ROOT as r
from PlotProducer import *
from TableMaker import *

### This is the config file for manipulating 
### histograms using the PlotProducer class.  

selection   = 'fakes'
cutList     = [
                '1_preselection'
                #'2_MET', 
                #'3_Z_veto', 
                #'4_top_veto', 
                #'5_dilepton_dR'

                #'2_Z_veto', 
                #'3_jet_cut', 
                #'4_soft_lepton_veto', 
                #'5_MET'
              ]

period      = '2012'
LUMIDATA    = 19.1

samples     = ['WZJets3LNu', 'ttbar', 'ZJets']#, 'WJets', 'ZJets_M-10To50']
cat3l       = ['3l_inclusive', '3l_OSSF', '3l_SSSF']

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)


r.gROOT.SetBatch()

plotter = PlotProducer(inputFile = 'histos/' + selection + '_cut1_' + period + '.root', savePath = '', scale = LUMIDATA, isAFS = False)
plotter.set_period(period)

### DATASETS ###
### Specify the datasets you wish to stack 
### and overlay accordingly. 

plotter.add_datasets(samples)
plotter._overlayList.extend(['fakes'])

plotter.get_scale_factors()

### VARIABLES ###
### First specify the directories in which your
### histograms are stored.  If directories are 
### not used enter '' as the only entry.  Then 
### list all of the variable names you wish to 
### plot while giving a key value which is the 
### directory that they are located in as a key.

plotter._directoryList1D            = ['Misc', 'Lepton', 'Dilepton', 'MET', 'Jet', 'GEN']
#plotter._directoryList2D            = ['2D_Dilepton']

plotter._variableDict['Misc']       = ['PvMult', 'YieldByCut', 'EventWeight', 'TriggerStatus']

plotter._variableDict['Lepton']     = ['LeptonCharge', 'LeptonFlavor', 
                                       'Lepton1Pt', 'Lepton2Pt','Lepton3Pt',
                                       'Lepton1Eta', 'Lepton2Eta', 'Lepton3Eta' ]
                                       #'Lepton1Phi', 'Lepton2Phi', 'Lepton3Phi']

plotter._variableDict['Dilepton']   = ['DileptonMass21', 'DileptonTransMass21', 'DileptonQt21',
                                       'DileptonDeltaPhi21', 'DileptonDeltaEta21', 'DileptonDeltaR21', 'DileptonDeltaPt21',
                                       'DileptonMass31', 'DileptonTransMass31', 'DileptonQt31', 
                                       'DileptonDeltaPhi31', 'DileptonDeltaEta31', 'DileptonDeltaR31', 'DileptonDeltaPt31',
                                       'DileptonMass32', 'DileptonTransMass32', 'DileptonQt32', 
                                       'DileptonDeltaPhi32', 'DileptonDeltaEta32', 'DileptonDeltaR32', 'DileptonDeltaPt32',
                                       'DileptonOSMass', 'DileptonOSTransMass']

plotter._variableDict['top']        = ['Top1Mass', 'Top1TransMass', 'Top1Pt', 'DeltaPhiTop1Met',
                                       'Top2TransMass', 'Top2Pt', 'DeltaPhiTop2Met', 'DeltaPhiTop1Top2Met']


plotter._variableDict['Lep+Jet']    = ['Lepton1BJetDeltaPhi', 'Lepton1BJetDeltaEta', 'Lepton1BJetDeltaR', 'Lepton1BJetDeltaPt',
                                       'Lepton2BJetDeltaPhi', 'Lepton2BJetDeltaEta', 'Lepton2BJetDeltaR', 'Lepton2BJetDeltaPt',
                                       'Lepton3BJetDeltaPhi', 'Lepton3BJetDeltaEta', 'Lepton3BJetDeltaR', 'Lepton3BJetDeltaPt'
                                       ]


plotter._variableDict['Jet']        = ['Jet1Pt', 'Jet2Pt', 'Jet3Pt',
                                       'Jet1Eta', 'Jet2Eta', 'Jet3Eta',
                                       #'Jet1Phi', 'Jet2Phi', 'Jet3Phi',
                                       'BJet1Discr', 'BJet1Pt', 'BJet1Eta', 'BJet1Phi', 
                                       'BJet2Discr', 'BJet2Pt', 'BJet2Eta', 'BJet2Phi',
                                       'HTs', 'EventBalance', 'JetMultCharge', 
                                       'JetMult', 'BJetMult']

plotter._variableDict['MET']        = ['Met', 'MetPhi', 'MetSumEt',
                                       'MetLepton1DeltaPhi', 'MetLepton2DeltaPhi'
                                       'MetLepDeltaPhiMin', 'nearLepIndex', 'ProjectedMet', 'MetLepton3DeltaPhi'] 

plotter._variableDict['GEN']        = ['GenChargeMisId', 'GenMisIdPt', 'GenMisIdEta',
                                       'GenDeltaR', 'GenBalance']

plotter._variableDict['2D_Dilepton']    = ['DileptonM12VsLepPt1', 'DileptonM12VsLepPt2', 'DileptonMVsDeltaPhi12', 'DileptonM12VsQt12',
                                           'DileptonM12VsM23', 'DileptonM12VsM13', 'DileptonM23VsM13', 'DileptonM12VsLepPt3', 
                                           'DileptonM13VsLepPt1', 'DileptonM13VsLepPt2', 'DileptonM13VsLepPt3',
                                           'DileptonM23VsLepPt1', 'DileptonM23VsLepPt2', 'DileptonM23VsLepPt3',
                                           'DileptonMos1VsMss', 'DileptonMos2VsMss']
                                           #'DileptonMVsDeltaPhi13', 'DileptonM12VsQt13', 'DileptonMVsDeltaPhi23', 'DileptonM12VsQt23'


 ###################   
 ### MAKE PLOTS! ###  
 ###################   

r.gROOT.SetStyle('Plain')
r.gStyle.SetOptStat(0)
#r.gROOT.ProcessLine('.L ./tdrStyle.C')
#r.setTDRStyle()


for i, cut in enumerate(cutList):
    print cut

    plotter.set_input_file('histos/' + selection + '_cut' + str(i+1) + '_' + period + '.root')
    plotter.set_save_path('../plots/fakes/' + currentDate + '/' + cut)

    for category in cat3l:
        print '\t' + category
        plotter._category = category
        plotter.make_overlays_1D(logScale = True, doRatio = True, doEff = False)


