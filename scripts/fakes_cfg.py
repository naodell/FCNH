#! /usr/bin/env python
import subprocess
import shlex
import ROOT as r
from PlotProducer import *
from TableMaker import *

### This is the config file for manipulating 
### histograms using the PlotProducer class.  

selection   = 'fakes'
period      = '2012'
LUMIDATA    = 20.

samples     = ['WZJets3LNu', 'ttbar', 'ZJets']#, 'WJets', 'ZJets_M-10To50']
eleFakes    = ['ele_v1', 'ele_v2', 'ele_v3', 'ele_v4']
muFakes     = ['mu_v1', 'mu_v2']

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)


r.gROOT.SetBatch()

plotter = PlotProducer(inputFile = 'histos/fakeHistograms.root', savePath = '', scale = LUMIDATA, isAFS = False)
plotter.set_period(period)

### DATASETS ###
### Specify the datasets you wish to stack 
### and overlay accordingly. 

plotter.add_datasets(samples)
plotter._overlayList.extend(['DATA_MUON'])

plotter.get_scale_factors()

### VARIABLES ###
### First specify the directories in which your
### histograms are stored.  If directories are 
### not used enter '' as the only entry.  Then 
### list all of the variable names you wish to 
### plot while giving a key value which is the 
### directory that they are located in as a key.

plotter._directoryList1D            = ['fakes']
#plotter._directoryList2D            = ['2D_Dilepton']

#plotter._variableDict['fakes']   = ['leptonMult', 'jetMult', 'bJetMult']
plotter._variableDict['fakes']   = ['DileptonMass']


 ###################   
 ### MAKE PLOTS! ###  
 ###################   

r.gROOT.SetStyle('Plain')
r.gStyle.SetOptStat(0)

plotter.set_input_file('histos/fakeHistograms.root')
plotter.set_save_path('../plots/fakes/' + currentDate)

for category in eleFakes:
    print '\t' + category
    plotter._category = category
    plotter.make_overlays_1D(logScale = False, doRatio = False, doEff = False)


