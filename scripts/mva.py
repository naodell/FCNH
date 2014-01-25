#! /usr/bin/env python
import subprocess, shlex, time, pickle, math, sys, os
from array import array
import ROOT as r

### Some useful functions ###

def add_scale_branch(inputFile, sampleList, scales, selection):
    '''
    Merges trees and adds a branch with weights based on cross-section.
    '''

    treeList    = []

    print 'Adding samples to', 

    for sample in sampleList:
        print sample,

        tree    = inputFile.Get('tree' + selection + '_' + sample)
        if sample not in ['QFlips', 'Fakes']:
            nInit   = inputFile.Get('inclusive/' + sample + '/h1_YieldByCut').GetBinContent(1)
            scale = array('f', [1e3*scales['2012'][sample]/nInit]) 
        else:
            scale = array('f', [1.])    

        entries = tree.GetEntries()

        if entries == 0:
            continue

        b_scale = tree.Branch('scale', scale, 'scale/F')

        for i in range(entries):
            tree.GetEntry(i)
            b_scale.Fill()

        treeList.append(tree)

    print ''

    return treeList


### Do MVA here ###

if __name__ == '__main__':

    ### Get command line arguements
    if len(sys.argv) > 1:
        batch       = sys.argv[1]
        selection   = sys.argv[2]
        flCat       = sys.argv[3]
    else:
        print 'Must provide information about input file!'
        exit()

    # Configuration parameters
    methods     = ['BDT']
    doGUI       = False
    selections  = ['3l', 'SS']

    flCats = {}
    flCats['3l'] = ['inclusive', 'eee', 'eemu', 'emumu', 'mumumu']
    flCats['SS'] = ['inclusive', 'ee', 'emu', 'mumu']

    flavorCuts = {}
    cut = '(jetMult + bJetMult) > 1'
    if selection == '3l':
        cut += ' && (dileptonMassOS < 76.2 || dileptonMassOS > 106.2)'
    if selection == 'SS':
        cut += ' && (dileptonMass < 76.2 || dileptonMass > 106.2)'

    flavorCuts['3l']  = {'inclusive': 'flavorCat > 0', 'eee':'flavorCat == 5', 'eemu':'flavorCat == 6 || flavorCat == 7 || flavorCat == 9', 'emumu':'flavorCat == 8 || flavorCat == 10 || flavorCat == 11', 'mumumu':'flavorCat == 12'}
    flavorCuts['SS'] = {'inclusive': 'flavorCat > 0', 'ee':'flavorCat == 1', 'emu':'flavorCat == 2 || flavorCat == 3', 'mumu':'flavorCat == 4'}

    # Scale factors
    paramFile = open('scripts/fcncParams.pkl', 'rb')
    scales    = pickle.load(paramFile)

    ### Start settting up MVA ###
    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    r.gROOT.LoadMacro("src/TCPhysObject.cc+");
    r.gROOT.SetMacroPath("${ROOTSYS}/tmva/test/.") 
    r.gROOT.Macro       ("${ROOTSYS}/tmva/test/TMVAlogon.C")
    if doGUI:
        r.gROOT.LoadMacro   ("${ROOTSYS}/tmva/test/TMVAGui.C")

    # Input file and tree merging
    inFile  = r.TFile('histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')

    print'\nCarrying out BDT optimization of {0} selection with category {1}\n'.format(selection, flCat)

    # Change output directory for weights
    weightDir = 'weights/{0}'.format(batch)
    if not os.path.exists(weightDir):
        os.system('mkdir -p {0}'.format(weightDir))
    #elif len(os.listdir(weightDir)) is not 0:
    #    os.system('rm -r {0}'.format(weightDir))

    r.TMVA.gConfig().GetIONames().fWeightFileDir = weightDir

    # Output file
    outputFile = r.TFile('mvaOutput/{0}_{1}_{2}.root'.format(batch, selection, flCat), 'RECREATE' )

    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = r.TMVA.Factory( 'TMVAClassification_{0}_{1}'.format(selection, flCat), outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose(False)

    factory.AddVariable('met', 'met', 'GeV', 'F')
    factory.AddVariable('HT', 'HT', 'GeV', 'F')
    factory.AddVariable('MT', 'MT', 'GeV', 'F')
    #factory.AddVariable('jetMult', 'jetMult', '', 'I')
    factory.AddVariable('bJetMult', 'bJetMult', '', 'I')

    if selection == '3l':
        #factory.AddVariable('trileptonMass', 'trileptonMass', 'GeV', 'F')
        factory.AddVariable('dileptonMassOS', 'dileptonMassOS', 'GeV', 'F')
        #factory.AddVariable('dileptonDROS', 'dileptonDROS', 'rad', 'F')

        factory.AddSpectator( 'flavorCat', 'flavorCat', '', 5., 12.) 
        factory.AddSpectator( 'chargeCat', 'chargeCat', '', 5, 12.)
    elif selection == 'SS':
        factory.AddVariable('dileptonMass', 'dileptonMass', 'GeV', 'F')
        #factory.AddVariable('dileptonDR', 'dileptonDR', 'rad', 'F')

        factory.AddSpectator( 'flavorCat', 'flavorCat', '', 1., 4.) 
        factory.AddSpectator( 'chargeCat', 'chargeCat', '', 1., 4.)

    factory.AddSpectator('jetMult', 'jetMult', '', 0., 15.)
    #factory.AddSpectator('bJetMult', 'bJetMult', '', 0., 5.)
    factory.AddSpectator( 'evtWeight', 'evtWeight', '', 0., 1e9)


    # add background trees
    bgList  = []
    #bgList.extend(['ZJets', 'ZJets_M-10To50', 'WJets']) # V+jets
    #bgList.extend(['tW', 'tbarW', 't_t-channel', 'tbar_t-channel', 'ttbar']) # Top
    #bgList.extend(['QCD_20_MU', 'QCD_20-30_EM', 'QCD_30-80_EM', 'QCD_80-170_EM', 'QCD_170-250_EM', 'QCD_250-350_EM', 'QCD_350_EM']) #QCD
    #bgList.extend(['WWJets2L2Nu', 'ZZJets2L2Nu']) #, 'ZZJets2L2Q', 'WZJets2L2Q']) # Diboson to 2l + X

    bgList.extend(['Fakes']) #FAKES
    bgList.extend(['WZJets3LNu']) # WZ to 3l+nu

    if selection == 'SS':
        bgList.extend(['QFlips']) # Charge flips
    elif selection == '3l':
        bgList.extend(['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau']) # ZZ to 4l
        bgList.extend(['ttZ', 'ttW', 'ttG']) # Top+V
        #bgList.extend(['ggHToZZ4L_M-125', 'WHToWWW3L_M-125']) # Higgs

    bgTrees = add_scale_branch(inFile, bgList, scales, selection)
    for tree in bgTrees:
        factory.AddBackgroundTree(tree, 1.)

    # add signal trees 
    sigList = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_M125_t_semilep', 'FCNC_M125_t_ZZ', 'FCNC_M125_t_TauTau']
    sigTrees = add_scale_branch(inFile, sigList, scales, selection)
    for tree in sigTrees:
        factory.AddSignalTree(tree, 1.)

    factory.SetBackgroundWeightExpression('evtWeight * scale')
    factory.PrepareTrainingAndTestTree(r.TCut('({0}) && ({1})'.format(cut, flavorCuts[selection][flCat])), r.TCut('({0}) && ({1})'.format(cut, flavorCuts[selection][flCat])), "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )

    if 'BDT' in methods:
        factory.BookMethod( r.TMVA.Types.kBDT, "BDTG",
                            '!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=2000:NNodesMax=8:IgnoreNegWeights:nEventsMin=80')
                            #'!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=SDivSqrtSPlusB:nCuts=20:MaxDepth=3' )

        #factory.BookMethod(r.TMVA.Types.kBDT, "BDT", 
        #                    "!H:!V:NTrees=1000:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=SDivSqrtSPlusB:nCuts=20:PruneMethod=NoPruning" )
        #factory.BookMethod(r.TMVA.Types.kBDT, "BDT", 
        #                    "!H:!V:NTrees=500:MinNodeSize=10:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=-1:PruneMethod=NoPruning" )


    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()    

    outputFile.Close()

    if doGUI:
        r.gROOT.ProcessLine('TMVAGui(\"mvaOutput/{0}.root\")'.format(batch))
        r.gApplication.Run() 

