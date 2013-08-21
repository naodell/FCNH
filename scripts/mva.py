#! /usr/bin/env python
import subprocess, shlex, time, pickle, math, sys, os
from array import array
import ROOT as r

### Some useful functions ###

def add_scale_branch(inputFile, sampleList, scales):
    '''
    Merges trees and adds a branch with weights based on cross-section.
    '''

    treeList    = []

    print 'Adding samples to', 

    for sample in sampleList:
        print sample,

        nInit   = inputFile.Get('inclusive/' + sample + '/h1_YieldByCut').GetBinContent(1)
        tree    = inputFile.Get('mvaTree_' + sample)
        scale   = array('f', [1e3*scales['2012'][sample]/nInit]) 
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
        batch = sys.argv[1]
    else:
        print 'Must provide information about input file!'
        exit()


    # Configuration parameters
    methods     = ['BDT']
    doGUI       = True
    selection   = '3l'

    # Scale factors
    paramFile = open('scripts/fcncParams.pkl', 'rb')
    scales    = pickle.load(paramFile)

    # Samples for merging
    bgList  = []
    bgList.extend(['ZJets', 'ZJets_M-10To50', 'WJets']) # V+jets
    bgList.extend(['WWJets2L2Nu', 'ZZJets2L2Nu']) #, 'ZZJets2L2Q', 'WZJets2L2Q']) # Diboson to 2l + X
    #bgList.extend(['WZJets3LNu']) # WZ to 3l+nu
    bgList.extend(['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau']) # ZZ to 4l
    bgList.extend(['tW', 'tbarW', 't_t-channel', 't_t-channel', 'ttbar']) # Top
    bgList.extend(['ttZ', 'ttW', 'ttG']) # Top+V
    #bgList.extend(['QCD_20_MU', 'QCD_20-30_EM', 'QCD_30-80_EM', 'QCD_80-170_EM', 'QCD_170-250_EM', 'QCD_250-350_EM', 'QCD_350_EM']) #QCD
    bgList.extend(['ggHToZZ4L_M-125', 'WHToWWW3L_M-125']) # Higgs

    sigList = ['FCNC_M125_t', 'FCNC_M125_tbar']

    # Input file and tree merging
    inFile  = r.TFile('histos/fcnh_cut7_2012_{0}.root'.format(batch), 'OPEN')

    bgTrees     = add_scale_branch(inFile, bgList, scales)
    sigTrees    = add_scale_branch(inFile, sigList, scales)

    ### Start settting up MVA ###

    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    r.gROOT.SetMacroPath("${ROOTSYS}/tmva/test/.") 
    r.gROOT.Macro       ("${ROOTSYS}/tmva/test/TMVAlogon.C")
    r.gROOT.LoadMacro   ("${ROOTSYS}/tmva/test/TMVAGui.C")

    # Change output directory for weights
    weightDir = 'weights/{0}'.format(batch)
    if not os.path.exists(weightDir):
        os.system('mkdir -p {0}'.format(weightDir))
    elif len(os.listdir(weightDir)) is not 0:
        os.system('rm -r {0}'.format(weightDir))

    r.TMVA.gConfig().GetIONames().fWeightFileDir = weightDir

    # Output file
    outputFile = r.TFile('mvaOutput/{0}.root'.format(batch), 'RECREATE' )

    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = r.TMVA.Factory( "TMVAClassification", outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose(True)

    factory.AddVariable('met', 'met', 'GeV', 'F')
    factory.AddVariable('metPhi', 'metPhi', 'GeV', 'F')
    factory.AddVariable('HT', 'HT', 'GeV', 'F')
    factory.AddVariable('MT', 'MT', 'GeV', 'F')

    factory.AddVariable('lep1Pt', 'lep1Pt', 'GeV', 'F')
    factory.AddVariable('lep1Eta', 'lep1Eta', '', 'F')
    factory.AddVariable('lep1Phi', 'lep1Phi', 'rad', 'F')
    factory.AddVariable('lep2Pt', 'lep2Pt', 'GeV', 'F')
    factory.AddVariable('lep2Eta', 'lep2Eta', '', 'F')
    factory.AddVariable('lep2Phi', 'lep2Phi', 'rad', 'F')

    factory.AddVariable('bJetPt', 'bJetPt', 'GeV', 'F')
    factory.AddVariable('bJetEta', 'bJetEta', '', 'F')
    factory.AddVariable('bJetPhi', 'bJetPhi', 'rad', 'F')

    factory.AddVariable('jetMult', 'jetMult', '', 'I')
    factory.AddVariable('bJetMult', 'bJetMult', '', 'I')

    if selection is '3l':
        factory.AddVariable('lep3Pt', 'lep3Pt', 'GeV', 'F')
        factory.AddVariable('lep3Eta', 'lep3Eta', '', 'F')
        factory.AddVariable('lep3Phi', 'lep3Phi', 'rad', 'F')

        factory.AddVariable('trileptonMass', 'trileptonMass', 'GeV', 'F')
        factory.AddVariable('dileptonMassOS', 'dileptonMassOS', 'GeV', 'F')
        factory.AddVariable('dileptonDROS', 'dileptonDROS', 'rad', 'F')


    for tree in sigTrees:
        factory.AddSignalTree(tree, 1.)

    for tree in bgTrees:
        factory.AddBackgroundTree(tree, 1.)

    factory.SetBackgroundWeightExpression('evtWeight * scale')

    factory.PrepareTrainingAndTestTree(r.TCut(''), r.TCut(''),
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )

    if 'CUTS' in methods:
        #factory.BookMethod( r.TMVA.Types.kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=100000:VarProp=FSmart" )
        factory.BookMethod( r.TMVA.Types.kCuts, "CutsGA",
                            "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" )

    if 'BDT' in methods:
        factory.BookMethod(r.TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=SDivSqrtSPlusB:nCuts=20:PruneMethod=NoPruning" )
        #factory.BookMethod(r.TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=150:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity" )


    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()    

    outputFile.Close()

    if doGUI:
        r.gROOT.ProcessLine('TMVAGui(\"mvaOutput/{0}.root\")'.format(batch))
        r.gApplication.Run() 

