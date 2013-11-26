============== Analysis_CMS ============

Tools for analyses of data from the CMS detector.  The tools here are designed
for a CMS based analysis, but do not depend on the CMSSW framework.  Rather it
    relies mostly on having CERN ROOT.  Below is a description of the directory
    contents and how they can be used.

-------------- ./fcncAnalysis --------------

Contains analysis code specifically for flavor changing higgs analysis.  The
analysis code is a ROOT TSelector which is designed for looping over a set of
ROOT ntuples (TTrees in our case).  The main piece of code is
fcncAnalyzer_Template.C which is run with the following command,

./fcncLocal.csh <suffix> <selection> <period>

where <suffix> is to identify the sample you are running over (creates a
directory in the output root file), <selection> specifies the trigger stream
that is desired (choice are electron, muon, mueg, or mc; mc doesn't apply any
trigger selection), and period is either 2011 or 2012 (currently only works for
2012).  This will output a number of histograms to the directory 'histos'.  The
script fcncLocal.csh is used to replace some information in the analysis code
and execute the run.C file.

Additionally, there are some python scripts for submitting jobs to the lpc
batch system.  These allow you to select which datasets to run over and how you
would like to split up the jobs for each data sample.  The executable is
batch_cfg.py.  The class for batch submission is BatchMaster.py, and the file
execBatch.csh is needed to configure the batch staging environment.

After the batch_cfg.py file has been modified to run over the desired datasets,
it can be run with the following command

./batch_cfg.py <period> <list of data>

Both of the arguements are based on how the datasets have been defined in the
config.  By default, <period> is assumed to be either 2011 or 2012 and the list
of data is 'data', 'bg', or 'signal'.  

-------------- ./fakeAnalysis --------------

Similar in structure to fcncAnalysis, the fakeAnalysis code is based off of the
ROOT TSelector class.  It is simpler and only requires that one execute the
run.C file in ROOT as follows,

root -l run.C

The output is a set of histograms for determining the probability of leptons
being faked by objects of the users definition.  Definitions of these objects
needs to be done in the relevant sections of the Selector class described
briefly below.

----------- ./interface -----------

Contains header files for a set of container classes for physics objects.  The
base class is TCPhysObject which inherits from the ROOT
TLorentzVector(http://root.cern.ch/root/html/TLorentzVector.html).  This allows
for saving kinematic variables and values that are common to various of physics
    objects (leptons, jets, etc.).  The current set of classes are

- TCJet
- TCMuon
- TCElectron
- TCTau (needs additional development)
- TCMet
- TCTriggerObject
- TCPhoton
- TCGenParticle
- TCGenJet

All classes inherit from TCPhysObject with the exception of TCMet which
inherits from TVector2 (http://root.cern.ch/root/html/TVector2.html).  


----- ./src -----

Source files for container classes.

--------- ./plugins ---------

Contains definitions of classes that are helpful in general for physics
classes.  This includes both header and source files of the following classes,

- HistManager: Allows for filling and booking of histograms in a single line.
  Thanks, Andy.
- Selector: For selecting physics objects (electrons, muons, jets, etc.).
  Carries out identification of leptons based on a set of working points and
  produces a set of selected objects which can be retrieved from the
  instantiation of the Selector class.
- TriggerSelector: Used for identifying whether a given event has passed the
  trigger requirement.  The trigger status is stored in 64-bit unsigned integer
  so there's some extra unpacking to figure out if a given trigger (specified
  by a string) has passed :|
- WeighUtils: Used for applying weights to a supplied event.  This is done
  based upon various event variables, e.g.the kinematics of the leptons.
- EGammaMvaEleEstimator: a tool for applying BDT based identification of
  electrons.  This is passed to the selector class.

--------- ./scripts ---------

A variety of scripts for doing such tasks as making overlays of histograms,
making tables, simple fits, producing efficiency plots, etc.

For producing plots comparing background, data, and signal one can make use of
the PlotProducer.py class. The working configuration for this is plot_cfg.py.
This is really most useful if you have many plots with many background
components because it requires a lot of configuration.

For producing tables the TableMaker class is useful. This is written with
producing tables displaying the cut flow for a given analysis, but can be
configured to output the entries of a histogram into nice LaTex or Twiki
formatted tables.

