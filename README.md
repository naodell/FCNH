Analysis_CMS
============

Tools for analyses of data from the CMS detector.  The tools here are designed for a CMS based analysis, but do not depend on the CMSSW framework.  Rather it relies mostly on having CERN ROOT.  Below is a description of the directory contents and how they can be used.

-----------
./interface
-----------

Contains header files for a set of container classes for physics objects.  The base class is TCPhysObject which inherits from the ROOT TLorentzVector(http://root.cern.ch/root/html/TLorentzVector.html).  This allows for saving kinematic variables and values that are common to various of physics objects (leptons, jets, etc.).  The current set of classes are

- TCJet
- TCMuon
- TCElectron
- TCTau (needs additional development)
- TCMet
- TCTriggerObject
- TCPhoton
- TCGenParticle
- TCGenJet

All classes inherit from TCPhysObject with the exception of TCMet which inherits from TVector2 (http://root.cern.ch/root/html/TVector2.html).  


-----
./src
-----

Source files for container classes.

---------
./plugins
---------

Contains definitions of classes that are helpful in general for physics classes.  This includes both header and source files of the following classes,

- HistManager: Allows for filling and booking of histograms in a single line.  Thanks, Andy.
- Selector: For selecting physics objects (electrons, muons, jets, etc.).  Carries out identification of leptons based on a set of working points and produces a set of selected objects which can be retrieved from the instantiation of the Selector class.
- TriggerSelector: Used for identifying whether a given event has passed the trigger requirement.  The trigger status is stored in 64-bit unsigned integer so there's some extra unpacking to figure out if a given trigger (specified by a string) has passed :|
- WeighUtils: Used for applying weights to a supplied event.  This is done based upon various event variables, e.g.the kinematics of the leptons.
- EGammaMvaEleEstimator: a tool for applying BDT based identification of electrons.  This is passed to the selector class.

---------
./scripts
---------
