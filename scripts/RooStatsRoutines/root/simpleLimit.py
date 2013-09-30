import ROOT as r

### Helper function for calculating binomial error ###
def weighted_binomial_error(nFin, nInit)

    eff = nFin[0]/nInit[0];
    sigma_eff = sqrt(((1 - eff)*nFin[1])**2 + nInit[0]*(1 - eff)*eff**2)/nInit[0];
    return sigma_eff;

######################################################


nCats = 2;
nToys = 10;

catNames = []

catNames.extend(["3l_inclusive", "ss_inclusive"])
    #"3l_OSSF",
    #"3l_SSSF",
    #"3l_mumumu", 
    #"3l_emumu", 
    #"3l_eemu",
    #"3l_eee",
    #"ss_mumu", 
    #"ss_emu", 
    #"ss_ee"

lumi       = 19.7e3;
lumiErr    = 0.04;
sigInit    = 291.6;
sigInitErr = 24.99;

sigList    = [6.2, 7.74]
sSigList   = [0.11, 0.13]
bckList    = [168.95, 1130.14]
sbckList   = [9.49, 35.80]

nObs       = [236, 1846]

r.gROOT.LoadMacro('roostats_cl95.C+')

for cat in catNames:
