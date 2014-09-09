import sys, os, shutil, pickle, datetime
from math import *
import ROOT as r

paramFile   = open('scripts/fcncParams.pkl', 'rb')
scales      = pickle.load(paramFile)
styles      = pickle.load(paramFile)
combos      = pickle.load(paramFile)
removes     = pickle.load(paramFile)
categories  = pickle.load(paramFile)
systematics = pickle.load(paramFile)


def bin_by_bin_hist_addition(hist1, hist2, constant = -1.):
    '''
    Custom histogram adder that does not allow negative bins
    '''

    for bin in range(hist1.GetNbinsX()):
        entries1, entries2  = hist1.GetBinContent(bin+1), hist2.GetBinContent(bin+1)
        error1, error2      = hist1.GetBinError(bin+1), hist2.GetBinError(bin+1)

        if hist1 > hist2:
            hist1.SetBinContent(bin+1, entries1 + constant*entries2)
            hist1.SetBinError(bin+1, sqrt(pow(error1, 2) + pow(constant*entries2, 2)))
        else:
            hist1.SetBinContent(bin+1, 0)
            hist1.SetBinError(bin+1, 0)

    return hist1

class AnalysisTools():
    '''
    Base class for analysis processes
    '''
    def __init__(self, inputFile, scale, savePath = ''):
        self._histFile      = r.TFile(inputFile, 'OPEN')
        #self._outFile       = r.TFile('test.root', 'OPEN')
        self._period        = '2012'
        self._savePath      = savePath
        self._rebinFactor   = 1
        self._scale         = scale
        self._scaleDict     = scales
        self._styleDict     = styles
        self._combineDict   = combos
        self._cleanDict     = removes
        self._cleanFakes    = True
        self._category      = ''
        self._datasets      = []

    def set_input_file(self, fileName):
       self._histFile = r.TFile(fileName, 'OPEN') 

    def set_save_path(self, savePath):
        self._savePath = savePath

    def set_category(self, category):
        self._category = category

    def set_period(self, period):
        self._period = period

    def set_clean_fakes(self, doClean):
        self._cleanFakes = doClean

    def set_rebin_factor(self, factor):
        self._rebinFactor = factor

    def get_category(self):
        return self._category

    def get_current_time(self):
        ''' 
        Returns a string of the current time with
        the format  
        '''

        now = datetime.datetime.now()
        currentTime = '{0:02d}{1:02d}{2:02d}_{3:02d}{4:02d}{5:02d}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
        return currentTime


    def add_datasets(self, newData, Clear=False):
        '''
        Add new data
        '''

        if Clear == True:
            self._datasets = []

        if type(newData) == list:
            self._datasets.extend(newData)
        elif type(newData) == str:
            self._datasets.append(newData)


    def make_save_path(self, filePath, clean = False):
        '''
        Create save path in case it doesn't already exist
        '''
        if not os.path.exists(filePath):
            os.system('mkdir -p '+filePath)
        elif len(os.listdir(filePath)) is not 0 and clean:
            os.system('rm -r {0}'.format(filePath))


    def binomial_error(self, efficiency, N0):
        '''
        calculates the error of a ratio
        '''

        err = sqrt(efficiency*(1 - efficiency)/N0)
        return err


    def get_scale_factors(self, addData = [], corrected = True):
        '''
        Gets the initial number of events for each
        dataset being run over
        '''

        for dataName in (self._datasets + addData):
            if dataName[:4] == 'DATA': continue

            if dataName in self._combineDict:
                print dataName, ':\t',
                for data in self._combineDict[dataName]:
                    print data,

                    if not self._histFile.GetDirectory('inclusive/' + data):
                        print '\nCould not find {0} in root file!'.format(data)
                        continue

                    yieldHist = self._histFile.GetDirectory('inclusive/' + data).Get('h1_YieldByCut')
                    nInit       = yieldHist.GetBinContent(1)
                    if corrected:
                        nRaw        = self._histFile.GetDirectory('inclusive/' + data).Get('h1_YieldByCutRaw').GetBinContent(6)
                        nWeighted   = yieldHist.GetBinContent(6)

                        nInit = nInit - (nRaw - nWeighted)

                    self._scaleDict[self._period][data] = 1e3*self._scaleDict[self._period][data]/nInit 

            else:
                if dataName in self._scaleDict[self._period]:
                    print dataName,

                    if not self._histFile.GetDirectory('inclusive/' + dataName):
                        print '\nCould not find {0} in root file!'.format(dataName)
                        continue

                    yieldHist = self._histFile.GetDirectory('inclusive/' + dataName).Get('h1_YieldByCut')
                    nInit       = yieldHist.GetBinContent(1)

                    if corrected:
                        nRaw        = self._histFile.GetDirectory('inclusive/' + dataName).Get('h1_YieldByCutRaw').GetBinContent(6)
                        nWeighted   = yieldHist.GetBinContent(6)

                        nInit = nInit - (nRaw - nWeighted)


                    self._scaleDict[self._period][dataName] = 1e3*self._scaleDict[self._period][dataName]/nInit 

                else:
                    print '{0} not found in scale dictionary; setting to 0'.format(dataName)
                    self._scaleDict[self._period][dataName] = 0.
                    continue

            print ''
        print '\n'


    def get_hist(self, var, dataName, histType = '1D', doScale = True):
        '''
        Get histogram from ROOT file
        '''
        histogramName = ''
        if histType == '1D':
            histogramName = 'h1_' + var  
        if histType == '2D':
            histogramName = 'h2_' + var  

        if not self._histFile.GetDirectory(self._category + '/' + dataName):
            #print self._category, dataName, histogramName
            return None

        # Get clone of histogram stored in root file
        inHist = self._histFile.GetDirectory(self._category + '/' + dataName).Get(histogramName)

        if not inHist:
            return None
        else:
            hist = inHist.Clone()

        if self._category in systematics:
            if dataName == 'QFlips':
                hist = self.add_systematic(hist, 'QFlips')
            elif dataName in self._combineDict['Fakes']:
                hist = self.add_systematic(hist, dataName)
            elif dataName in self._combineDict['Irreducible']:
                hist = self.add_systematic(hist, 'Irreducible')

        if dataName.split('_')[0] in ['Fakes', 'eFakes', 'muFakes', 'llFakes'] and len(dataName.split('_')) > 1:
            dataName = dataName.split('_', 1)[1]

        if doScale:
            if dataName[:4] == 'DATA' or dataName in ['Fakes', 'eFakes', 'muFakes', 'llFakes', 'QFlips', 'AIC', 'eeeAIC', 'eemuAIC', 'emumuAIC', 'mumumuAIC']:
                if self._category == 'ss_ee' and dataName == 'eFakes':
                    hist.Scale(self._scaleDict[self._period][dataName])
                else:
                    hist.Scale(self._scaleDict[self._period][dataName])

            else:
                hist.Scale(self._scale*self._scaleDict[self._period][dataName]) 

        
        #if histType == '1D' and hist.GetNbinsX()%self._rebinFactor is 0:
        #    hist.Rebin(self._rebinFactor)

        return hist

    def combine_samples(self, var, dataName, histType = '1D'):
        '''
        Combines histograms from different samples into one histogram
        '''

        outHist = None
        doFakes = False
        fakeCats = ['eFakes', 'muFakes', 'llFakes', 'Fakes']

        if dataName.split('_')[0] in fakeCats: # Treat fakes separately 
            if dataName is 'Fakes': # fakes from data
                for data in self._combineDict['Fakes']:
                    hist = self.get_hist(var, data, histType)

                    if self._cleanFakes and hist:
                        for mc in self._combineDict['Remove_{0}'.format(self._category.split('_')[0])]:
                            mc_hist = self.get_hist(var, data + '_' + mc, histType)
                            if mc_hist is None:
                                continue
                            else:
                                hist.Add(mc_hist, -1.)
                                #bin_by_bin_hist_addition(outHist, mc_hist, -1.)

                    if not outHist:
                        outHist = hist
                    elif hist:
                        outHist.Add(hist, 1)

            elif dataName in fakeCats: # fakes from data
                outHist = self.get_hist(var, dataName, histType)
                if self._cleanFakes and outHist:
                    for mc in self._combineDict['Remove_{0}'.format(self._category.split('_')[0])]:
                        mc_hist = self.get_hist(var, dataName + '_' + mc, histType)

                        if not mc_hist:
                            continue
                        elif mc_hist:
                            outHist.Add(mc_hist, -1.)
                            #bin_by_bin_hist_addition(outHist, mc_hist, -1.)

            else: # MC fakes
                if dataName.split('_')[1] not in self._combineDict:
                    outHist = self.get_hist(var, dataName, histType)
                else:
                    for data in self._combineDict[dataName.split('_')[1]]:
                        hist = self.get_hist(var, dataName.split('_')[0] + '_' + data, histType)

                        if not outHist:
                            outHist = hist
                        elif hist:
                            outHist.Add(hist, 1)

        else: # Non-fake datasets
            if dataName not in self._combineDict:
                outHist = self.get_hist(var, dataName, histType)
            else:
                for data in self._combineDict[dataName]:
                    hist = self.get_hist(var, data, histType)

                    if outHist is None and hist is not None:
                        outHist = hist
                    elif hist is not None:
                        outHist.Add(hist, 1)

            if dataName in self._cleanDict:
                for data in self._cleanDict[dataName]:
                    hist = self.get_hist(var, data, histType)

                    if outHist is None and hist is not None:
                        outHist = hist
                    elif hist is not None:
                        outHist.Add(hist, -1.)

        return outHist


    def add_systematic(self, hist, dataset):
        '''
        Add a systematic uncertainty to a histogram based the data sample and the
        category.
        '''

        if dataset not in systematics[self._category]:
            return hist

        for bin in range(hist.GetNbinsX()):

            entries = hist.GetBinContent(bin+1)
            if entries == 0:
                continue

            errorSq = pow(hist.GetBinError(bin+1), 2)
            for syst in systematics[self._category][dataset]:
                errorSq += pow(syst*entries, 2)

            error = sqrt(errorSq)
            #print self._category, dataset, entries, error, hist.GetBinError(bin+1)
            hist.SetBinError(bin+1, error)

        return hist
