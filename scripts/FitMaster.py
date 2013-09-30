import os, sys, subprocess
import ROOT as r

class FitMaster():
    '''
    This is going to be a general purpose fitting class.  For now,
    it will fit a single 1 dimensional histogram, but will be expanded
    to work for multiple histograms, 2D histograms, etc.  
    '''
    def __init__(self, inputFile, inHist, dataName = 'TEST', path = '.', directory = '', options = 'RM', doLog = False):
        self._inFile   = r.TFile(inputFile, 'OPEN')
        self._hist     = inHist
        self._dir      = directory
        self._path     = path
        self._dataName = dataName
        self._options  = options
        self._doLog    = doLog

    def GetHist(self):
        '''Retreive  histogram from root file'''
        if self._dir == '':
            hist = self._inFile.Get(self._hist)
        else:
            hist = self._inFile.Get(self._dir+'/'+self._hist)
        hist.Print()
        return hist

    def DrawHist(self, hist, option = ''):
        '''For drawing histograms'''
        r.gStyle.SetOptStat(0)
        canvas = r.TCanvas('canvas','canvas', 700, 600)
        if self._doLog:
            canvas.SetLogy()
        canvas.SetGridy()
        canvas.SetGridx()
        hist.Draw()
        canvas.SaveAs(self._path+'/'+self._dataName+'.png')

    def DoFit(self, fitFunction, xLow, xHigh, paramInit = [], draw=True):
        '''Fitting procedure'''
        hist = self.GetHist()    
        fit  = r.TF1('fit', fitFunction, xLow, xHigh);
        for i, init in enumerate(paramInit):
            fit.SetParameter(i, init)

        fitResults = hist.Fit('fit', self._options)
        if draw:
            self.DrawHist(hist)

