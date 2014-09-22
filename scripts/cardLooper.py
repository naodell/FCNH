#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle,re
from array import array
import ROOT as r


if __name__ == '__main__':

    inputDir = 'data/fcnh/dataCards/'
    dataCards = os.listdir(inputDir)

    h_limits = r.TH1D('h1_limits', 'LIMITS!!!;BR(t#rightarrow hc);Entries', 50, 0.7, 1.3)

    bestLimit = [10., '']
    limitDict = {}
    for dataCard in dataCards[:]:

        ### parse datacard name:
        ### format is <variable>_<ybin>_<xbin> 
        parsedName = re.split('_|\.', str(dataCard))
        variable = parsedName[0]
        yBin = parsedName[1]
        xBin = parsedName[2]

        if yBin not in limitDict.keys(): limitDict[yBin] = [10., '']
        
        #print 'combine -M Asymptotic {0}{1}'. format(inputDir, dataCard)
        os.system('combine -M Asymptotic {0}{1} > output_tmp'. format(inputDir, dataCard))
        limitFile = open('output_tmp', 'r')

        for line in limitFile:
            if line.find('Expected 50.0%') is not -1:
                expLimit = float(line.split()[-1])
                h_limits.Fill(expLimit)

                if expLimit < limitDict[yBin][0]:
                    print 'Expected limit at 50 % for MET = {0} and HT = {1}: {2}'.format(xBin, yBin, expLimit)
                    limitDict[yBin] = [expLimit, xBin]


    dataFile = open('MetHT_mask.txt', 'w')
    dataFile.write('#yBin, xRangeLow, xRangeHigh, expLimit\n')
    for key in sorted(limitDict.keys()):
        expLimit, xRange = limitDict[key][0], limitDict[key][1]
        xRangeLow, xRangeHigh = xRange.split('-')
        dataFile.write('{0} {1} {2} {3} \n'.format(key, xRangeLow, xRangeHigh, expLimit))

    dataFile.close()
