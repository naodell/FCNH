from AnalysisTools import *

class TableMaker(AnalysisTools):
    '''Class for producing tables'''
    def __init__(self, inputFile, outFile, scale = 1, delimiter='|', doSumBG = False):
        AnalysisTools.__init__(self, inputFile, scale, savePath = '')
        self._outFile       = outFile
        self._rowList       = []
        self._columnList    = []
        self._doSumBG       = doSumBG
        self._delimiter     = delimiter


    def get_hist_dict(self, variable, histType = '1D'):
        '''
        Fill a dictionary with histograms for yields
        '''

        histDict = {}
        for data in self._datasets:
            if data in self._combineDict:
                hist = self.combine_samples(variable, data, histType)
            else:
                hist =  self.get_hist(variable, data, histType)

            #hist.Print('all')

            if hist is None: continue

            if histType == '1D':
                histDict[data] = hist
            elif histType == '2D':
                histDict[data] = hist.ProjectionX('h1_' + variable + str(level) + '_' + data, level, level+1, 'e') 

        return histDict



    def print_header(self):
        '''
        Prints header to table
        '''

        if self._delimiter == '|':

            self._outFile.write('\n   * *{0}*\n\n   {1} ** '.format(self._category, self._delimiter)),
            for data in self._columnList:

                if data in self._datasets:
                    self._outFile.write('{0} *{1}* '.format(self._delimiter, self._styleDict[data][4]))
                else:
                    self._outFile.write('{0} *{1}* '.format(self._delimiter, data))
            self._outFile.write('{0}\n'.format(self._delimiter))

        elif self._delimiter == '&':

            self._outFile.write('\n\n\\begin{table}\n\t\\ \\resizebox{\\linewidth}{!}{\\centering \n\t\\begin{tabular}{| l |' + len(self._columnList)*'| c ' + '|}\n')
            self._outFile.write('\t\hline \n\t')

            for data in self._columnList:
                if data in self._datasets:
                    #if not self._doSumBG or data in ['FCNC_M125_t', 'FCNC_M145_t', 'Signal']:
                    self._outFile.write('& ${0}$ '.format(self._styleDict[data][4]).replace('#', '\\'))
                else:
                    self._outFile.write('& {0} '.format(data))

            self._outFile.write('\\\\ \\hline \\hline \n')


    def print_table(self, histDict, doErrors = False, doEff = False, startBin=0):
        '''
        Prints a table from the entries in histograms
        '''

        delimiter = self._delimiter
        self.print_header()

        pStatement = ' {3} {0:.1f}'
        if doErrors:
            if self._delimiter == '|':
                pStatement += ' &plusmn; {1:.1f}'
            elif self._delimiter == '&':
                pStatement += ' $\pm$ {1:.1f}'
        if doEff:
            pStatement += ' ({2:.1f} %)'

        for i, row in enumerate(self._rowList):

            if row is '.': continue

            if self._delimiter == '|':
                self._outFile.write('\t'+delimiter+' '+row+' '+delimiter+' ')
            elif self._delimiter == '&':
                self._outFile.write('\t'+row+' ')

            count    = i+startBin
            totalSig = 0.
            totalBG  = 0.
            sigErr2  = 0.
            bgErr2   = 0.

            # Calculate total background and signal
            for dataset in self._datasets:

                value = histDict[dataset].GetBinContent(count)
                error = histDict[dataset].GetBinError(count)
                eff   = 1. #value/histDict[dataset].GetBinContent(1)

                if dataset in ['DATA']:
                    continue

                if dataset in ['FCNH', 'HIGGS', 'SIGNAL']: # Add extra signal names as needed
                    totalSig  += value
                    sigErr2   += pow(error, 2)

                else:
                    totalBG += value
                    bgErr2  += pow(error, 2)


            for column in self._columnList:

                if column in self._datasets:
                    value = histDict[column].GetBinContent(count)
                    error = histDict[column].GetBinError(count)
                    eff   = 1. #value/histDict[column].GetBinContent(1)

                if column in ['FCNH', 'HIGGS', 'SIGNAL'] or column[:4] == 'FCNC':
                    if doErrors:
                        #self._outFile.write(' {3} {0:.1f} $\pm$ {1:.1f} ({2:.1f} \%) '.format(value, error, eff*100, delimiter))
                        self._outFile.write(' {3} {0:.1f} $\pm$ {1:.1f} '.format(value, error, eff*100, delimiter))
                    else:
                        #self._outFile.write(' {3} {0:.1f} ({2:.1f} \%) '.format(value, error, eff*100, delimiter))
                        self._outFile.write('{3} {0:.1f}'.format(value, error, eff*100, delimiter))

                elif column is 'DATA':
                    value = histDict['DATA'].GetBinContent(count)
                    self._outFile.write(' {0} {1}'.format(delimiter, int(value)))

                elif column == 'BG':
                    self._outFile.write(pStatement.format(totalBG, sqrt(bgErr2), eff*100, delimiter))

                elif column == 'S/B':
                    if totalBG is not 0:
                        SB = totalSig/totalBG
                    else:
                        pStatement = 'nan'

                    sbErr = SB*sqrt(pow(sqrt(sigErr2)/totalSig, 2) + pow(sqrt(bgErr2)/totalBG, 2)) 
                    self._outFile.write('{0} {1} '.format(delimiter, SB*100))

                elif column == 'Significance':
                    significance    = totalSig/sqrt(totalBG+totalSig)
                    significanceErr = pow((totalSig + totalBG), -3/2)*sqrt(pow((totalSig + totalBG),2)*sigErr2 + pow(totalSig,2)*bgErr2)
                    self._outFile.write('{0} {1:.2f} '.format(delimiter, significance))

                else:
                    self._outFile.write(pStatement.format(value, error, eff*100, delimiter))


            if self._delimiter == '|':
                self._outFile.write('| \n')
            elif self._delimiter == '&':
                self._outFile.write(' \\\\ \n')

        if self._delimiter == '&':
            self._outFile.write('\t\\hline \n\t\\end{tabular}} \n')

            caption = '\t\\caption{' + self._category + '} \n'
            self._outFile.write(caption.replace('_', ' '))

            self._outFile.write('\\end{table} \n')


