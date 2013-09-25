from AnalysisTools import *

def make_index_afs(self, filePath):

    if not os.path.isfile(filePath+'/../writeIndexHTML.py'):
        os.system('cp ~/afs/public_html/writeIndexHTML.py '+filePath+'/..')
    os.system('cd '+filePath+'/..; ./writeIndexHTML.py')


def scale_to_pad(histogram):

    histMax = histogram.GetMaximum()
    scale   = r.gPad.GetUymax()/(histMax)
    histogram.Scale(scale)

    return histMax

def prep_hist(hist, yRange = (0,1)):

    hist.SetStats(0)
    hist.GetYaxis().SetLabelSize(0.08);
    hist.GetYaxis().SetTitleSize(0.09);
    hist.GetYaxis().SetTitleOffset(0.44);
    hist.GetYaxis().SetNdivisions(5);
    hist.GetYaxis().CenterTitle();
    hist.GetXaxis().SetLabelSize(0.08);
    hist.GetXaxis().SetTitleSize(0.09);
    hist.GetXaxis().SetTitleOffset(0.90);
    hist.SetMinimum(yRange[0]);
    hist.SetMaximum(yRange[1]);


def format_axis(axis, offset, title, color):

    axis.SetLineColor(color)
    axis.SetLabelColor(color)
    axis.SetTextColor(color)
    axis.SetTitleSize(0.1)
    axis.SetLabelSize(0.08)
    axis.SetTitle(title)
    axis.SetTitleOffset(offset)
    axis.CenterTitle()
    axis.Draw()


def set_hist_style(hist, dataType, styleDict, histType = '1D'):
    '''
    Set styles of histograms to be stacked
    '''

    hist.SetMarkerStyle(styleDict[dataType][3])
    hist.SetMarkerColor(styleDict[dataType][1])
    #hist.SetMarkerSize(styleDict[dataType][4])
    hist.SetMarkerSize(0.8)
    hist.SetFillStyle(styleDict[dataType][2])
    hist.SetFillColor(styleDict[dataType][1])
    hist.SetLineWidth(styleDict[dataType][0])
    hist.SetLineColor(styleDict[dataType][1])


def build_legend(hists, dataList, var, styleDict):

    legend = r.TLegend(0.91,0.45,1.0,0.89)
    legend.SetFillColor(0)
    legend.SetFillStyle(3001)
    legend.SetLineWidth(0)
    legend.SetLineColor(0)
    legend.SetTextSize(0.045)

    for data in dataList[::-1]:
        legend.AddEntry(hists[data], styleDict[data][4])

    return legend

####                                      ####  
#### Beginning of PlotProducer definition ####
####                                      ####  

class PlotProducer(AnalysisTools):
    '''For manipulating histograms'''
    def __init__(self, inputFile, scale = 1., savePath = '', isAFS = False):
        AnalysisTools.__init__(self, inputFile, scale, savePath)
        self._plotType          = '.png'
        self._makeIndex         = isAFS
        self._overlayList       = []
        self._directoryList1D   = []
        self._directoryList2D   = []
        self._variableDict      = {}

    def set_output_type(self, plotType):

        self._plotType = plotType

    def get_hist_dict(self, directory, histType = '1D'):
        '''
        Returns a dictionary with keys that are variable names and values 
        are a list of histograms.
        '''

        histDict = {}
        for var in self._variableDict[directory]:
            histList = []

            for data in self._overlayList:
                if data in self._combineDict:
                    hist = self.combine_samples(var, data, histType)
                else:
                    hist =  self.get_hist(var, data, histType)

                if hist is None: continue

                set_hist_style(hist, data, self._styleDict, histType)
                histList.append((hist, data))

            if len(histList) is 0: continue

            histDict[var] = histList

        return histDict

    def get_stack_dict(self, directory, histType = '1D'):
        '''
        Returns a dictionary with keys that are variable 
        names and values are stacks.  
        '''

        stackDict = {}
        sumDict = {}

        for var in self._variableDict[directory]:
            histList = []

            for data in self._datasets:
                if data in self._combineDict:
                    hist = self.combine_samples(var, data, histType)
                else:
                    hist =  self.get_hist(var, data, histType)


                if hist is None: continue

                set_hist_style(hist, data, self._styleDict, histType)

                histList.append((hist, data))
                
            if len(histList) is 0: continue

            stackDict[var] = self.build_stack(histList)
            
            sumHist = histList[0][0].Clone()
            for hist,data in histList[1:]:
                sumHist.Add(hist)
            set_hist_style(sumHist, 'BGERROR', self._styleDict, histType)
            sumDict[var] = sumHist

        return stackDict, sumDict



    def build_stack(self, histList):
        '''
        Combine histograms into a stack
        '''
        hStack = r.THStack('hStack', '')
        for i, (hist, data) in enumerate(histList):
            hStack.Add(hist)
            if i == 0:
                hStack.SetTitle(hist.GetTitle() + ';' + hist.GetXaxis().GetTitle() + ';' + hist.GetYaxis().GetTitle())
        return hStack



    def get_ratio(self, hist1, hist2):
        '''
        Makes a ratio of two histograms.
        '''

        nBins = hist1.GetNbinsX()
        xAxisName = hist1.GetXaxis().GetTitle()
        hRatio = r.TH1D('hRatio', ';'+xAxisName+';Data/BG', nBins, hist1.GetBinLowEdge(1), hist1.GetBinLowEdge(nBins + 1))
        hRatio.Divide(hist1, hist2)

        set_hist_style(hRatio, 'RATIO', self._styleDict)
        prep_hist(hRatio, (0, 2.499))

        ### a fudge to get labels for category hists
        if nBins == 12:
            hRatio.GetXaxis().SetLabelSize(0.15);
            for i in range(nBins):
                hRatio.GetXaxis().SetBinLabel(i+1, hist1.GetXaxis().GetBinLabel(i+1))
        ###
        
        return hRatio


    def get_cut_efficiency(self, hist, dataType):
        '''
        Generates a histogram with the cut efficiency
        for a 1-D plot.  Intended for comparison of
        background and signal.
        '''

        nBins = hist.GetNbinsX()
        xAxisName = hist.GetXaxis().GetTitle()
        hEff = r.TH1D('h1_Eff_'+dataType, ';' + xAxisName + ';#varepsilon_{cut}', nBins, hist.GetBinLowEdge(1), hist.GetBinLowEdge(nBins + 1))

        for i in range(nBins):
            sumError = r.Double(0)
            eff = hist.IntegralAndError(i+1, nBins, sumError)/hist.Integral()
            hEff.SetBinContent(i + 1, eff)
            hEff.SetBinError(i + 1, 0.025*eff)

            #hEff.SetBinError(i + 1, self.binomial_error(eff, hist.Integral()))


        set_hist_style(hEff, dataType, self._styleDict)
        prep_hist(hEff, (-0.05, 1.1))

        return hEff


    def make_stacks_by_category(self, categoryList = '', logScale = False):
        '''
        Builds and plots a stack of variables by category
        '''
        canvas = r.TCanvas('canvas', 'canvas', 800, 600)
        canvas.SetGridx()
        canvas.SetGridy()

        if logScale:
            canvas.SetLogy()

        ### Build the legend from the list of samples
        tmpHists = {}
        for data in self._datasets + self._overlayList:
            tmpHists[data] = r.TH1D('h1_tmp_' + data, ';;', 1, 0, 1)
            tmpHists[data].Fill(1)
            set_hist_style(tmpHists[data], data, self._styleDict)

        legend = build_legend(tmpHists, self._datasets, 'YieldByCut', self._styleDict)

        for directory in self._directoryList1D:
            stacks, sums = self.get_stack_dict(directory)

            self.make_save_path(self._savePath + '/' + self._category + '/' + directory)

            for var in self._variableDict[directory]:
                if logScale:
                    stacks[var].SetMaximum(stacks[var].GetMaximum()*10)
                else:
                    stacks[var].SetMaximum(stacks[var].GetMaximum()*1.5)
                stacks[var].SetMinimum(0.01)
                stacks[var].Draw('HIST')

                legend.Draw()

                canvas.SaveAs(self._savePath + '/' + self._category + '/' + directory + '/' + var + self._plotType)


    def make_overlays_1D(self, logScale = False, doRatio = True, doEff = False):
        '''
        Process to produce overlays and stacks from 1D histograms.
        '''

        ### Setting up the canvas and splitting
        ### if doing complimentary plotting
        canvas = r.TCanvas('canvas', 'canvas', 650, 700)

        if (doRatio or doEff):
            pad1 = r.TPad('pad1', '', 0.02, 0.34, 0.89, 0.98, 0)
            pad2 = r.TPad('pad2', '', 0.02, 0.02, 0.89, 0.35, 0)

            pad1.SetBottomMargin(0.)
            pad2.SetTopMargin(0.)
            pad2.SetBottomMargin(0.2)
            pad1.Draw()
            pad2.Draw()
            pad2.SetGridx()
            pad2.SetGridy()
        else:
            pad1 = r.TPad('pad1', '', 0.06, 0.02, 0.89, 0.98, 0)
            pad1.Draw()
            pad1.SetGridx()
            pad1.SetGridy()

        if logScale:
            pad1.SetLogy()

        ### Build the legend from the list of samples
        tmpHists = {}
        for data in self._datasets + self._overlayList:
            tmpHists[data] = r.TH1D('h1_tmp_' + data, ';;', 1, 0, 1)
            tmpHists[data].Fill(1)
            set_hist_style(tmpHists[data], data, self._styleDict)

        legend = build_legend(tmpHists, self._datasets, 'YieldByCut', self._styleDict)

        for data in self._overlayList:
            legend.AddEntry(tmpHists[data], self._styleDict[data][4])

        ### Starting loop over directories in histogram file
        for directory in self._directoryList1D:
            hists        = self.get_hist_dict(directory)
            stacks, sums = self.get_stack_dict(directory)

            if directory is self._directoryList1D[0]: 
                legend.AddEntry(sums[sums.keys()[0]], 'BG error')

            self.make_save_path(self._savePath + '/' + self._category + '/' + directory)

            for var in self._variableDict[directory]:

                if var not in hists.keys() or var not in stacks.keys(): continue

                pad1.cd()

                if logScale:
                    stacks[var].SetMaximum(max(stacks[var].GetMaximum(), hists[var][0][0].GetMaximum())*5)
                    stacks[var].SetMinimum(0.005)
                else:
                    stacks[var].SetMaximum(max(stacks[var].GetMaximum(), hists[var][0][0].GetMaximum())*1.25)
                    stacks[var].SetMinimum(0.00001)

                if doRatio or doEff: 
                    legend.SetX1(0.91)
                    legend.SetX2(1.0)
                    legend.SetY1(0.25)
                    legend.SetY2(0.89)


                stacks[var].Draw('HIST')
                sums[var].Draw('E2 SAME')

                if not logScale:
                    stacks[var].GetYaxis().SetTitleOffset(1.00);

                for (hist, data) in hists[var]:
                    hist.Draw('E SAME')

                legend.Draw()

                ## Draw info box ##

                r.gStyle.SetOptTitle(0)

                textBox = r.TPaveText(0.09, 0.91, 0.79, 0.96, 'NDC')
                textBox.SetFillColor(0)
                textBox.SetFillStyle(0)
                textBox.SetLineWidth(0)
                textBox.SetLineColor(0)

                if self._period is '2011':
                    textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 7 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
                elif self._period is '2012':
                    textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 8 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')

                textBox.Draw('same')

                if doRatio and not doEff:
                    pad2.cd()
                    hRatio = self.get_ratio(hists[var][0][0], sums[var])
                    hRatio.Draw("E2")

                elif doEff and not doRatio:
                    pad2.cd()
                    hEffBG  = self.get_cut_efficiency(sums[var], 'SUM_EFF')
                    hEffBG.Draw("E3")

                    for hist in hists[var]:
                        if hist[1] in ['SIGNAL', 'FCNC_M125_t']:
                            hEffSig = self.get_cut_efficiency(hist[0], 'SIG_EFF')
                            hEffSig.Draw("E3 SAME")
                            continue

                elif doEff and doRatio:
                    pad2.cd()
                    hRatio = self.get_ratio(hists[var][0][0], sums[var])
                    hRatio.Draw("E2")

                    hEffBG  = self.get_cut_efficiency(sums[var], 'SUM_EFF')
                    canvas.Update()
                    scale_to_pad(hEffBG)
                    hEffBG.Draw("E3 SAME")

                    axisEff = r.TGaxis(r.gPad.GetUxmax(), r.gPad.GetUymin(), r.gPad.GetUxmax(), r.gPad.GetUymax(), 0.0, 1.10, 510, '+L')
                    format_axis(axisEff, 0.5, '#varepsilon_{cut}', r.kBlack)

                    for hist in hists[var]:
                        if hist[1] in ['SIGNAL', 'FCNC_M125_t']:
                            hEffSig = self.get_cut_efficiency(hist[0], 'SIG_EFF')
                            canvas.Update()
                            scale_to_pad(hEffSig)
                            hEffSig.Draw("E3 SAME")
                            break

                canvas.SaveAs(self._savePath + '/' + self._category + '/' + directory + '/' + var + self._plotType)

                if doRatio:
                    hRatio.Delete()

                if doEff:
                    hEffBG.Delete()
                    hEffSig.Delete()
                

        if self._makeIndex: make_index_afs(self._savePath)


    def make_overlays_2D(self, logScale = False, doProjection = False):
        '''
        Makes an overlay of 2D histograms.  As in the 1D case,
        there are samples that are summed into a 'stack' and samples
        to be overlayed.  Adding projection functionality would be
        nice...
        '''

        canvas = r.TCanvas('canvas', 'canvas', 650, 700)

        if doProjection:
            pad1 = r.TPad('pad1', '', 0.02, 0.34, 0.90, 0.98, 0)
            pad2 = r.TPad('pad2', '', 0.02, 0.02, 0.90, 0.35, 0)
            pad1.SetBottomMargin(0.)
            pad2.SetTopMargin(0.)
            pad2.SetBottomMargin(0.2)
            pad1.Draw()
            pad2.Draw()
            pad2.SetGridx()
            pad2.SetGridy()
        else:
            pad1 = r.TPad('pad1', '', 0.02, 0.65, 0.96, 0.96, 0)
            pad2 = r.TPad('pad2', '', 0.02, 0.335, 0.96, 0.645, 0)
            pad3 = r.TPad('pad3', '', 0.02, 0.02, 0.96, 0.33, 0)

            pad2.SetTopMargin(0.)
            pad3.SetTopMargin(0.)
            pad3.SetBottomMargin(0.2)

            pad1.Draw()
            pad2.Draw()
            pad3.Draw()

            pad1.SetGridx()
            pad1.SetGridy()
            pad2.SetGridx()
            pad2.SetGridy()
            pad3.SetGridx()
            pad3.SetGridy()

        if logScale:
            pad1.SetLogz()
            pad2.SetLogz()
            pad3.SetLogz()


        textBox = r.TPaveText(0.12, 0.95, 0.67, 0.98, 'NDC')
        textBox.SetFillColor(0)
        textBox.SetFillStyle(0)
        textBox.SetLineWidth(0)
        textBox.SetLineColor(0)

        if self._period is '2011':
            textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 7 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
        elif self._period is '2012':
            textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 8 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')

        for directory in self._directoryList2D:
            hists           = self.get_hist_dict(directory, '2D')
            stacks, sums    = self.get_stack_dict(directory, '2D')

            #print hists, stacks, sums

            self.make_save_path(self._savePath + '/' + self._category + '/' + directory)

            for var in self._variableDict[directory]:

                ## Draw info box ##

                if var not in hists.keys() and var not in stacks.keys(): continue

                r.gStyle.SetOptTitle(0)

                textBox.Draw('same')
                  
                #sigBox = r.TPaveText(0.95, 0.925, 0.7, 0.96, 'NDC')
                #sigBox.SetFillColor(0)
                #sigBox.SetFillStyle(0)
                #sigBox.SetLineWidth(0)
                #sigBox.SetLineColor(0)

                for (hist, data) in hists[var]:
                    if data in ['Signal', 'FCNH']:
                        pad1.cd()
                        hist.GetYaxis().SetTitle('')
                        hist.GetYaxis().SetLabelOffset(0.01)
                        hist.GetYaxis().SetLabelSize(0.06)

                        hist.GetXaxis().SetTitle('')
                        hist.GetXaxis().SetLabelSize(0.)
                        hist.Draw('COLZ')

                    if data in ['DATA', 'DATA_MUON','DATA_ELECTRON', 'DATA_MUEG']:
                        pad2.cd()
                        hist.GetYaxis().CenterTitle()
                        hist.GetYaxis().SetTitleSize(0.09)
                        hist.GetYaxis().SetTitleOffset(0.6)
                        hist.GetYaxis().SetLabelOffset(0.01)
                        hist.GetYaxis().SetLabelSize(0.06)

                        hist.GetXaxis().SetTitle('')
                        hist.GetXaxis().SetLabelSize(0.)
                        hist.Draw('COLZ')

                pad3.cd()
                sums[var].GetYaxis().SetTitle('')
                sums[var].GetYaxis().SetLabelOffset(0.01)
                sums[var].GetYaxis().SetLabelSize(0.06)

                sums[var].GetXaxis().SetLabelOffset(0.03)
                sums[var].GetXaxis().SetLabelSize(0.065)
                sums[var].GetXaxis().SetTitleSize(0.09)
                sums[var].GetXaxis().SetTitleOffset(1.)
                sums[var].Draw('COLZ')


                canvas.SaveAs(self._savePath + '/' + self._category + '/' + directory + '/' + var + self._plotType)

                #legend.Draw()

####                                      ####  
#### End of PlotProducer class definition ####
####                                      ####  

def plotter_wrapper(plotter, category, inputPath, outputPath, do1D, do2D, doLog):

    plotter.set_input_file(inputPath)
    plotter.set_save_path(outputPath)
    plotter._category = category

    if do1D:
        plotter.make_overlays_1D(logScale = doLog, doRatio = True, doEff = False)
    if do2D:
        plotter.make_overlays_2D(logScale = False, doProjection = False)

