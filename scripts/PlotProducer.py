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
    hist.GetYaxis().SetLabelSize(0.08)
    hist.GetYaxis().SetTitleSize(0.09)
    hist.GetYaxis().SetTitleOffset(0.44)
    hist.GetYaxis().SetNdivisions(5)
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetLabelSize(0.08)
    hist.GetXaxis().SetTitleSize(0.09)
    hist.GetXaxis().SetTitleOffset(0.90)
    hist.SetMinimum(yRange[0])
    hist.SetMaximum(yRange[1])


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

    if dataType.split('_')[0] in ['Fakes', 'eFakes', 'muFakes', 'llFakes'] and dataType not in ['Fakes', 'eFakes', 'muFakes', 'llFakes']:
        dataType = dataType.split('_')[1]

    hist.SetMarkerStyle(styleDict[dataType][3])
    hist.SetMarkerColor(styleDict[dataType][1])
    #hist.SetMarkerSize(styleDict[dataType][4])
    hist.SetMarkerSize(0.8)
    hist.SetFillStyle(styleDict[dataType][2])
    hist.SetFillColor(styleDict[dataType][1])
    hist.SetLineWidth(styleDict[dataType][0])
    hist.SetLineColor(styleDict[dataType][1])


def build_legend(hists, dataList, styleDict):

    for data in dataList:
        hists[data] = r.TH1D('h1_tmp_' + data, ';;', 1, 0, 1)
        hists[data].Fill(1)
        set_hist_style(hists[data], data, styleDict)

    legend = r.TLegend(0.6,0.55,0.79,0.88)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineWidth(0)
    legend.SetLineColor(0)
    legend.SetTextSize(0.03)

    for data in dataList[::-1]:
        if data.split('_')[0] in ['Fakes', 'eFakes', 'muFakes', 'llFakes'] and data not in ['Fakes', 'eFakes', 'muFakes', 'llFakes']:
            dataName = data.split('_')[1]
        else:
            dataName = data

        legend.AddEntry(hists[data], styleDict[dataName][4])

    return legend


####                                      ####  
#### Beginning of PlotProducer definition ####
####                                      ####  

class PlotProducer(AnalysisTools):
    '''For manipulating histograms'''
    def __init__(self, inputFile, scale = 1., savePath = '', isAFS = False, drawNorm = False):
        AnalysisTools.__init__(self, inputFile, scale, savePath)
        self._plotType          = '.png'
        self._makeIndex         = isAFS
        self._drawNormalized    = drawNorm
        self._overlayList       = []
        self._directoryList1D   = []
        self._directoryList2D   = []
        self._variableDict      = {}

    def draw_normalized(self, drawNorm):
        self._drawNormalized = drawNorm

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
                hist = self.combine_samples(var, data, histType)

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
                hist = self.combine_samples(var, data, histType)

                if hist is None: continue

                if var == 'HT':
                    hist.GetXaxis().SetTitle('H_{T} [GeV]')
                if var == 'Met':
                    hist.GetXaxis().SetTitle('MET [GeV]')
                if var in ['DileptonMass21', 'DileptonZMass21', 'DileptonHiggsMass21']:
                    #hist.GetYaxis().SetTitle('Entries / 5 GeV/c^{2}')
                    hist.GetXaxis().SetTitle('M_{ll} [GeV/c^{2}]')
                if var == 'TrileptonMVsDileptonMOS':
                    hist.GetXaxis().SetTitle('M_{OS} [GeV/c^{2}]')
                    hist.GetYaxis().SetTitle('M_{lll} [GeV/c^{2}]')

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


    def get_difference(self, hist1, hist2):
        '''
        Makes a difference of two histograms.
        '''

        nBins = hist1.GetNbinsX()
        xAxisName = hist1.GetXaxis().GetTitle()
        hDiff = r.TH1D('hDiff', ';{0};Data/BG'.format(xAxisName), nBins, hist1.GetBinLowEdge(1), hist1.GetBinLowEdge(nBins + 1))
        hDiff.Add(hist1, hist2, -1.)

        set_hist_style(hDiff, 'RATIO', self._styleDict)
        prep_hist(hDiff, (0, 2.499))

        ### a fudge to get labels for category hists
        if hist1.GetName() in ['h1_LeptonFlavor', 'h1_LeptonCharge']:
            hDiff.GetXaxis().SetLabelSize(0.15)
            for i in range(nBins):
                hDiff.GetXaxis().SetBinLabel(i+1, hist1.GetXaxis().GetBinLabel(i+1))
        ###
        
        return hDiff

    def get_ratio(self, hist1, hist2):
        '''
        Makes a ratio of two histograms.
        '''

        nBins = hist1.GetNbinsX()
        xAxisName = hist1.GetXaxis().GetTitle()
        hRatio = hist1.Clone()
        hRatio.Divide(hist2)
        hRatio.GetYaxis().SetTitle('Data/BG')

        #hRatio = r.TH1D('hRatio', ';{0};Data/BG'.format(xAxisName), nBins, hist1.GetBinLowEdge(1), hist1.GetBinLowEdge(nBins + 1))
        #hRatio.Divide(hist1, hist2)

        set_hist_style(hRatio, 'RATIO', self._styleDict)
        prep_hist(hRatio, (0.5, 2.))

        ### a fudge to get labels for category hists ###
        if hist1.GetName() in ['h1_LeptonFlavor', 'h1_LeptonCharge']:
            hRatio.GetXaxis().SetLabelSize(0.15)
            for i in range(nBins):
                hRatio.GetXaxis().SetBinLabel(i+1, hist1.GetXaxis().GetBinLabel(i+1))
        ###
        
        return hRatio

    def get_significance(self, signal, background, dataType):
        '''
        Generates a histogram with the cut efficiency
        for a 1-D plot.  Intended for comparison of
        background and signal.
        '''

        nBins = signal.GetNbinsX()
        xAxisName = signal.GetXaxis().GetTitle()
        hSig = r.TH1D('h1_Sig_'+dataType, ';' + xAxisName + ';#varepsilon_{cut}', nBins, signal.GetBinLowEdge(1), signal.GetBinLowEdge(nBins + 1))

        for i in range(nBins):
            sig     = signal.Integral(i+1, nBins)
            bg      = background.Integral(i+1, nBins)

            #print sig + bg
            if sig + bg > 0.:
                signif  = sig/sqrt(sig + bg)
            else:
                signif  = 0.

            hSig.SetBinContent(i + 1, signif)
            hSig.SetBinError(i + 1, 0.025*signif)

        set_hist_style(hSig, dataType, self._styleDict)
        prep_hist(hSig, (-0.05, hSig.GetMaximum()*1.1))

        return hSig


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

    def make_stacks_by_category(self, logScale = False):
        '''
        Builds and plots a stack of variables by category
        '''
        canvas = r.TCanvas('canvas', 'canvas', 650, 700)
        #canvas.SetGridx()
        #canvas.SetGridy()
        canvas.SetRightMargin(0.18)
        canvas.SetLeftMargin(0.15)

        if logScale:
            canvas.SetLogy()

        ### Build the legend from the list of samples
        tmpHists = {}
        legend = build_legend(tmpHists, self._datasets, self._styleDict)
        legend.SetX1(0.83)
        legend.SetX2(0.99)
        legend.SetY1(0.45)
        legend.SetY2(0.89)
        legend.SetTextSize(0.03)

        for directory in self._directoryList1D:
            stacks, sums = self.get_stack_dict(directory)

            if directory is self._directoryList1D[0]: 
                legend.AddEntry(sums.values()[0], 'BG uncertainty')

            self.make_save_path(self._savePath + '/' + self._category + '/' + directory)

            for var in self._variableDict[directory]:
                if var not in stacks.keys(): continue

                if logScale:
                    stacks[var].SetMaximum(stacks[var].GetMaximum()*10)
                    stacks[var].SetMinimum(0.09)
                else:
                    stacks[var].SetMaximum(stacks[var].GetMaximum()*1.3)

                #stacks[var].GetYaxis().SetTitleOffset(1.3);
                #stacks[var].GetYaxis().SetTitleSize(0.04);

                stacks[var].Draw('HIST')
                sums[var].Draw('E2 SAME')

                legend.Draw()

                ## Draw info box ##
                r.gStyle.SetOptTitle(0)
                textBox = r.TPaveText(0.09, 0.91, 0.81, 0.98, 'NDC')
                textBox.SetFillColor(0)
                textBox.SetFillStyle(0)
                textBox.SetLineWidth(0)
                #textBox.SetLineColor(0)
                textBox.SetTextSize(0.025)

                if self._period is '2011':
                    textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 7 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
                elif self._period is '2012':
                    textBox.AddText('#scale[1.2]{CMS preliminary, 8 TeV, ' + '{0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')

                textBox.Draw('same')

                canvas.SaveAs('{0}/{1}/{2}/{3}{4}'.format(self._savePath, self._category, directory, var, self._plotType))


    def make_overlays_diff(self, hists, directory, outputName, logScale = False, doRatio = False, doEff = False):
        '''
        Process to produce overlays of different histograms.
        '''

        ### Setting up the canvas and splitting
        ### if doing complimentary plotting
        canvas = r.TCanvas('canvas', 'canvas', 650, 700)

        if doRatio:
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

        self.make_save_path('{0}/{1}/{2}'.format(self._savePath, self._category, directory))

        histList = []
        dataNames = []
        for samples, var in hists:

            drawOption = 'E'
            if len(samples) == 1:
                hist = self.combine_samples(var[0], samples[0])
                if hist is None: 
                    continue
                else:
                    dataNames.append(samples[0])
                    set_hist_style(hist, samples[0], self._styleDict)
                    histList.append((hist, drawOption))
            else:
                drawOption = 'HIST'
                stackList = []
                for j,sample in enumerate(samples):
                    tmpHist = self.combine_samples(var[j], sample)

                    if tmpHist is None: 
                        continue
                    else:
                        dataNames.append(sample)
                        set_hist_style(tmpHist, sample, self._styleDict)
                        stackList.append((tmpHist, sample))

                hist = self.build_stack(stackList)

                if hist is None: 
                    continue
                else:
                    histList.append((hist, drawOption))

        tmpHists = {}
        legend = build_legend(tmpHists, list(set(dataNames)), self._styleDict)

        # Get maximum
        max = 0.
        for (hist, opt) in histList:
            histMax = hist.GetMaximum()
            if histMax > max:
                max = histMax

        pad1.cd()
        isBlank = True
        for (hist, opt) in histList:
            if isBlank:
                hist.SetMaximum(max*1.25)
                hist.Draw(opt)

                if hist.GetYaxis():
                    hist.GetYaxis().SetTitleOffset(1.3)
                    hist.GetYaxis().SetTitleSize(0.04)
                    hist.GetXaxis().SetTitleOffset(0.9)
                    hist.GetXaxis().SetTitleSize(0.04)

                isBlank = False
            else:
                hist.Draw('{0} SAME'.format(opt))

        ## Draw info box ##
        r.gStyle.SetOptTitle(0)
        textBox = r.TPaveText(0.09, 0.91, 0.91, 0.98, 'NDC')
        textBox.SetFillColor(0)
        textBox.SetFillStyle(0)
        textBox.SetLineWidth(0)
        textBox.SetLineColor(0)
        textBox.SetTextSize(0.03)

        if self._period is '2011':
            textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 7 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
        elif self._period is '2012':
            textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 8 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')

        legend.Draw('SAME')
        textBox.Draw('SAME')

        canvas.SaveAs('{0}/{1}/{2}/{3}{4}'.format(self._savePath, self._category, directory, outputName, self._plotType))

    def make_overlays_1D_simple(self, logScale = False):
        '''
        Process to produce simple overlays from 1D histograms.
        '''

        ### Setting up the canvas and splitting
        ### if doing complimentary plotting
        canvas = r.TCanvas('canvas', 'canvas', 650, 700)
        canvas.SetRightMargin(0.18)
        canvas.SetLeftMargin(0.15)

        if logScale:
            canvas.SetLogy()

        ### Build the legend from the list of samples
        tmpHists = {}
        legend = build_legend(tmpHists, self._overlayList, self._styleDict)

        ### Starting loop over directories in histogram file
        for directory in self._directoryList1D:
            hists        = self.get_hist_dict(directory)

            self.make_save_path('{0}/{1}/{2}'.format(self._savePath, self._category, directory))

            for var in self._variableDict[directory]:
                if var not in hists.keys(): continue

                #if not logScale or not doRatio:
                    #stacks[var].GetYaxis().SetTitleOffset(1.3);
                    #stacks[var].GetYaxis().SetTitleSize(0.04);
                    #stacks[var].GetXaxis().SetTitleOffset(0.9);
                    #stacks[var].GetXaxis().SetTitleSize(0.04);
                    #stacks[var].Draw('HIST')

                hist_max = 0
                for (hist, data) in hists[var]:
                    if hist.GetMaximum() > hist_max:
                        hist_max = hist.GetMaximum()

                if logScale:
                    hists[var][0][0].SetMaximum(2.*hist_max)
                    hists[var][0][0].SetMinimum(0.1)
                else:               
                    hists[var][0][0].SetMaximum(1.25*hist_max)
                    hists[var][0][0].SetMinimum(0.00001)

                for i,(hist, data) in enumerate(hists[var]):
                    if i == 0:
                        drawOpt = 'HIST'
                    else:
                        drawOpt = 'HIST SAME'

                    hist.SetFillStyle(0)

                    if self._drawNormalized:
                        hist.DrawNormalized(drawOpt)
                    else:
                        hist.Draw(drawOpt)


                legend.SetX1(0.85)
                legend.SetX2(1.0)
                legend.SetY1(0.55)
                legend.SetY2(0.89)
                legend.SetTextSize(0.03)
                legend.Draw()

                ## Draw info box ##
                r.gStyle.SetOptTitle(0)
                textBox = r.TPaveText(0.09, 0.91, 0.91, 0.98, 'NDC')
                textBox.SetFillColor(0)
                textBox.SetFillStyle(0)
                textBox.SetLineWidth(0)
                textBox.SetLineColor(0)
                textBox.SetTextSize(0.03)

                if self._period is '2011':
                    textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 7 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
                elif self._period is '2012':
                    textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 8 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
                textBox.Draw('same')

                canvas.SaveAs(self._savePath + '/' + self._category + '/' + directory + '/' + var + self._plotType)

    def make_overlays_1D(self, logScale = False, doRatio = True, doEff = False, doDiff = False, rebin = 1.):
        '''
        Process to produce overlays and stacks from 1D histograms.
        '''

        ### Setting up the canvas and splitting
        ### if doing complimentary plotting
        canvas = r.TCanvas('canvas', 'canvas', 650, 700)

        if (doRatio or doEff or doDiff):
            pad1 = r.TPad('pad1', '', 0.00, 0.35, 0.99, 0.99, 0)
            pad2 = r.TPad('pad2', '', 0.00, 0.02, 0.99, 0.35, 0)

            #pad1.SetBottomMargin(0.02)
            #pad2.SetTopMargin(0.02)
            pad1.SetBottomMargin(0.)
            pad2.SetTopMargin(0.)
            pad2.SetBottomMargin(0.2)
            pad1.Draw()
            pad2.Draw()
            pad2.SetGridx()
            pad2.SetGridy()
        else:
            pad1 = r.TPad('pad1', '', 0.00, 0.02, 1., 1., 0)
            pad1.Draw()
            pad1.SetLeftMargin(0.15)
            pad1.SetTickx()
            pad1.SetTicky()
            #pad1.SetGridx()
            #pad1.SetGridy()

        if logScale:
            pad1.SetLogy()

        ### Build the legend from the list of samples
        tmpHists = {}
        legend = build_legend(tmpHists, self._overlayList+self._datasets, self._styleDict)

        ### Starting loop over directories in histogram file
        for directory in self._directoryList1D:
            hists        = self.get_hist_dict(directory)
            stacks, sums = self.get_stack_dict(directory)

            if directory is self._directoryList1D[0]: 
                legend.AddEntry(sums.values()[0], 'BG uncertainty')

            self.make_save_path('{0}/{1}/{2}'.format(self._savePath, self._category, directory))

            for var in self._variableDict[directory]:
                if var not in hists.keys() or var not in stacks.keys(): continue

                pad1.cd()

                if logScale:
                    stacks[var].SetMaximum(max(stacks[var].GetMaximum(), hists[var][0][0].GetMaximum())*5)
                    stacks[var].SetMinimum(0.2)
                else:
                    stacks[var].SetMaximum(max(stacks[var].GetMaximum(), hists[var][0][0].GetMaximum())*1.3)
                    stacks[var].SetMinimum(0.00001)

                #if doRatio or doEff: 
                #    legend.SetX1(0.91)
                #    legend.SetX2(1.0)
                #    legend.SetY1(0.25)
                #    legend.SetY2(0.89)

                if stacks[var]:
                    stacks[var].Draw('HIST')    
                else:
                    continue

                ### Hack for flavor histograms ###
                if var == 'LeptonFlavor':
                    if self._category[:2] == 'ss':
                        stacks[var].GetXaxis().SetRange(1, 4)
                    if self._category[:2] == '3l':
                        stacks[var].GetXaxis().SetRange(5, 12)

                if not logScale or not doRatio:

                    if var == 'HT':
                        stacks[var].GetXaxis().SetRangeUser(0., 1500.)
                        stacks[var].GetXaxis().SetTitle('H_{T} [GeV]')

                    stacks[var].GetYaxis().SetTitleOffset(1.95)
                    stacks[var].GetYaxis().SetTitleSize(0.04)
                    stacks[var].GetXaxis().SetTitleOffset(0.9)
                    stacks[var].GetXaxis().SetTitleSize(0.04)
                    stacks[var].Draw('HIST')

                sums[var].Draw('E2 SAME')

                for (hist, data) in hists[var]:
                    if data == 'FCNH':
                        hist.Draw('H SAME')
                    else:
                        hist.Draw('E SAME')

                legend.Draw()

                ## Draw info box ##
                r.gStyle.SetOptTitle(0)
                cmsBox = r.TPaveText(0.16, 0.75, 0.51, 0.87, 'NDC')
                cmsBox.SetFillColor(0)
                cmsBox.SetFillStyle(0)
                cmsBox.SetLineWidth(0)
                cmsBox.SetLineColor(0)
                cmsBox.SetTextSize(0.028)

                cmsBox.AddText('#scale[2.5]{CMS}')
                cmsBox.AddText('#it{Preliminary}')
                cmsBox.AddText('#bf{#color[2]{' + categories[self._category] + '}}')

                cmsBox.Draw('same')

                energy = ''
                if self._period is '2011':
                    energy = '7'
                elif self._period is '2012':
                    energy = '8'

                lumiBox = r.TPaveText(0.51, 0.91, 0.89, 0.96, 'NDC')
                lumiBox.SetFillColor(0)
                lumiBox.SetFillStyle(0)
                lumiBox.SetLineWidth(0)
                lumiBox.SetLineColor(0)
                lumiBox.SetTextSize(0.035)

                lumiBox.AddText('{0:.1f}'.format(self._scale) + ' fb^{-1} (' + energy + ' TeV)')       

                lumiBox.Draw('same')

                if doDiff:
                    doRatio = False
                    doEff   = False
                    pad2.cd()
                    hDiff = hists[var][0][0].Clone()
                    hDiff.Add(sums[var], -1.)
                    hDiff.Draw("E2")

                if doRatio and not doEff:
                    pad2.cd()
                    hRatio = self.get_ratio(hists[var][0][0], sums[var])
                    hRatio.Draw("E2")

                    ### Hack for flavor histograms ###
                    if var == 'LeptonFlavor':
                        if self._category[:2] == 'ss':
                            hRatio.GetXaxis().SetRange(1, 4)
                        if self._category[:2] == '3l':
                            hRatio.GetXaxis().SetRange(5, 12)

                elif doEff and not doRatio:
                    pad2.cd()
                    hEffBG  = self.get_cut_efficiency(sums[var], 'SUM_EFF')
                    hEffBG.Draw("E3")

                    hEffSig = self.get_cut_efficiency(hists[var][1][0], 'SIG_EFF')
                    canvas.Update()
                    scale_to_pad(hEffSig)
                    hEffSig.Draw("E3 SAME")

                    hSignif = self.get_significance(hists[var][1][0], sums[var], 'SIGNIFICANCE')
                    canvas.Update()
                    scale_to_pad(hSignif)
                    hSignif.Draw("E3 SAME")

                    axisEff = r.TGaxis(r.gPad.GetUxmax(), r.gPad.GetUymin(), r.gPad.GetUxmax(), r.gPad.GetUymax(), 0.0, 1.10, 510, '+L')
                    format_axis(axisEff, 0.5, '#varepsilon_{cut}', r.kBlack)


                elif doEff and doRatio:
                    pad2.cd()
                    hRatio = self.get_ratio(hists[var][0][0], sums[var])
                    hRatio.Draw("E2")

                    hEffBG  = self.get_cut_efficiency(sums[var], 'SUM_EFF')
                    canvas.Update()
                    scale_to_pad(hEffBG)
                    hEffBG.Draw("L SAME")

                    hEffSig = self.get_cut_efficiency(hists[var][1][0], 'SIG_EFF')
                    canvas.Update()
                    scale_to_pad(hEffSig)
                    hEffSig.Draw("L SAME")

                    hSignif = self.get_significance(hists[var][1][0], sums[var], 'SIGNIFICANCE')
                    canvas.Update()
                    scale_to_pad(hSignif)
                    hSignif.Draw("L SAME")

                    axisEff = r.TGaxis(r.gPad.GetUxmax(), r.gPad.GetUymin(), r.gPad.GetUxmax(), r.gPad.GetUymax(), 0.0, 1.10, 510, '+L')
                    format_axis(axisEff, 0.5, '#varepsilon_{cut}', r.kBlack)

                    cutBox = r.TPaveText(0.75, 0.35, 0.89, 0.75, 'NDC')
                    cutBox.SetFillColor(0)
                    cutBox.SetFillStyle(0)
                    cutBox.SetLineWidth(1)
                    cutBox.SetLineColor(1)

                    cutBox.AddText('sig_{max} = ' + '{0:.1f}'.format(hSignif.GetMaximum()))
                    cutBox.AddText('eff_{bg} = ' + '{0:.1f}'.format(hEffBG.GetBinContent(hSignif.GetMaximumBin())))
                    cutBox.AddText('eff_{sig} = ' + '{0:.1f}'.format(hEffSig.GetBinContent(hSignif.GetMaximumBin())))

                    cutBox.Draw("SAME")

                canvas.SaveAs(self._savePath + '/' + self._category + '/' + directory + '/' + var + self._plotType)

                if doDiff:
                    hDiff.Delete()
                if doRatio:
                    hRatio.Delete()

                if doEff: 
                    hEffBG.Delete()
                    hEffSig.Delete()
                    hSignif.Delete()

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
            pad1 = r.TPad('pad1', '', 0.0, 0.34, 0.90, 0.98, 0)
            pad2 = r.TPad('pad2', '', 0.0, 0.0, 0.90, 0.35, 0)
            pad1.SetBottomMargin(0.)
            pad2.SetTopMargin(0.)
            pad2.SetBottomMargin(0.2)
            pad1.Draw()
            pad2.Draw()
            #pad2.SetGridx()
            #pad2.SetGridy()
        else:
            pad1 = r.TPad('pad1', '', 0.0, 0.65, 0.99, 0.96, 0)
            pad2 = r.TPad('pad2', '', 0.0, 0.335, 0.99, 0.645, 0)
            pad3 = r.TPad('pad3', '', 0.0, 0.0, 0.99, 0.33, 0)

            pad1.SetRightMargin(0.15)
            pad2.SetRightMargin(0.15)
            pad3.SetRightMargin(0.15)

            pad2.SetTopMargin(0.)
            pad3.SetTopMargin(0.)
            pad3.SetBottomMargin(0.25)

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


        textBox = r.TPaveText(0.09, 0.94, 0.89, 0.98, 'NDC')
        textBox.SetFillColor(0)
        textBox.SetFillStyle(0)
        textBox.SetLineWidth(0)
        textBox.SetLineColor(0)
        textBox.SetTextSize(0.023)

        if self._period is '2011':
            textBox.AddText('#scale[1.2]{CMS preliminary, #sqrt{s} = 7 TeV, #it{L}_{int}' + ' = {0:.1f}'.format(self._scale) + ' fb^{-1}       #bf{#color[2]{' + categories[self._category] + '}}}')
        elif self._period is '2012':
            textBox.AddText('#scale[1.2]{CMS preliminary, ' + '{0:.1f}'.format(self._scale) + ' fb^{-1} (8 TeV),        #bf{#color[2]{' + categories[self._category] + '}}}')

        for directory in self._directoryList2D:
            hists           = self.get_hist_dict(directory, '2D')
            stacks, sums    = self.get_stack_dict(directory, '2D')

            #print hists, stacks, sums

            self.make_save_path(self._savePath + '/' + self._category + '/' + directory)

            for var in self._variableDict[directory]:

                ## Draw info box ##

                if var not in hists.keys() or var not in stacks.keys(): continue

                r.gStyle.SetOptTitle(0)

                textBox.Draw('same')
                  
                idBox = r.TPaveText(0.7, 0.24, 0.89, 0.4, 'NDC')
                idBox.SetFillColor(0)
                #idBox.SetFillStyle(1)
                idBox.SetLineWidth(1)
                idBox.SetLineColor(1)
                idBox.SetTextSize(0.1)
                idBox.SetTextColor(r.kBlue)

                for (hist, data) in hists[var]:
                    if data in ['Signal', 'FCNH']:
                        pad1.cd()
                        hist.GetYaxis().SetTitle('')
                        hist.GetYaxis().SetLabelOffset(0.01)
                        hist.GetYaxis().SetLabelSize(0.06)

                        hist.GetXaxis().SetTitle('')
                        hist.GetXaxis().SetLabelSize(0.)

                        hist.GetZaxis().SetLabelSize(0.09)

                        # Hacks
                        if var == 'MetVsHT':
                            hist.GetXaxis().SetRangeUser(0., 500.)
                            hist.GetYaxis().SetRangeUser(0., 200.)

                        hist.Draw('COLZ')

                        sigBox = idBox.Clone()
                        sigBox.AddText('signal')
                        #sigBox.Draw('SAME')

                    if data in ['DATA', 'DATA_MUON','DATA_ELECTRON', 'DATA_MUEG']:
                        pad2.cd()
                        hist.GetYaxis().CenterTitle()
                        hist.GetYaxis().SetTitleSize(0.13)
                        hist.GetYaxis().SetTitleOffset(0.35)
                        hist.GetYaxis().SetLabelOffset(0.01)
                        hist.GetYaxis().SetLabelSize(0.06)

                        # Hacks
                        if var == 'TrileptonMVsDileptonMOS':
                            hist.GetYaxis().SetTitle('M_{lll} [GeV/c^{2}]')

                        if var == 'MetVsHT':
                            hist.GetYaxis().SetTitle('MET [GeV]')
                            hist.GetXaxis().SetRangeUser(0., 500.)
                            hist.GetYaxis().SetRangeUser(0., 200.)


                        hist.GetXaxis().SetTitle('')
                        hist.GetXaxis().SetLabelSize(0.)

                        hist.GetZaxis().SetLabelSize(0.09)

                        hist.Draw('COLZ')

                        dataBox = idBox.Clone()
                        dataBox.AddText('data')
                        #dataBox.Draw('SAME')

                pad3.cd()
                sums[var].GetYaxis().SetTitle('')
                sums[var].GetYaxis().SetLabelOffset(0.01)
                sums[var].GetYaxis().SetLabelSize(0.06)

                sums[var].GetXaxis().SetLabelOffset(0.03)
                sums[var].GetXaxis().SetLabelSize(0.065)
                sums[var].GetXaxis().SetTitleSize(0.10)
                sums[var].GetXaxis().SetTitleOffset(1.)

                if var == 'MetVsHT':
                    sums[var].GetXaxis().SetRangeUser(0., 500.)
                    sums[var].GetYaxis().SetRangeUser(0., 200.)
                    sums[var].GetXaxis().SetTitle('H_{T} [GeV]')

                sums[var].GetZaxis().SetLabelSize(0.09)

                sums[var].Draw('COLZ')

                bgBox = idBox.Clone()
                bgBox.AddText('background')
                #bgBox.Draw('SAME')


                #canvas.SaveAs(self._savePath + '/' + self._category + '/' + directory + '/' + var + self._plotType)
                canvas.SaveAs('{0}/{1}/{2}/{3}{4}'.format(self._savePath, self._category, directory, var, self._plotType))

                #legend.Draw()

    def make_overlay_2D_projections(self, varName, samples, directory, axis = 'X', projection = 'normal'):
        ### Makes 1D histograms by projecting onto one of the axes of a 2D histogram ###

        canvas = r.TCanvas('canvas', 'canvas', 650, 700)
        canvas.SetLogy()

        self.make_save_path(self._savePath + '/' + self._category + '/' + directory)

        hist = self.combine_samples(varName, samples[0], '2D')
        for sample in samples[1:]:
            hist.Add(self.combine_samples(varName, sample, '2D'))

        legend = r.TLegend(0.55,0.6,.95,0.93)
        legend.SetFillColor(0)
        #legend.SetFillStyle(1)
        legend.SetLineWidth(1)
        legend.SetLineColor(1)
        legend.SetTextSize(0.03)

        nBins   = hist.GetYaxis().GetNbins()
        xMin    = hist.GetYaxis().GetXmax()
        xMax    = hist.GetYaxis().GetXmax()

        for i in range(nBins):
            if axis == 'X':
                if projection == 'normal':
                    h_Proj = hist.ProjectionX('h_{0}_{1}'.format(varName, i+1), i, nBins-1)
                elif projection == 'reverse':
                    h_Proj = hist.ProjectionX('h_{0}_{1}'.format(varName, i+1), 0, nBins-i)
                elif projection == 'exclusive':
                    h_Proj = hist.ProjectionX('h_{0}_{1}'.format(varName, i+1), i, i+1)
            if axis == 'Y':
                if projection == 'normal':
                    h_Proj = hist.ProjectionY('h_{0}_{1}'.format(varName, i+1), i, nBins-1)
                elif projection == 'reverse':
                    h_Proj = hist.ProjectionY('h_{0}_{1}'.format(varName, i+1), 0, nBins-i)
                elif projection == 'exclusive':
                    h_Proj = hist.ProjectionY('h_{0}_{1}'.format(varName, i+1), i, i+1)

            #h_Proj.SetTitle('')
            h_Proj.SetLineWidth(3)
            if i != 4:
                h_Proj.SetLineColor(i+1)
            else:
                h_Proj.SetLineColor(38)

            if projection == 'normal':
                if i is 0:
                    legend.AddEntry(h_Proj, 'IsoRel > 0.12')
                else:
                    legend.AddEntry(h_Proj, 'IsoRel > {0:.2f}'.format(i*xMax/nBins))
            elif projection == 'reverse':
                if i is 0:
                    legend.AddEntry(h_Proj, 'IsoRel > 0.12')
                else:
                    legend.AddEntry(h_Proj, '0.12 < IsoRel < {0:.2f}'.format(xMax*(1 - float(i)/nBins)))
            elif projection == 'exclusive':
                if i is 0:
                    legend.AddEntry(h_Proj, '0.12 < IsoRel < {0:.2f}'.format(xMax/nBins))
                else:
                    legend.AddEntry(h_Proj, '{0:.2f} < IsoRel < {1:.2f}'.format(i*xMax/nBins, (i + 1)*xMax/nBins))

            if i is 0:
                h_Proj.SetMaximum(1.50*h_Proj.GetMaximum())
                h_Proj.DrawNormalized('HIST')
            else:
                h_Proj.DrawNormalized('HIST SAME')

            legend.Draw()

        if projection == 'normal':
            canvas.Print('{0}/{1}/{2}/{3}_normal.png'.format(self._savePath, self._category, directory, varName))
        elif projection == 'reverse':
            canvas.Print('{0}/{1}/{2}/{3}_reverse.png'.format(self._savePath, self._category, directory, varName))
        elif projection == 'exclusive':
            canvas.Print('{0}/{1}/{2}/{3}_exclusive.png'.format(self._savePath, self._category, directory, varName))

####                                      ####  
#### End of PlotProducer class definition ####
####                                      ####  

def plotter_wrapper(plotter, category, inputPath, outputPath, do1D, do2D, doSimple, log, ratios, eff, diff):

    plotter.set_input_file(inputPath)
    plotter.set_save_path(outputPath)
    plotter.set_category(category)

    if do1D:
        if doSimple:
            plotter.make_overlays_1D_simple(logScale = log)
        else:
            plotter.make_overlays_1D(logScale = log, doRatio = ratios, doEff = eff, doDiff = diff)
    if do2D:
        plotter.make_overlays_2D(logScale = log, doProjection = False)

