#include "HistManager.h"

HistManager::HistManager()
{
    _evtWeight  = 1.;
    _weightErr  = 0.;
    _fileNumber = 0;
    _directory  = "";
}


HistManager::~HistManager() 
{
}

void HistManager::SetWeight(float w)
{ 
    _evtWeight = w;
}

void HistManager::SetWeightError(float e)
{ 
    _weightErr = e;
}

void HistManager::SetFileNumber(unsigned f)
{
    _fileNumber = f;
}

void HistManager::SetDirectory(string d)
{
    _directory = d;
}

void HistManager::AddFile(TFile* f)
{
    _files.push_back(f);

    map<string,TH1F*>       h1Map;
    map<string,TH2F*>       h2Map;
    map<string,TProfile*>   p1Map;

    histMap1D.push_back(h1Map);
    histMap2D.push_back(h2Map);
    histMapProf.push_back(p1Map);
}

TFile* HistManager::GetFile()
{
    return _files[_fileNumber];
}

void HistManager::Fill1DHist(float x,
        string name, string title,
        int bins, float xmin, float xmax)
{

    string key = _directory + "/" + name;
    map<string,TH1F*>::iterator it = histMap1D[_fileNumber].find(key);

    if (it == histMap1D[_fileNumber].end()){
        _files[_fileNumber]->cd(_directory.c_str());
        histMap1D[_fileNumber][key] = new TH1F(name.c_str(),title.c_str(),bins,xmin,xmax);
        _files[_fileNumber]->cd();
    }

    histMap1D[_fileNumber][key]->Fill(x, _evtWeight);

    unsigned histBin = histMap1D[_fileNumber][key]->FindBin(x);
    float errStat = histMap1D[_fileNumber][key]->GetBinError(histBin);
    float errComb = sqrt(pow(errStat, 2) + pow(_weightErr, 2));
    histMap1D[_fileNumber][key]->SetBinError(histBin, errComb);

}

void HistManager::Fill1DHistUnevenBins(float x, 
        string name, string title,
        int bins, float *binEdges)
{
    string key = _directory + "/" + name;
    map<string,TH1F*>::iterator it = histMap1D[_fileNumber].find(key);

    if (it == histMap1D[_fileNumber].end()){
        _files[_fileNumber]->cd(_directory.c_str());
        histMap1D[_fileNumber][key] = new TH1F(name.c_str(),title.c_str(),bins,binEdges);
        _files[_fileNumber]->cd();
    }

    histMap1D[_fileNumber][key]->Fill(x, _evtWeight);
}

void HistManager::Fill2DHist(float x, float y, 
        string name, string title,
        int binsx, float xmin, float xmax,
        int binsy, float ymin, float ymax)
{

    string key = _directory + "/" + name;
    map<string,TH2F*>::iterator it = histMap2D[_fileNumber].find(key);

    if (it == histMap2D[_fileNumber].end()){
        _files[_fileNumber]->cd(_directory.c_str());
        histMap2D[_fileNumber][key] = new TH2F(name.c_str(),title.c_str(),binsx,xmin,xmax,binsy,ymin,ymax);
        _files[_fileNumber]->cd();
    }

    histMap2D[_fileNumber][key]->Fill(x, y, _evtWeight);
}

void HistManager::Fill2DHistUnevenBins(float x, float y, 
        string name, string title,
        int binsx, float *binEdgesx,
        int binsy, float *binEdgesy)
{

    string key = _directory + "/" + name;
    map<string,TH2F*>::iterator it = histMap2D[_fileNumber].find(key);

    if (it == histMap2D[_fileNumber].end()){
        _files[_fileNumber]->cd(_directory.c_str());
        histMap2D[_fileNumber][key] = new TH2F(name.c_str(),title.c_str(),binsx,binEdgesx,binsy,binEdgesy);
        _files[_fileNumber]->cd();
    }

    histMap2D[_fileNumber][key]->Fill(x, y, _evtWeight);
}


void HistManager::FillProfile(float x, float y, 
        string name, string title,
        int binsx, float xmin, float xmax,
        float ymin, float ymax)
{

    string key = _directory + "/" + name;
    map<string,TProfile*>::iterator it = histMapProf[_fileNumber].find(key);

    if (it == histMapProf[_fileNumber].end()){
        _files[_fileNumber]->cd(_directory.c_str());
        histMapProf[_fileNumber][key] = new TProfile(name.c_str(),title.c_str(),binsx,xmin,xmax,ymin,ymax);
        _files[_fileNumber]->cd();
    }

    histMapProf[_fileNumber][key]->Fill(x, y, _evtWeight);
}


void HistManager::FillProfile(float x, float y, 
        string name, string title,
        int binsx, float* bins)
{

    string key = _directory + "/" + name;
    map<string,TProfile*>::iterator it = histMapProf[_fileNumber].find(key);

    if (it == histMapProf[_fileNumber].end()){
        _files[_fileNumber]->cd(_directory.c_str());
        histMapProf[_fileNumber][key] = new TProfile(name.c_str(),title.c_str(), binsx, bins);
        _files[_fileNumber]->cd();
    }

    histMapProf[_fileNumber][key]->Fill(x, y, _evtWeight);
}


void HistManager::SetBinNames1D(string histName, string* binNames, unsigned nBins)
{

    for (unsigned i = 0; i < _files.size(); ++i) {
        this->SetFileNumber(i);

        string key = _directory + "/" + histName;
        map<string,TH1F*>::iterator it = histMap1D[i].find(key);

        if (it != histMap1D[i].end()){
            cout << ".";
            for (unsigned j = 0; j < nBins; ++j) histMap1D[i][key]->GetXaxis()->SetBinLabel(j+1, binNames[j].c_str());
            histMap1D[i][key]->GetXaxis()->SetLabelSize(0.06);
        }
    }
}

void HistManager::SetBinNames2D(string histName, string* binNamesX, unsigned nBinsX, string* binNamesY, unsigned nBinsY)
{
    for (unsigned i = 0; i < _files.size(); ++i) {
        this->SetFileNumber(i);

        string key = _directory + "/" + histName;
        map<string,TH2F*>::iterator it = histMap2D[i].find(key);

        if (it != histMap2D[i].end()){
            for (unsigned j = 0; j < nBinsX; ++j) 
                histMap2D[i][key]->GetXaxis()->SetBinLabel(j+1, binNamesX[j].c_str());
            for (unsigned j = 0; j < nBinsY; ++j) 
                histMap2D[i][key]->GetYaxis()->SetBinLabel(j+1, binNamesY[j].c_str());
        }
    }
}
