#ifndef _HistManager_H
#define _HistManager_H

// c++ libraries
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

// ROOT libraries
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"

using namespace std;

class HistManager : public TObject {

    public:

        // constructor
        HistManager();

        // destructor
        ~HistManager();

        virtual string  str(int i) {return static_cast<ostringstream*>( &(ostringstream() << i) )->str();}

        void SetWeight(float);
        void SetWeightError(float);
        void SetFileNumber(unsigned);
        void SetDirectory(string);
        void AddFile(TFile*);

        TFile* GetFile();

        // fill 1D histogram 
        void Fill1DHist(float x, 
                std::string name, std::string title,
                int bins, float xmin, float xmax);

        void Fill1DHistUnevenBins(float x, 
                std::string name, std::string title,
                int bins, float *binEdges);


        // fill 2D histogram
        void Fill2DHist(float x, float y, 
                std::string name, std::string title,
                int binsx, float xmin, float xmax,
                int binsy, float ymin, float ymax);

        void Fill2DHistUnevenBins(float x, float y, std::string name, std::string title,
                int binsx, float *binEdgesx,
                int binsy, float *binEdgesy);


        // make a profile histogram
        void FillProfile(float x, float y, std::string name, std::string title,
                int binsx, float xmin, float xmax,
                float ymin, float ymax);
        void FillProfile(float x, float y, std::string name, std::string title,
                int binsx, float* bins);

        // methods for manipulating cached histograms
        void SetBinNames1D(string, string*, unsigned);
        void SetBinNames2D(string, string*, unsigned, string*, unsigned);

        ClassDef(HistManager, 1);

    private:

        vector<TFile*> _files;
        unsigned    _fileNumber;
        string      _directory;  
        float       _evtWeight;
        float       _weightErr;

        // map to hold histograms
        vector<map<string,TH1F*> >      histMap1D;
        vector<map<string,TH2F*> >      histMap2D;
        vector<map<string,TProfile*> >  histMapProf;

};

#endif


#if !defined(__CINT__)
ClassImp(HistManager);
#endif
