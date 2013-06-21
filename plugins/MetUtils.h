/*
    Utilities for retrieving weights for PU,etc.
*/

#ifndef _MetUtils_H
#define _MetUtils_H

#include <string>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

using namespace std;

class MetUtils: public TObject {
    public:
        MetUtils() {};
        virtual ~MetUtils() {};
        MetUtils(string selection, string dataPeriod);

        ClassDef(MetUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _selection;

};

#endif

#if !defined(__CINT__)
ClassImp(MetUtils);
#endif
