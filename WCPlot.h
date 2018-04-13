#ifndef WCReader_H
#define WCReader_H
#pragma link C++ class vector<vector<int> >+;

#include "TString.h"
class TFile;
class TCanvas;

class WCPlot {
public:
    TFile *rootFile;
    int runNo;
    int subRunNo;
    int eventNo;
    TString eventIdStr;
    TString plotDir;

    WCPlot();
    WCPlot(const char* filename);
    virtual ~WCPlot();

    void ReadRunInfo();

    // plots
    void PlotInBeamLight();
    void PlotInBeamCharge2D();

    // utilities
    void SaveAs(TCanvas *c1, const char* name, const char* ext="png");
    void PlotAll();
    bool SetAspectRatio(TCanvas* const c, const int axis=1);


private:
    void InitStyle();
    void PaletteRainbow();
};


#endif
