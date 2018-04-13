#include "WCPlot.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TUnixSystem.h"
#include "TGaxis.h"
#include "TMath.h"
// #include "TEntryList.h"
// #include "TDirectory.h"

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

WCPlot::WCPlot(){}

//----------------------------------------------------------------
WCPlot::WCPlot(const char* filename)
{
    InitStyle();

    rootFile = new TFile(filename);
    ReadRunInfo();

    TString cmd;
    plotDir = "tmp_plots";
    cmd = TString::Format("mkdir -p %s", plotDir.Data());
    gSystem->Exec(cmd.Data());
}

//----------------------------------------------------------------
WCPlot::~WCPlot()
{
    rootFile->Close();
    delete rootFile;
}

//----------------------------------------------------------------
void WCPlot::ReadRunInfo()
{
    TTree *t = (TTree*)rootFile->Get("Trun");
    if (!t) {
        cout << "root file doesn't have Trun info" << endl;
        exit(1);
    }
    t->SetBranchAddress("runNo", &runNo);
    t->SetBranchAddress("subRunNo", &subRunNo);
    t->SetBranchAddress("eventNo", &eventNo);
    t->GetEntry(0);
    eventIdStr = TString::Format("%d-%d-%d", runNo, subRunNo, eventNo);
    cout << "processing event " << eventIdStr << endl;
}

//----------------------------------------------------------------
void WCPlot::PlotInBeamLight()
{
    TCanvas *c1 = 0;
    const float CONV = 0.09375;
    const float timemax = 10; // us to show
    TLine *llow  = new TLine(3/CONV, 0, 3/CONV, 32);
    TLine *lhigh = new TLine(5/CONV, 0, 5/CONV, 32);
    llow->SetLineWidth(1);
    lhigh->SetLineWidth(1);
    llow->SetLineColor(kWhite);
    lhigh->SetLineColor(kWhite);
    llow->SetLineStyle(2);
    lhigh->SetLineStyle(2);

    TH2F *hdecon = (TH2F*)rootFile->Get("hdecon");
    c1 = new TCanvas("c1","c1",800,600);
    hdecon->SetTitle("Deconvoluted In-beam Light");
    hdecon->GetXaxis()->SetTitle("Time [#mus]");
    hdecon->GetYaxis()->SetTitle("PMT #");
    hdecon->GetZaxis()->SetRangeUser(-1, 25);
    hdecon->GetXaxis()->SetRangeUser(0, timemax/CONV);
    hdecon->GetXaxis()->SetLabelOffset(999);
    hdecon->GetXaxis()->SetTickLength(0);
    hdecon->Draw("colz");
    llow->Draw();
    lhigh->Draw();
    TGaxis *axis1 = new TGaxis(0,0,timemax/CONV,0,0,timemax,505,"");
    // axis1->SetLabelOffset(-0.04);
    axis1->SetLabelFont(42);
    axis1->SetLabelSize(0.04);
    axis1->Draw();
    SaveAs(c1, "inBeamLight2D_decon");
    delete c1;

    TH2F *hl1 = (TH2F*)rootFile->Get("hl1");
    c1 = new TCanvas("c1","c1",800,600);
    hl1->SetTitle("L1-fitted In-beam Light");
    hl1->GetXaxis()->SetTitle("Time [#mus]");
    hl1->GetYaxis()->SetTitle("PMT #");
    hl1->GetZaxis()->SetRangeUser(-1, 25);
    hl1->GetXaxis()->SetRangeUser(0, timemax/CONV);
    hl1->GetXaxis()->SetLabelOffset(999);
    hl1->GetXaxis()->SetTickLength(0);
    hl1->Draw("colz");
    llow->Draw();
    lhigh->Draw();
    axis1->Draw();
    SaveAs(c1, "inBeamLight2D_l1");
    delete c1;


}

//----------------------------------------------------------------
void WCPlot::PlotInBeamCharge2D()
{
    TTree *T_flash = (TTree*)rootFile->Get("T_flash");
    int nFlash = T_flash->Draw("flash_id","time<6 && time>2", "goff");
    cout << "found " << nFlash << " in beam flashes" << endl;
    if (nFlash==0) { return;}

    TString sel = TString::Format("flash_id==%i", int(T_flash->GetV1()[0]));
    for (int i=1; i<nFlash; i++) {
        sel += TString::Format(" || flash_id==%i", int(T_flash->GetV1()[i]));
    }
    cout << sel << endl;
    TTree *T_match = (TTree*)rootFile->Get("T_match");
    int nMCluster = T_match->Draw("tpc_cluster_id", sel.Data());
    cout << "found " << nMCluster << " in-beam clusters" << endl;

    if (nMCluster==0) { return;}

    TTree *T = (TTree*)rootFile->Get("T_proj");
    TCanvas *c1 = 0;

    vector<int> *cluster_id = new std::vector<int>;
    vector<vector<int> > *channel = new std::vector<std::vector<int> >;
    vector<vector<int> > *time_slice =new std::vector<std::vector<int> >;
    vector<vector<int> > *charge =new std::vector<std::vector<int> >;

    T->SetBranchAddress("cluster_id", &cluster_id);
    T->SetBranchAddress("channel", &channel);
    T->SetBranchAddress("time_slice", &time_slice);
    T->SetBranchAddress("charge", &charge);
    T->GetEntry(0);
    int nClusters = cluster_id->size();
    const int MARGIN = 10;
    const int SCALE = 1.;
    const int PIXEL_X = 3*SCALE;
    const int PIXEL_Y = 2*1.1*SCALE;
    const float ZMAX = 20.0;

    for (int mcluster=0; mcluster<nMCluster; mcluster++) {
        int cluster_inbeam_id = T_match->GetV1()[mcluster];
        // if (cluster_inbeam_id > 20) { continue; }
        int cluster_ind = -1;
        for (int i=0; i< nClusters; i++) {
            if (cluster_id->at(i) == cluster_inbeam_id) { cluster_ind = i; break; }
        }
        if (cluster_inbeam_id == -1 ) { continue; }

        vector<int> &time_slice_inbeam = time_slice->at(cluster_ind);
        vector<int> &channel_inbeam = channel->at(cluster_ind);
        vector<int> &charge_inbeam = charge->at(cluster_ind);
        int n = charge_inbeam.size();

        int t_min = 10000;
        int t_max = 0;
        int u_min = 2400;
        int u_max = 0;
        int v_min = 4800;
        int v_max = 2400;
        int w_min = 8256;
        int w_max = 4800;

        for (int i=0; i<n; i++) {
            int t = time_slice_inbeam[i];
            int ch = channel_inbeam[i];
            int q = charge_inbeam[i];
            // if (q/10000<0.1) continue;

            if (t>t_max) { t_max = t; }
            if (t<t_min) {t_min = t; }

            if (ch<2400 && ch<u_min) {u_min = ch; }
            if (ch<2400 && ch>u_max) {u_max = ch;}

            if (ch<4800 && ch>=2400 && ch<v_min) {v_min = ch; }
            if (ch<4800 && ch>=2400 && ch>v_max) {v_max = ch; }

            if (ch<8256 && ch>=4800 && ch<w_min) {w_min = ch; }
            if (ch<8256 && ch>=4800 && ch>w_max) {w_max = ch; }
        }

        if (t_min >= t_max) { continue; }

        // cout << t_min << ", " << t_max << endl;
        // cout << u_min << ", " << u_max << endl;
        // cout << v_min << ", " << v_max << endl;
        // cout << w_min << ", " << w_max << endl;

        // return;

        TH2F *hw = new TH2F("hw","hw",w_max-w_min+2*MARGIN, w_min-MARGIN,w_max+MARGIN, t_max-t_min+2*MARGIN,t_min-MARGIN,t_max+MARGIN);
        TH2F *hu = new TH2F("hu","hu",u_max-u_min+2*MARGIN, u_min-MARGIN,u_max+MARGIN, t_max-t_min+2*MARGIN,t_min-MARGIN,t_max+MARGIN);
        TH2F *hv = new TH2F("hv","hv",v_max-v_min+2*MARGIN, v_min-MARGIN,v_max+MARGIN, t_max-t_min+2*MARGIN,t_min-MARGIN,t_max+MARGIN);
        for (int i=0; i<n; i++) {
            int t = time_slice_inbeam[i];
            int ch = channel_inbeam[i];
            float q = charge_inbeam[i] / 1000.;
            hw->Fill(ch, t, q);
            hu->Fill(ch, t, q);
            hv->Fill(ch, t, q);
        }
        c1 = new TCanvas("c1","c1",(w_max-w_min+2*MARGIN)*PIXEL_X+200,(t_max-t_min+2*MARGIN)*PIXEL_Y+200);
        hw->SetTitle(TString::Format("Y Plane, Cluster %i", cluster_inbeam_id));
        hw->GetXaxis()->SetTitle("channel #");
        hw->GetYaxis()->SetTitle("time slice (x2 #mus)");
        hw->GetZaxis()->SetRangeUser(-0.1, ZMAX);
        hw->Draw("colz");
        SetAspectRatio(c1, 1);
        SaveAs(c1, TString::Format("c%i_charge2D_w", cluster_inbeam_id));
        delete c1;

        c1 = new TCanvas("c1","c1",(u_max-u_min+2*MARGIN)*PIXEL_X+200,(t_max-t_min+2*MARGIN)*PIXEL_Y+200);
        hu->SetTitle(TString::Format("U Plane, Cluster %i", cluster_inbeam_id));
        hu->GetXaxis()->SetTitle("channel #");
        hu->GetYaxis()->SetTitle("time slice (x2 #mus)");
        hu->GetZaxis()->SetRangeUser(-0.1, ZMAX);
        hu->Draw("colz");
        SetAspectRatio(c1, 1);
        SaveAs(c1, TString::Format("c%i_charge2D_u", cluster_inbeam_id));
        delete c1;

        c1 = new TCanvas("c1","c1",(v_max-v_min+2*MARGIN)*PIXEL_X+200,(t_max-t_min+2*MARGIN)*PIXEL_Y+200);
        hv->SetTitle(TString::Format("V Plane, Cluster %i", cluster_inbeam_id));
        hv->GetXaxis()->SetTitle("channel #");
        hv->GetYaxis()->SetTitle("time slice (x2 #mus)");
        hv->GetZaxis()->SetRangeUser(-0.1, ZMAX);
        hv->Draw("colz");
        SetAspectRatio(c1, 1);
        SaveAs(c1, TString::Format("c%i_charge2D_v", cluster_inbeam_id));
        delete c1;

    }

    // TH2F *hw = new TH2F("hw","hw",8400,0,8400,2400,0,2400);
    // T->Project("hw","time_slice:channel","charge*(cluster_id==19)");

}

//----------------------------------------------------------------
void WCPlot::SaveAs(TCanvas *c1, const char* name, const char* ext)
{
    TString filename = TString::Format("%s/%s_%s.%s",
        plotDir.Data(), name, eventIdStr.Data(), ext);
    c1->SaveAs(filename);
}

//----------------------------------------------------------------
void WCPlot::PlotAll()
{
    PlotInBeamLight();
    PlotInBeamCharge2D();
}

//----------------------------------------------------------------
void WCPlot::InitStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTitleFontSize(0.045);
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetHistLineWidth(2);
    gStyle->SetLegendBorderSize(0);

    // gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.15);

    gStyle->SetMarkerSize(0.3);

    PaletteRainbow();

    gROOT->ForceStyle();
}

//----------------------------------------------------------------
void WCPlot::PaletteRainbow()
{
    gStyle->SetFrameFillColor(TColor::GetColor(float(0.1), float(0.1), float(0.1)));

    // http://diana.parno.net/thoughts/?p=28
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

}


bool WCPlot::SetAspectRatio(TCanvas* const c, const int axis)
{
    if(!c){
        cout << "Error in SetRealAspectRatio: canvas is NULL" << endl;
        return false;
    }

    {
        //Get the current min-max values if SetUserRange has been called
        c->Update();
        const Double_t xmin = c->GetUxmin();
        const Double_t xmax = c->GetUxmax();
        const Double_t ymin = c->GetUymin();
        const Double_t ymax = c->GetUymax();

        //Get the length of zoomed x and y axes
        const Double_t xlength = xmax - xmin;
        const Double_t ylength = ymax - ymin;
        const Double_t ratio = xlength/ylength;

        //Get how many pixels are occupied by the canvas
        const Int_t npx = c->GetWw();
        const Int_t npy = c->GetWh();

        //Get x-y coordinates at the edges of the canvas (extrapolating outside the axes, NOT at the edges of the histogram)
        const Double_t x1 = c->GetX1();
        const Double_t y1 = c->GetY1();
        const Double_t x2 = c->GetX2();
        const Double_t y2 = c->GetY2();

        //Get the length of extrapolated x and y axes
        const Double_t xlength2 = x2 - x1;
        const Double_t ylength2 = y2 - y1;
        const Double_t ratio2 = xlength2/ylength2;

        //Get now number of pixels including window borders
        const Int_t bnpx = c->GetWindowWidth();
        const Int_t bnpy = c->GetWindowHeight();

        // cout << "WindX\tWindY\tCanvX\tCanvY\tx1\ty1\tx2\ty2\tratiox/y\tCanvX/CanvY" << endl;
        // cout << bnpx << "\t" << bnpy << "\t" << npx << "\t" << npy << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << ratio2 << "\t" << (double)npx/npy << "\tOriginal Canvas" << endl;

        if(axis==1)
        {
            c->SetCanvasSize(TMath::Nint(npy*ratio2), npy);
            c->SetWindowSize((bnpx-npx)+TMath::Nint(npy*ratio2), bnpy);
        }
        else if(axis==2)
        {
            c->SetCanvasSize(npx, TMath::Nint(npx/ratio2));
            c->SetWindowSize(bnpx, (bnpy-npy)+TMath::Nint(npx/ratio2));
        }
        else
        {
            cout << "Error in SetRealAspectRatio: axis value " << axis << " is neither 1 (resize along x-axis) nor 2 (resize along y-axis).";
            return false;
        }
    }
}