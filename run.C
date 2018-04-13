void run(TString filename)
{
    gROOT->Reset();
    gROOT->ProcessLine(".x loadClasses.C" );
    gErrorIgnoreLevel=2001;

    WCPlot *p = new WCPlot(filename);
    p->PlotInBeamLight();
    p->PlotInBeamCharge2D();
}