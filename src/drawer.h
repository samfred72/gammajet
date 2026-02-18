#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include <string>
#include <vector>
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
using namespace std;

class drawer {
  public :
    drawer() {
      TFile * f = TFile::Open(         "/home/samson72/sphnx/gammajet/hists/histsData.root");         
      TFile * f05_p = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Photon5"));
      TFile * f10_p = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Photon10"));
      TFile * f20_p = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Photon20"));
      TFile * f05_j = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Jet5"));
      TFile * f10_j = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Jet10"));
      TFile * f20_j = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Jet20"));
      TFile * f30_j = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Jet30"));
      TFile * f50_j = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Jet50"));
      TFile * f70_j = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists%s.root","Jet70"));
      dfiles[0] = f;

      pfiles[0] = f05_p;
      pfiles[1] = f10_p;
      pfiles[2] = f20_p;

      jfiles[0] = f05_j;
      jfiles[1] = f10_j;
      jfiles[2] = f20_j;
      jfiles[3] = f30_j;
      jfiles[4] = f50_j;
      jfiles[5] = f70_j;
    }
    ~drawer(); 

    void drawLine(float x1, float y1, float x2, float y2);
    void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=18);
    void drawAll(vector<string> samples, vector<string> features, float drawx, float drawy, int fontsize, float csize);
    void drawMany(vector<string> features, float drawx, float drawy, int fontsize, float csize);
    void scale(TH1D * h, float low = 0, float high = 2);
    void format(TH1D * h, int type);
    void format(TF1 * h, int type);
    TF1 * fit(TH1D * h, float low, float high, const char * options);
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> collect_hists(const char * histname, int type);
    TH1D * combine_hists(TH1D * A, TH1D * B, TH1D * C, TH1D * D, int ipt, string name); 
    TH1D * combineMC(const char * histname, bool isphoton); 
    TH2D * combineMC2d(const char * histname, bool isphoton); 
    TH1D * get(const char * histname, int type);
    TH2D * get2d(const char * histname, int type);

  private :

    const static int ndfiles = 1;
    const static int npfiles = 3;
    const static int njfiles = 6;
    vector<TFile*> dfiles{ndfiles};
    vector<TFile*> pfiles{npfiles};
    vector<TFile*> jfiles{njfiles};

    vector<int> psamples = {5,10,20};
    vector<int> jsamples = {5,10,20,30,50,70};

    map<bool,map<int,double>> scalemap = {
      {0,{{5,1.369e+08},{10,3.997e+06},{15,4.073e+05},{20,6.218e+04},{30,2.502e+03},{50,7.2695},{70,1.034e-02}}},
      {1,{{5,146359.3},{10,6944.675},{20,130.4461}}}
    };
};
