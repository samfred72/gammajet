#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

drawer d;

void draw_many(TCanvas * c, const char * cname, const char * info1, const char * info2, const char * variation, 
    TH1D * H1[], TH1D * H2[],int ir = 2) {
  c->cd();
  gStyle->SetStripDecimals(0);

  float drawx = 0.75;
  float drawy = 0.88;
  float fontsize = 60;
  int offset = 0;
  float scale = 1.15;

  TPad * p = new TPad("p","",0,0,1,1);
  p->SetLeftMargin(0.25);
  p->SetBottomMargin(0.25);
  int nx = 3;
  int ny = 2;
  p->Divide(nx+1,ny,0,0);
  p->Draw();

  for (int ipt = 0; ipt < ana::nTrijetPtBins; ipt++) {
    TH1D * h1 = H1[ipt];
    TH1D * h2 = H2[ipt];

    int index = ipt + 1 + offset;
    if (index % (nx+1) == 0) {index++; offset++;}
    p->cd(index);
    gPad->SetRightMargin(0.01);
    if ((ipt+1) % 3 == 0) gPad->SetRightMargin(0.057);
    gPad->SetTicks(1,1);
    if ((ipt) / nx == 0) gPad->SetBottomMargin(0.04);
    gPad->SetTopMargin(0.01);
    if ((ipt) % nx != 0) gPad->SetLeftMargin(0.056);
    if (index == (ny*(nx+1) - 1)) h1->GetXaxis()->SetTitle("x_{J}^{-1}"); // bottom right plot
    else h1->GetXaxis()->SetTitle("");
    h1->GetXaxis()->SetTitleSize(0.10*scale);
    h1->GetXaxis()->SetLabelSize(0.072*scale);
    h1->GetXaxis()->SetLabelOffset(0.018);
    if (index == (nx+1)*(ny-1)+1) h1->GetXaxis()->SetLabelSize(0.062*scale); // bottom left plot
    if (index == (nx+1)*(ny-1)+1) h1->GetXaxis()->SetLabelOffset(0.026); // bottom left plot
    if (ipt < 3 ) h1->GetXaxis()->SetLabelSize(0.00);
    h1->GetYaxis()->SetNdivisions(605);
    //if (index == (nx+1)*(ny-1)+1) h1->GetYaxis()->ChangeLabel(-1,-1,0);
    //if (index != (ny*(nx+1) - 1)) h1->GetXaxis()->ChangeLabel(-1,-1,0);
    if (index == 1) h1->GetYaxis()->SetTitle("Normalized Counts"); // top left plot
    else h1->GetYaxis()->SetTitle("");
    h1->GetYaxis()->SetTitleSize(0.10);
    h1->GetYaxis()->SetLabelSize(0.084*scale);
    if (ipt != 0 && ipt != 3 ) h1->GetYaxis()->SetLabelSize(0.00);
    if (ipt == 3) h1->GetYaxis()->SetLabelSize(0.064*scale); // bottom left plot
    h1->GetYaxis()->SetLabelOffset(0.04);
    if (ipt == 5) {
      h1->Scale(0.3);
      h2->Scale(0.3);
    } 
    h1->GetYaxis()->SetMaxDigits(3);
    h1->GetYaxis()->SetDecimals(3);

    h1->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0.0");
    h1->GetYaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "1.0");
    h1->GetYaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "2.0");
    h1->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "0.5");
    h1->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "1.0");
    h1->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "1.5");
    h1->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "2.0");
    h1->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "2.5");
    h1->GetYaxis()->SetRangeUser(0,2.70);

    h1->Draw("hist e1");
    h2->Draw("hist same");
    
    float drawx = 0.15;
    float drawy_temp = 0.88;
    if ((index - 1) % (nx+1) == 0) drawx = 0.35; // on the left
    if (ipt < nx) drawy_temp = 0.83;
    d.drawMany({
        Form("#bf{%.0f GeV < p_{T,1} < %.0f GeV}",ana::trijetPtBins[ipt],ana::trijetPtBins[ipt+1]),
        //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}","Data",  h1->GetMean(), h1->GetMeanError()),
        //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}","MC",  h2->GetMean(), h2->GetMeanError())
        //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info1,  h1->GetMean(), h1->GetMeanError()),
        //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info2,  h2->GetMean(), h2->GetMeanError())
        },drawx,drawy_temp,42,c->GetWh()/3.0);
    if (ipt == 5) {
      d.drawMany({
          "Scaled",
          "by 0.3",
          },0.64,drawy_temp-0.13,42,c->GetWh()/2.0);
    }
    float err = TMath::Sqrt(h1->GetMeanError()*h1->GetMeanError()+h2->GetMeanError()*h2->GetMeanError());
    cout << ana::trijetPtBins[ipt] << "-" << ana::trijetPtBins[ipt+1] << std::fixed << std::setprecision(3) << " \\GeV & " 
      //<< h1->GetMean() << " $\\pm$ " << h1->GetMeanError() << " & " 
      //<< h2->GetMean() << " $\\pm$ " << h2->GetMeanError() << " & "
      << h1->GetMean()/h2->GetMean() << " $\\pm$ " << err << " \\\\ " << std::defaultfloat << endl;
      
    
    TLine * mline1 = new TLine(h1->GetMean(),0,h1->GetMean(),h1->GetMaximum()*0.8);
    mline1->SetLineColor(h1->GetLineColor());
    mline1->SetLineStyle(8);
    TLine * mline2 = new TLine(h2->GetMean(),0,h2->GetMean(),h1->GetMaximum()*0.8);
    mline2->SetLineColor(h2->GetLineColor());
    mline2->SetLineStyle(7);
    mline1->Draw();
    mline2->Draw();
  }

  p->cd();
  d.drawAll({
      info1, 
      //info2
      },
      {
      //"Analysis cuts",
      Form("Jet R=%1.1f",ana::JetRs[ir]),
      Form("|#eta| < %1.1f", ana::etacut - ana::JetRs[ir]),
      //variation,
      },
      drawx,drawy,fontsize,c->GetWh());
  TLegend * l2 = new TLegend(drawx,drawy-.43,0.99,drawy-.23);
  //TLegend * l2 = new TLegend(drawx,drawy-.20,0.99,drawy-.10);
  l2->SetLineWidth(0);
  l2->AddEntry(H1[0],"Data");
  l2->AddEntry(H2[0],info2);
    
  TLine * mline1 = new TLine(0,0,1,1);
  mline1->SetLineColor(H1[0]->GetLineColor());
  mline1->SetLineStyle(8);
  TLine * mline2 = new TLine(0,0,1,1);
  mline2->SetLineColor(H2[0]->GetLineColor());
  mline2->SetLineStyle(7);
  
  l2->AddEntry(mline1,Form("Mean %s","Data"),"l");
  l2->AddEntry(mline2,Form("Mean %s",info2),"l");
  
  l2->SetLineWidth(0);
  l2->Draw();
  c->SaveAs(cname);
  return;
}

void draw_one(TCanvas * c, const char * cname, const char * info1, const char * info2,
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H1,          
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H2,          
    vector<vector<vector<vector<vector<vector<TF1*>>>>>> f1,          
    vector<vector<vector<vector<vector<vector<TF1*>>>>>> f2,          
    vector<vector<vector<vector<vector<vector<float>>>>>> m1,          
    vector<vector<vector<vector<vector<vector<float>>>>>> m2,          
    int ipt1, int ir1, int icalib1, int ibdt1, int i3jet1, int iabcd1,
    int ipt2, int ir2, int icalib2, int ibdt2, int i3jet2, int iabcd2) {
  float drawx = 0.6;
  float drawy = 0.85;
  string calibstring = (icalib1 == 1 ? "JES Calibrated" : "Uncalibrated Jets");
  vector<string> t1 = {"#bf{Analysis region A}","#bf{Analysis region B}","#bf{Analysis region C}","#bf{Analysis region D}", "Combined ABCD"};
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.15);
  
  TH1D * h1 = H1[ipt1][ir1][icalib1][ibdt1][i3jet1][iabcd1];
  TH1D * h2 = H2[ipt2][ir2][icalib2][ibdt2][i3jet2][iabcd2];
  h1->GetXaxis()->SetTitle("p_{T,max}^{Jet}/p_{T,max}^{cluster}");
  h1->GetXaxis()->SetTitleSize(0.04);
  h1->GetXaxis()->SetTitleOffset(1.5);
  h1->GetXaxis()->SetLabelSize(0.04);

  h1->GetYaxis()->SetTitle("Arbitrary Units");
  h1->GetYaxis()->SetTitleSize(0.04);
  h1->GetYaxis()->SetTitleOffset(2);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetMaxDigits(3);
  h1->GetYaxis()->SetDecimals(2);
  h1->Draw("hist e");
  h2->Draw("hist e same");
  //h3->Draw("hist e same");

  float lowcluster = ana::ptBins[ipt1];
  float highcluster = ana::ptBins[ipt1+1];
  float minjet = (icalib1 ? ana::jet_calib_pt_cut[ir1] : ana::jet_pt_cut[ir1]);
  TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,h1->GetMaximum());
  line->SetLineStyle(8);
  line->Draw();
  TLine * line1 = new TLine(m1[ipt1][ir1][icalib1][ibdt1][i3jet1][iabcd1],0,m1[ipt1][ir1][icalib1][ibdt1][i3jet1][iabcd1],h1->GetMaximum());
  line1->SetLineStyle(9);
  line1->SetLineColor(kBlue);
  line1->Draw();
  TLine * line2 = new TLine(m2[ipt2][ir2][icalib2][ibdt2][i3jet2][iabcd2],0,m2[ipt2][ir2][icalib2][ibdt2][i3jet2][iabcd2],h1->GetMaximum());
  line2->SetLineStyle(9);
  line2->SetLineColor(kSpring-1);
  line2->Draw();
  TLegend * l = new TLegend(.25,.75,.55,.87);
  l->SetLineWidth(0);
  l->SetFillStyle(0);
  l->AddEntry(h1,("data"));
  l->AddEntry(h2,("MC Photon reco"));
  //l->AddEntry(h3,("MC Jet reco"));
  l->Draw();
  if (iabcd1 == 0 || iabcd1 == 4) {
    f1[ipt1][ir1][icalib1][ibdt1][i3jet1][iabcd1]->Draw("same");
    f2[ipt2][ir2][icalib2][ibdt2][i3jet2][iabcd2]->Draw("same");
    d.drawAll({
        "run 47289-53864",
        //"MC run28 Jet",
        "MC run28 Photon"},
        {Form("%0.0f GeV < p_{T}^{cluster} < %0.0f GeV",ana::ptBins[ipt1],ana::ptBins[ipt1+1]),
        "p_{T}^{jet} > 3 GeV", 
        Form("Jet R=%0.1f",ana::JetRs[ir1]), 
        "Analysis cuts",
        calibstring,
        t1[iabcd1].c_str(),
        Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info1,f1[ipt1][ir1][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1),f1[ipt1][ir1][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1)),
        Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info2,f2[ipt2][ir2][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1),f2[ipt2][ir2][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1))
        },
        drawx,drawy,15,c->GetWh());
  }
  else {
    d.drawAll({
        "run 47289-53864",
        //"MC run28 Jet",
        "MC run28 Photon"},
        {Form("%0.0f GeV < p_{T}^{cluster} < %0.0f GeV",ana::ptBins[ipt1],ana::ptBins[ipt1+1]),
        "p_{T}^{jet} > 3 GeV", 
        Form("Jet R=%0.1f",ana::JetRs[ir1]), 
        "Analysis cuts",
        calibstring,
        t1[iabcd1].c_str()
        },
        drawx,drawy,15,c->GetWh());
  }
  return;
}

void draw_multijet() {
  gStyle->SetOptStat(0);
  if (!gROOT->IsBatch()) {
    cout << "Run with -b flag or else!!" << endl;
    gROOT->SetBatch(kTRUE);
    //return;
  }
  
  cout << "collecting/formatting hists" << endl;
  const int ntypes = 5;
  const char * typenames[ntypes] = {"","_JERhigh","_JERlow","_JERreco","_HERWIG"};
  TH1D * hd[ana::nTrijetPtBins][ana::nJetR][ntypes]; // ntypes is for nominal, JERhigh, JERlow, JERreco, HERWIG
  TH1D * hj[ana::nTrijetPtBins][ana::nJetR][ntypes]; // ntypes is for nominal, JERhigh, JERlow, JERreco, HERWIG
  for (int k = 0; k < ntypes; k++) {
    for (int j = 0; j < ana::nJetR; j++) {
      TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet%s_%s.root",typenames[k],ana::rnames[j]),"READ");
      for (int i = 0; i < ana::nTrijetPtBins; i++) {
        hd[i][j][k] = (TH1D*)f->Get(Form("hdata%i",i));
        hj[i][j][k] = (TH1D*)f->Get(Form("htrijet%i",i));
        
        hd[i][j][k]->SetLineColor(kBlue);
        hd[i][j][k]->Scale(1.0/hd[i][j][k]->Integral("width"));
        hd[i][j][k]->GetYaxis()->SetRangeUser(0,hd[i][j][k]->GetMaximum()*1.5);
        hd[i][j][k]->SetLineWidth(2);
        
        hj[i][j][k]->SetLineColor(kMagenta);
        if (strcmp(typenames[k],"_HERWIG")==0) hj[i][j][k]->SetLineColor(kSpring-1);
        hj[i][j][k]->Scale(1.0/hj[i][j][k]->Integral("width"));
        hj[i][j][k]->GetYaxis()->SetRangeUser(0,hj[i][j][k]->GetMaximum()*1.5);
        hj[i][j][k]->SetLineWidth(2);
        
        hd[i][j][k]->SetDirectory(0);
        hj[i][j][k]->SetDirectory(0);
      }
      f->Close();
    }
  }

  int ijet = 2;
  cout << "Drawing many bins..." << endl;
  int csize = 700;
  TCanvas * cmany1 = new TCanvas("cmany1","",csize*4,csize*2);
  TCanvas * cmany2 = new TCanvas("cmany2","",csize*4,csize*2);
  TCanvas * cmany3 = new TCanvas("cmany3","",csize*4,csize*2);
  TCanvas * cmany4 = new TCanvas("cmany4","",csize*4,csize*2); 
  TH1D * htemp1[ana::nTrijetPtBins];
  TH1D * htemp2[ana::nTrijetPtBins];
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    htemp1[i] = hd[i][ijet][0];
    htemp2[i] = hj[i][ijet][0];
  }
  draw_many(
    cmany1, Form("/home/samson72/sphnx/gammajet/pdfs/note/multijet_%s_nominal.pdf",ana::rnames[ijet]), "p+p #sqrt{s}=200 GeV", "Pythia 8", "Nominal",
    htemp1, htemp2);
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    htemp1[i] = hd[i][ijet][1];
    htemp2[i] = hj[i][ijet][1];
  }
  draw_many(
    cmany2, Form("/home/samson72/sphnx/gammajet/pdfs/note/multijet_%s_JERhigh.pdf",ana::rnames[ijet]), "p+p #sqrt{s}=200 GeV", "Pythia 8", "JER High",
    htemp1, htemp2);
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    htemp1[i] = hd[i][ijet][2];
    htemp2[i] = hj[i][ijet][2];
  }
  draw_many(
    cmany3, Form("/home/samson72/sphnx/gammajet/pdfs/note/multijet_%s_JERlow.pdf",ana::rnames[ijet]), "p+p #sqrt{s}=200 GeV", "Pythia 8", "JER Low",
    htemp1, htemp2);
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    htemp1[i] = hd[i][ijet][3];
    htemp2[i] = hj[i][ijet][3];
  }
  draw_many(
    cmany3, Form("/home/samson72/sphnx/gammajet/pdfs/note/multijet_%s_JERreco.pdf",ana::rnames[ijet]), "p+p #sqrt{s}=200 GeV", "Pythia 8", "JER Reco",
    htemp1, htemp2);
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    htemp1[i] = hd[i][ijet][4];
    htemp2[i] = hj[i][ijet][4];
  }
  draw_many(
    cmany4, Form("/home/samson72/sphnx/gammajet/pdfs/note/multijet_%s_HERWIG.pdf",ana::rnames[ijet]), "p+p #sqrt{s}=200 GeV", "Herwig 7.3", "HERWIG",
    htemp1, htemp2);

  /*
  cout << "Drawing <xj>..." << endl;
  TCanvas * caxj1 = new TCanvas("caxj1","",1000,1000);
  TCanvas * caxj2 = new TCanvas("caxj2","",1000,1000);
  TCanvas * caxj3 = new TCanvas("caxj3","",1000,1000);
  TCanvas * caxj4 = new TCanvas("caxj4","",1000,1000);
  TCanvas * caxj5 = new TCanvas("caxj5","",1000,1000);
  TCanvas * caxj6 = new TCanvas("caxj6","",1000,1000);
  TCanvas * caxj7 = new TCanvas("caxj7","",1000,1000);
  TCanvas * caxj8 = new TCanvas("caxj8","",1000,1000);
  TCanvas * caxj9 = new TCanvas("caxj9","",1000,1000);
  TCanvas * caxj10 = new TCanvas("caxj10","",1000,1000);
  TCanvas * caxj11 = new TCanvas("caxj11","",1000,1000);
  TCanvas * caxj12 = new TCanvas("caxj12","",1000,1000);
  TCanvas * caxj13 = new TCanvas("caxj13","",1000,1000);
  comp_axj(caxj1,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit.pdf","Nominal", 1, // 1 means use fit
    0, 2, 0, 0, 0, "Data", kBlue, 
    1, 2, 0, 0, 0, "MC Photon", kMagenta+1);
  comp_axj(caxj2,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit_narrowBDT.pdf","Narrow BDT score", 1,
    0, 2, 1, 0, 0, "Data", kBlue,
    1, 2, 1, 0, 0, "MC Photon", kMagenta+1);
  comp_axj(caxj2,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit_3jetCut.pdf", "Third Jet Cut", 1,
    0, 2, 0, 1, 0, "Data", kBlue,
    1, 2, 0, 1, 0, "MC Photon", kMagenta+1);
  comp_axj(caxj4,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit_smear.pdf", "JER smeared", 1,
    1, 1, 0, 0, 0, "MC Photon", kMagenta+1,
    1, 2, 0, 0, 0, "MC Photon smeared", kMagenta+4);
  comp_axj(caxj9,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit_narrowIso.pdf", "Narrow Iso cut", 1,
    0, 2, 2, 0, 0, "Data", kBlue,
    1, 2, 2, 0, 0, "MC Photon", kMagenta+4);
  comp_axj(caxj11,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit_Herwig.pdf", "HERWIG-7.3", 1,
    0, 2, 0, 0, 0, "Data", kBlue,
    3, 2, 0, 0, 0, "MC Photon", kMagenta+4);
  comp_axj(caxj5,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean.pdf","Nominal", 0, // 0 means use mean
    0, 2, 0, 0, 0, "Data", kBlue, 
    1, 2, 0, 0, 0, "MC Photon", kMagenta+1);
  comp_axj(caxj6,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_narrowBDT.pdf","Narrow BDT score", 0,
    0, 2, 1, 0, 0, "Data", kBlue,
    1, 2, 1, 0, 0, "MC Photon", kMagenta+1);
  comp_axj(caxj7,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_3jetCut.pdf", "Third Jet Cut", 0,
    0, 2, 0, 1, 0, "Data", kBlue,
    1, 2, 0, 1, 0, "MC Photon", kMagenta+1);
  comp_axj(caxj8,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_smear.pdf", "JER smeared", 0,
    1, 1, 0, 0, 0, "MC Photon", kMagenta+1,
    1, 2, 0, 0, 0, "MC Photon smeared", kMagenta+4);
  comp_axj(caxj10,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_narrowIso.pdf", "Narrow Iso cut", 0,
    0, 2, 2, 0, 0, "Data", kBlue,
    1, 2, 2, 0, 0, "MC Photon", kMagenta+4);
  comp_axj(caxj12,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_Herwig.pdf", "HERWIG-7.3", 0,
    0, 2, 0, 0, 0, "Data", kBlue,
    3, 2, 0, 0, 0, "MC Photon", kMagenta+4);
  comp_axj(caxj13,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_clusterSmear.pdf", "Cluster Smearing", 0,
    0, 2, 0, 0, 0, "Data", kBlue,
    1, 4, 0, 0, 0, "MC Photon", kMagenta+4);
*/  
  
}
