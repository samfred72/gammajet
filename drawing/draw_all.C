#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
//drawer d("herwig");
drawer d("pythia");
drawer dh("herwig");
vector<vector<vector<vector<vector<vector<TF1*>>>>>> fitd   = d.get_empty_TF1();
vector<vector<vector<vector<vector<vector<TF1*>>>>>> fitp   = d.get_empty_TF1(); 
vector<vector<vector<vector<vector<vector<TF1*>>>>>> fitj   = d.get_empty_TF1();
vector<vector<vector<vector<vector<vector<TF1*>>>>>> fith   = d.get_empty_TF1();
vector<vector<vector<vector<vector<vector<float>>>>>> meand = d.get_empty_float();
vector<vector<vector<vector<vector<vector<float>>>>>> meanp = d.get_empty_float(); 
vector<vector<vector<vector<vector<vector<float>>>>>> meanj = d.get_empty_float();
vector<vector<vector<vector<vector<vector<float>>>>>> meanh = d.get_empty_float();
vector<vector<vector<vector<vector<vector<float>>>>>> merrd = d.get_empty_float();
vector<vector<vector<vector<vector<vector<float>>>>>> merrp = d.get_empty_float(); 
vector<vector<vector<vector<vector<vector<float>>>>>> merrj = d.get_empty_float();
vector<vector<vector<vector<vector<vector<float>>>>>> merrh = d.get_empty_float();
const char * histname = "hratio";
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratiod = d.collect_hists(histname,0);
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratiop = d.collect_hists(histname,1);
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratioj = d.get_empty_TH1D();//d.collect_hists(histname,2);
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratioh = dh.collect_hists(histname,1);
float xjm[ana::nJetR];
float xjme[ana::nJetR];
float xjmeb[ana::nJetR];
float xjme3[ana::nJetR];
float xjmei[ana::nJetR];
float xjmeh[ana::nJetR];
float xjf[ana::nJetR];
float xjfe[ana::nJetR];
float xjfeb[ana::nJetR];
float xjfe3[ana::nJetR];
float xjfei[ana::nJetR];
float xjfeh[ana::nJetR];


void draw_wall(TCanvas * c, int iabcd, int k) {

  float drawx = 0.80;
  float drawy = 0.92;
  float fontsize = 60;
  int offset = 0;
  string calibstring = (k ? "JES Calibrated" : "Uncalibrated Jets");
  vector<string> t1 = {"#bf{Analysis region A}:","#bf{Analysis region B}:","#bf{Analysis region C}:","#bf{Analysis region D}:", "Combined ABCD"};

  TPad * p = new TPad("p","",0,0,1,1);
  p->SetLeftMargin(0.25);
  p->SetBottomMargin(0.25);
  p->Divide(ana::nJetR+1,ana::nPtBins,0,0);
  p->Draw();

  for (int i = 0; i < ana::nPtBins; i++) {
    for (int j = 0; j < ana::nJetR; j++) {
      TH1D * h1 = hratiod[i][j][k][0][0][iabcd];
      TH1D * h2 = hratiop[i][j][k][0][0][iabcd];
      TH1D * h3 = hratioj[i][j][k][0][0][iabcd];

      int index = i*ana::nJetR + j + 1 + offset;
      if (index % (ana::nJetR+1) == 0) {index++; offset++;}
      p->cd(index);
      if ((index + 1) % (ana::nJetR+1) == 0) gPad->SetRightMargin(0.01);
      gPad->SetTicks(1,1);
      if (index == (ana::nPtBins*(ana::nJetR+1) - 1)) h1->GetXaxis()->SetTitle("p_{T,max}^{Jet}/p_{T,max}^{cluster}");
      else h1->GetXaxis()->SetTitle("");
      h1->GetXaxis()->SetTitleSize(0.10);
      h1->GetXaxis()->SetLabelSize(0.08);
      if (index == (ana::nJetR+1)*(ana::nPtBins-1)+1) h1->GetXaxis()->SetLabelSize(0.07);
      h1->GetXaxis()->SetNdivisions(405,false);
      h1->GetXaxis()->ChangeLabel(-1,-1,0);
      if (index == 1) h1->GetYaxis()->SetTitle("Arbitrary Units");
      else h1->GetYaxis()->SetTitle("");
      h1->GetYaxis()->SetTitleSize(0.10);
      h1->GetYaxis()->SetLabelSize(0.08);
      if (index == (ana::nJetR+1)*(ana::nPtBins-1)+1) h1->GetYaxis()->SetLabelSize(0.07);
      h1->GetYaxis()->SetLabelOffset(0.04);
      h1->GetYaxis()->SetMaxDigits(3);
      h1->GetYaxis()->SetDecimals(2);
      
      h1->Draw("hist same");
      h2->Draw("hist same");
      //h3->Draw("hist same");
     
      if (iabcd == 0 || iabcd == 4) {
        fitd[i][j][k][0][0][iabcd]->Draw("same");
        fitp[i][j][k][0][0][iabcd]->Draw("same");
        //fitj[i][j][k][0][0][iabcd]->Draw("same");
      }

      float lowcluster = ana::ptBins[i];
      float highcluster = ana::ptBins[i+1];
      float minjet = (k ? ana::jet_calib_pt_cut[j] : ana::jet_pt_cut[j]);
      float drawx = 0.15;
      if ((index - 1) % 5 == 0) drawx = 0.35;
      d.drawMany({
          Form("%.0f GeV < p_{T}^{cluster} < %.0f GeV",lowcluster,highcluster),
          Form("Jet R=%1.1f",ana::JetRs[j])},
          drawx,0.85,42,c->GetWh());
      TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,1);
      line->SetLineStyle(8);
      line->Draw();
    }
  }
  p->cd();
  d.drawAll({
      "run 47289-53864",
      //"MC run28 Jet",
      "MC run28 Photon"},
      {
      Form("|vz| < %.0f cm",ana::vzcut),
      Form("|#eta|_{jet} < %.1f - R",ana::etacut),
      Form("|#eta|_{cluster} < %.1f",ana::etacut),
      Form("#Delta#phi > %.0f#pi/%.0f",ana::oppnum,ana::oppden),
      Form("p_{T}^{jet} > %.0f GeV",(k == 1 ? ana::jet_calib_pt_cut[0] : ana::jet_pt_cut[0])),
      Form("|#DeltaT|_{jet,cluster} < %.0f",ana::tcut),
      calibstring,
      t1[iabcd].c_str()},
      drawx,drawy,fontsize,c->GetWh());
  TLegend * l2 = new TLegend(drawx,drawy-.7,0.99,drawy-.55);
  l2->SetLineWidth(0);
  l2->AddEntry(hratiod[0][0][0][0][0][iabcd],("data"));
  l2->AddEntry(hratiop[0][0][0][0][0][iabcd],("MC Photon reco"));
  //l2->AddEntry(hratioj[0][0][l][ia],("MC Jet reco"));
  l2->Draw();
  return;
}

void draw_many(TCanvas * c, const char * cname, const char * info1, const char * info2, const char * variation,
    int type1, int ir1, int icalib1, int ibdt1, int i3jet1, int iabcd1,
    int type2, int ir2, int icalib2, int ibdt2, int i3jet2, int iabcd2) {
  c->cd();
  vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H1;          
  vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H2;          
  vector<vector<vector<vector<vector<vector<TF1*>>>>>> f1;
  vector<vector<vector<vector<vector<vector<TF1*>>>>>> f2; 
  vector<vector<vector<vector<vector<vector<float>>>>>> mean1;
  vector<vector<vector<vector<vector<vector<float>>>>>> mean2;
  vector<vector<vector<vector<vector<vector<float>>>>>> merr1;
  vector<vector<vector<vector<vector<vector<float>>>>>> merr2;
  if (type1 == 0) {
    H1 = hratiod;
    f1 = fitd;
    mean1 = meand;
    merr1 = merrd;
  }
  else if (type1 == 1) {
    H1 = hratiop;
    f1 = fitp;
    mean1 = meanp;
    merr1 = merrp;
  }
  else if (type1 == 2) {
    H1 = hratioj;
    f1 = fitj;
    mean1 = meanj;
    merr1 = merrj;
  }
  else if (type1 == 3) {
    H1 = hratioh;
    f1 = fith;
    mean1 = meanh;
    merr1 = merrh;
  }
  if (type2 == 0) {
    H2 = hratiod;
    f2 = fitd;
    mean2 = meand;
    merr2 = merrd;
  }
  else if (type2 == 1) {
    H2 = hratiop;
    f2 = fitp;
    mean2 = meanp;
    merr2 = merrp;
  }
  else if (type2 == 2) {
    H2 = hratioj;
    f2 = fitj;
    mean2 = meanj;
    merr2 = merrj;
  }
  else if (type2 == 3) {
    H2 = hratioh;
    f2 = fith;
    mean2 = meanh;
    merr2 = merrh;
  }


  float drawx = 0.77;
  float drawy = 0.92;
  float fontsize = 50;
  int offset = 0;
  int ir = 1;
  string calibstring = (icalib1 ? "JES Calibrated" : "Uncalibrated Jets");
  vector<string> t1 = {"#bf{Analysis region A}","#bf{Analysis region B}","#bf{Analysis region C}","#bf{Analysis region D}", "Combined ABCD"};

  TPad * p = new TPad("p","",0,0,1,1);
  p->SetLeftMargin(0.25);
  p->SetBottomMargin(0.25);
  p->Divide(4,3,0,0);
  p->Draw();

  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * h1 = H1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1];
    TH1D * h2 = H2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2];

    int index = ipt + 1 + offset;
    if (index % (3+1) == 0) {index++; offset++;}
    p->cd(index);
    if ((index + 1) % (3+1) == 0) gPad->SetRightMargin(0.01);
    gPad->SetTicks(1,1);
    if (index == (3*(3+1) - 1)) h1->GetXaxis()->SetTitle("p_{T,max}^{Jet}/p_{T,max}^{cluster}");
    else h1->GetXaxis()->SetTitle("");
    h1->GetXaxis()->SetTitleSize(0.10);
    h1->GetXaxis()->SetLabelSize(0.08);
    if (index == (3+1)*(3-1)+1) h1->GetXaxis()->SetLabelSize(0.07);
    h1->GetXaxis()->SetNdivisions(405,false);
    h1->GetXaxis()->ChangeLabel(-1,-1,0);
    if (index == 1) h1->GetYaxis()->SetTitle("Arbitrary Units");
    else h1->GetYaxis()->SetTitle("");
    h1->GetYaxis()->SetTitleSize(0.10);
    h1->GetYaxis()->SetLabelSize(0.08);
    if (index == (3+1)*(ana::nPtBins-1)+1) h1->GetYaxis()->SetLabelSize(0.07);
    h1->GetYaxis()->SetLabelOffset(0.04);
    h1->GetYaxis()->SetMaxDigits(3);
    h1->GetYaxis()->SetDecimals(2);

    h1->Draw("hist same");
    h2->Draw("hist same");

    if (iabcd1 == 0 || iabcd1 == 4) {
      f1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1]->Draw("same");
    }
    if (iabcd2 == 0 || iabcd2 == 4) {
      f2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2]->Draw("same");
    }

    float lowcluster = ana::ptBins[ipt];
    float highcluster = ana::ptBins[ipt+1];
    float minjet = (icalib1 ? ana::jet_calib_pt_cut[ir1] : ana::jet_pt_cut[ir1]);
    float drawx = 0.15;
    if ((index - 1) % 4 == 0) drawx = 0.35;
    d.drawMany({
        Form("#bf{%.0f GeV < p_{T}^{cluster} < %.0f GeV}",lowcluster,highcluster),
        //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info1,  f1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1),f1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1)),
        //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info2,  f2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1),f2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1))
        Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info1,  mean1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1],merr1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1]),
        Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info2,  mean2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2],merr2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2])
        },drawx,0.88,42,c->GetWh()/3.0);
    TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,1);
    line->SetLineStyle(8);
    line->Draw();
    TLine * mline1 = new TLine(mean1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1],0,mean1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1],1);
    mline1->SetLineColor(h1->GetLineColor());
    mline1->SetLineStyle(9);
    TLine * mline2 = new TLine(mean2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2],0,mean2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2],1);
    mline2->SetLineColor(h2->GetLineColor());
    mline2->SetLineStyle(9);
    mline1->Draw();
    mline2->Draw();
  }
  p->cd();
  d.drawAll({
      info1, 
      info2},
      {
      "Analysis cuts",
      Form("Jet R=%1.1f",ana::JetRs[ir1]),
      variation,
      calibstring,
      t1[iabcd1].c_str()},
      drawx,drawy,fontsize,c->GetWh());
  TLegend * l2 = new TLegend(drawx,drawy-.3,0.99,drawy-.2);
  l2->SetLineWidth(0);
  l2->AddEntry(H1[0][ir1][icalib1][ibdt1][i3jet1][iabcd1],info1);
  l2->AddEntry(H2[0][ir2][icalib2][ibdt2][i3jet2][iabcd2],info2);
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
  int colors[6] = {kSpring + 2, kBlue, kGreen + 3, kMagenta+1, kTeal, kSpring+2};
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
  line2->SetLineColor(kMagenta+1);
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


void comp_axj(TCanvas * c, const char * cname, const char * info, bool usefit,
              int type1, int icalib1, int ibdt1, int i3jet1, int iabcd1, string info1, int color1, 
              int type2, int icalib2, int ibdt2, int i3jet2, int iabcd2, string info2, int color2) {
  c->cd();
  vector<vector<vector<vector<vector<vector<TF1 *>>>>>> f1;
  vector<vector<vector<vector<vector<vector<TF1 *>>>>>> f2; 
  vector<vector<vector<vector<vector<vector<float>>>>>> mean1;
  vector<vector<vector<vector<vector<vector<float>>>>>> mean2;
  vector<vector<vector<vector<vector<vector<float>>>>>> merr1;
  vector<vector<vector<vector<vector<vector<float>>>>>> merr2;
  if (type1 == 0) {
    f1 = fitd;
    mean1 = meand;
    merr1 = merrd;
  }
  else if (type1 == 1) {
    f1 = fitp;
    mean1 = meanp;
    merr1 = merrp;
  }
  else if (type1 == 2) {
    f1 = fitj;
    mean1 = meanj;
    merr1 = merrj;
  }
  else if (type1 == 3) {
    f1 = fith;
    mean1 = meanh;
    merr1 = merrh;
  }
  if (type2 == 0) {
    f2 = fitd;
    mean2 = meand;
    merr2 = merrd;
  }
  else if (type2 == 1) {
    f2 = fitp;
    mean2 = meanp;
    merr2 = merrp;
  }
  else if (type2 == 2) {
    f2 = fitj;
    mean2 = meanj;
    merr2 = merrj;
  }
  else if (type2 == 3) {
    f2 = fith;
    mean2 = meanh;
    merr2 = merrh;
  }

  TFile * ftemp = 0;
  if (usefit == 0 && strcmp(info,"Nominal") == 0) {
    ftemp = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu.root"),"RECREATE");
  }
  c->SaveAs(Form("%s[",cname));
  for (int ir = 0; ir < ana::nJetR; ir++) {
    TPad * p1 = new TPad("p1","",0,.6,1,1);
    TPad * p2 = new TPad("p2","",0,0,1,0.6);
    p1->Draw();
    p2->Draw();
    p1->SetBottomMargin(0);
    p2->SetTopMargin(0);
    p1->cd();
    gPad->SetTicks(1,1);
    TLegend * lxj = new TLegend(.15,.65,.4,.85);
    lxj->SetLineWidth(0);
    TH1D * h1 = new TH1D(Form("h_%s_%s_%s_%i",info1.c_str(),info,f1[0][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetName(),usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * h2 = new TH1D(Form("h_%s_%s_%s_%i",info2.c_str(),info,f1[0][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetName(),usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      if (usefit) {
        h1->SetBinContent(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1));
        h2->SetBinContent(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1));
        h1->SetBinError(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1));
        h2->SetBinError(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1));
      }
      else {
        h1->SetBinContent(ipt+1,mean1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]);
        h2->SetBinContent(ipt+1,mean2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]);
        h1->SetBinError(ipt+1,merr1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]);
        h2->SetBinError(ipt+1,merr2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]);
      }

    }
    
    h1->GetYaxis()->SetRangeUser(0,2);
    h1->GetYaxis()->SetTitleSize(0.06);
    h1->GetYaxis()->SetLabelSize(0.06);
    h1->GetYaxis()->SetTitleOffset(0.5);
    h1->SetLineColor(color1);
    h1->SetMarkerColor(color1);
    h1->SetMarkerStyle(20);
    h2->SetLineColor(color2);
    h2->SetMarkerColor(color2);
    h2->SetMarkerStyle(20);
    h1->Draw("p");
    h2->Draw("p same");
    lxj->AddEntry(h1,info1.c_str());
    lxj->AddEntry(h2,info2.c_str());
    lxj->Draw("same");
    const char * fittext = (usefit ? "Fit method" : "Mean method");
    d.drawAll(
        {
        info1,
        info2
        },
        {
        "Analysis cuts",
        fittext,
        Form("Jet R=%0.1f",ana::JetRs[ir]), 
        },
        .5, .8, 20, c->GetWh()*0.6);

    TBox *texclude1 = new TBox(10,0,13.0,2);
    texclude1->SetFillColorAlpha(kGray,0.3);
    //texclude1->Draw("same");
    p2->cd();
    gPad->SetBottomMargin(0.2);
    gPad->SetTicks(1,1);
    TH1D * dummy = (TH1D*)h1->Clone();
    dummy->SetName(Form("d%i",ir));
    dummy->Reset("ICES");

    dummy->GetYaxis()->SetTitle("Ratio data/MC");
    dummy->GetYaxis()->SetTitleSize(0.04);
    dummy->GetYaxis()->SetLabelSize(0.04);
    dummy->GetYaxis()->SetTitleOffset(1);
    dummy->GetYaxis()->SetRangeUser(.8,1.2);

    dummy->GetXaxis()->SetTitle("Cluster p_{T}");
    dummy->GetXaxis()->SetLabelSize(0.06);
    dummy->GetXaxis()->SetTitleSize(0.06);
    dummy->GetXaxis()->SetTitleOffset(1);

    dummy->Draw();
    TH1D * hdivide = (TH1D*)h1->Clone();
    hdivide->SetName(Form("hdivide%i",ir));
    hdivide->Divide(h1,h2);
    hdivide->SetMarkerColor(kBlack);
    hdivide->SetLineColor(kBlack);
    hdivide->SetMarkerStyle(20);
    hdivide->SetMarkerSize(1);
    hdivide->Draw("same");
    d.drawLine(ana::ptBins[0],1,ana::ptBins[ana::nPtBins],1); 
    if (usefit == 0 && strcmp(info,"Nominal") == 0) {
      hdivide->Write();
    }
    
    TBox *texclude2 = new TBox(10,0.8,13.0,1.2);
    texclude2->SetFillColorAlpha(kGray,0.3);
    //texclude2->Draw("same");

    TF1 * fline = new TF1(Form("func_%s",cname),"pol0",10,35);
    hdivide->Fit(fline,"RQIM0");
    fline->SetLineStyle(9);
    fline->SetLineWidth(3);
    fline->Draw("same");
    
    float flow = fline->GetParameter(0) - fline->GetParError(0);
    float fhigh = fline->GetParameter(0) + fline->GetParError(0);
    TBox *terr = new TBox(10,flow,35,fhigh);
    terr->SetFillColorAlpha(kRed,0.3);
    terr->Draw("same");

    if (strcmp(info, "Nominal")) d.drawText(Form("Variation: %s",info),.25,.35,kRed,40);
    d.drawText(Form("#bf{#it{in situ} JES = %2.4f #pm %2.4f}",fline->GetParameter(0), fline->GetParError(0)),.25,.28,kRed,40);

    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
  if (ftemp) ftemp->Close();
}

void comp_comp_axj(TCanvas * c, const char * cname, const char * variation, bool usefit,
              int type1, int icalib1, int ibdt1, int i3jet1, int iabcd1, string info1, int color1, 
              int type2, int icalib2, int ibdt2, int i3jet2, int iabcd2, string info2, int color2) {
  c->cd();
  vector<vector<vector<vector<vector<vector<TF1 *>>>>>> f1;
  vector<vector<vector<vector<vector<vector<TF1 *>>>>>> f2; 
  vector<vector<vector<vector<vector<vector<float>>>>>> mean1;
  vector<vector<vector<vector<vector<vector<float>>>>>> mean2;
  vector<vector<vector<vector<vector<vector<float>>>>>> merr1;
  vector<vector<vector<vector<vector<vector<float>>>>>> merr2;
  if (type1 == 0) {
    f1 = fitd;
    mean1 = meand;
    merr1 = merrd;
  }
  else if (type1 == 1) {
    f1 = fitp;
    mean1 = meanp;
    merr1 = merrp;
  }
  else if (type1 == 2) {
    f1 = fitj;
    mean1 = meanj;
    merr1 = merrj;
  }
  else if (type1 == 3) {
    f1 = fith;
    mean1 = meanh;
    merr1 = merrh;
  }
  if (type2 == 0) {
    f2 = fitd;
    mean2 = meand;
    merr2 = merrd;
  }
  else if (type2 == 1) {
    f2 = fitp;
    mean2 = meanp;
    merr2 = merrp;
  }
  else if (type2 == 2) {
    f2 = fitj;
    mean2 = meanj;
    merr2 = merrj;
  }
  else if (type2 == 3) {
    f2 = fith;
    mean2 = meanh;
    merr2 = merrh;
  }
  c->SaveAs(Form("%s[",cname));
  for (int ir = 0; ir < ana::nJetR; ir++) {
    gPad->SetTicks(1,1);
    TH1D * h1 = new TH1D(Form("hc_%s_%s_%s_%i",info1.c_str(),variation,f1[0][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetName(),usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * h2 = new TH1D(Form("hc_%s_%s_%s_%i",info2.c_str(),variation,f1[0][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetName(),usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      if (usefit) {
        h1->SetBinContent(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1));
        h2->SetBinContent(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1));
        h1->SetBinError(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1));
        h2->SetBinError(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1));
      }
      else {
        h1->SetBinContent(ipt+1,mean1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]);
        h2->SetBinContent(ipt+1,mean2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]);
        h1->SetBinError(ipt+1,merr1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]);
        h2->SetBinError(ipt+1,merr2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]);
      }
    }
    // c stands for "compare"
    TH1D * h1c = new TH1D(Form("hcc_%s_%s_%s_%i",info1.c_str(),variation,f1[0][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetName(),usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * h2c = new TH1D(Form("hcc_%s_%s_%s_%i",info2.c_str(),variation,f1[0][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetName(),usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      if (usefit) {
        h1c->SetBinContent(ipt+1,fitd[ipt][ir][2][0][0][0]->GetParameter(1));
        h2c->SetBinContent(ipt+1,fitp[ipt][ir][2][0][0][0]->GetParameter(1));
        h1c->SetBinError(ipt+1,fitd[ipt][ir][2][0][0][0]->GetParError(1));
        h2c->SetBinError(ipt+1,fitp[ipt][ir][2][0][0][0]->GetParError(1));
      }
      else {
        h1c->SetBinContent(ipt+1,meand[ipt][ir][2][0][0][0]);
        h2c->SetBinContent(ipt+1,meanp[ipt][ir][2][0][0][0]);
        h1c->SetBinError(ipt+1,merrd[ipt][ir][2][0][0][0]);
        h2c->SetBinError(ipt+1,merrp[ipt][ir][2][0][0][0]);
      }
    }
    
    TH1D * dummy = (TH1D*)h1->Clone();
    dummy->SetName(Form("d%i",ir));
    dummy->Reset("ICES");

    gPad->SetBottomMargin(.15);
    dummy->GetYaxis()->SetTitle(Form("Double Ratio Data/MC to %s Data/MC",variation));
    dummy->GetYaxis()->SetTitleSize(0.04);
    dummy->GetYaxis()->SetLabelSize(0.04);
    dummy->GetYaxis()->SetTitleOffset(1);
    dummy->GetYaxis()->SetRangeUser(.8,1.2);

    dummy->GetXaxis()->SetTitle("Cluster p_{T}");
    dummy->GetXaxis()->SetLabelSize(0.04);
    dummy->GetXaxis()->SetTitleSize(0.04);
    dummy->GetXaxis()->SetTitleOffset(1);

    dummy->Draw();
    TH1D * hdivide1 = (TH1D*)h1->Clone();
    hdivide1->SetName(Form("hdivide1%i",ir));
    TH1D * hdivide2 = (TH1D*)h1->Clone();
    hdivide2->SetName(Form("hdivide2%i",ir));
    TH1D * hdivide = (TH1D*)h1->Clone();
    hdivide->SetName(Form("hdivide%i",ir));

    hdivide1->Divide(h1,h2);
    hdivide2->Divide(h1c,h2c);
    hdivide->Divide(hdivide1,hdivide2);
    hdivide->SetMarkerColor(kBlack);
    hdivide->SetLineColor(kBlack);
    hdivide->SetMarkerStyle(21);
    hdivide->SetMarkerSize(1);
    hdivide->Draw("same");
    d.drawLine(ana::ptBins[0],1,ana::ptBins[ana::nPtBins],1); 
    
    d.drawAll(
        {
        },
        {
        "Analysis cuts",
        Form("Jet R=%0.1f",ana::JetRs[ir]), 
        },
        .25, .8, 20, c->GetWh());
    
    TBox *texclude = new TBox(10,0.8,13.0,1.2);
    texclude->SetFillColorAlpha(kGray,0.3);
    //texclude->Draw("same");

    TF1 * fline = new TF1(Form("func_%s",variation),"pol0",10,35);
    hdivide->Fit(fline,"RQIM0");
    fline->SetLineColor(kSpring+2);
    fline->SetLineStyle(9);
    fline->SetLineWidth(3);
    fline->Draw("same");
    
    float flow = fline->GetParameter(0) - fline->GetParError(0);
    float fhigh = fline->GetParameter(0) + fline->GetParError(0);
    TBox *terr = new TBox(10,flow,35,fhigh);
    terr->SetFillColorAlpha(kBlue,0.3);
    terr->Draw("same");

    if (strcmp(variation, "Nominal")) d.drawText(Form("Variation: %s",variation),.20,.25,kBlue,40);
    d.drawText(Form("#bf{Relative #it{in situ} JES = %2.4f #pm %2.4f}",fline->GetParameter(0), fline->GetParError(0)),.20,.18,kBlue,40);

    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
}

void fillxj(bool usefit) {
  for (int ir = 0; ir < ana::nJetR; ir++) {
    TH1D * hd = new TH1D(Form("hxjd_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hp = new TH1D(Form("hxjp_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hh = new TH1D(Form("hxjh_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hdb = new TH1D(Form("hxjbd_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hpb = new TH1D(Form("hxjbp_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hhb = new TH1D(Form("hxjbh_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hd3 = new TH1D(Form("hxj3d_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hp3 = new TH1D(Form("hxj3p_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hh3 = new TH1D(Form("hxj3h_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hdi = new TH1D(Form("hxjid_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hpi = new TH1D(Form("hxjip_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * hhi = new TH1D(Form("hxjih_%i_%i_fit",ir, usefit),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      if (usefit) {
        hd->SetBinContent(ipt+1,fitd[ipt][ir][2][0][0][0]->GetParameter(1));
        hp->SetBinContent(ipt+1,fitp[ipt][ir][2][0][0][0]->GetParameter(1));
        hh->SetBinContent(ipt+1,fith[ipt][ir][2][0][0][0]->GetParameter(1));
        hd->SetBinError(  ipt+1,fitd[ipt][ir][2][0][0][0]->GetParError(1));
        hp->SetBinError(  ipt+1,fitp[ipt][ir][2][0][0][0]->GetParError(1));
        hh->SetBinError(  ipt+1,fith[ipt][ir][2][0][0][0]->GetParError(1));

        hdb->SetBinContent(ipt+1,fitd[ipt][ir][2][1][0][0]->GetParameter(1));
        hpb->SetBinContent(ipt+1,fitp[ipt][ir][2][1][0][0]->GetParameter(1));
        hhb->SetBinContent(ipt+1,fith[ipt][ir][2][1][0][0]->GetParameter(1));
        hdb->SetBinError(  ipt+1,fitd[ipt][ir][2][1][0][0]->GetParError(1));
        hpb->SetBinError(  ipt+1,fitp[ipt][ir][2][1][0][0]->GetParError(1));
        hhb->SetBinError(  ipt+1,fith[ipt][ir][2][1][0][0]->GetParError(1));

        hd3->SetBinContent(ipt+1,fitd[ipt][ir][2][0][1][0]->GetParameter(1));
        hp3->SetBinContent(ipt+1,fitp[ipt][ir][2][0][1][0]->GetParameter(1));
        hh3->SetBinContent(ipt+1,fith[ipt][ir][2][0][1][0]->GetParameter(1));
        hd3->SetBinError(  ipt+1,fitd[ipt][ir][2][0][1][0]->GetParError(1));
        hp3->SetBinError(  ipt+1,fitp[ipt][ir][2][0][1][0]->GetParError(1));
        hh3->SetBinError(  ipt+1,fith[ipt][ir][2][0][1][0]->GetParError(1));
        
        hdi->SetBinContent(ipt+1,fitd[ipt][ir][2][2][0][0]->GetParameter(1));
        hpi->SetBinContent(ipt+1,fitp[ipt][ir][2][2][0][0]->GetParameter(1));
        hhi->SetBinContent(ipt+1,fith[ipt][ir][2][2][0][0]->GetParameter(1));
        hdi->SetBinError(  ipt+1,fitd[ipt][ir][2][2][0][0]->GetParError(1));
        hpi->SetBinError(  ipt+1,fitp[ipt][ir][2][2][0][0]->GetParError(1));
        hhi->SetBinError(  ipt+1,fith[ipt][ir][2][2][0][0]->GetParError(1));
      }
      else {
        hd->SetBinContent(ipt+1,meand[ipt][ir][2][0][0][0] );
        hp->SetBinContent(ipt+1,meanp[ipt][ir][2][0][0][0] );
        hh->SetBinContent(ipt+1,meanh[ipt][ir][2][0][0][0] );
        hd->SetBinError(  ipt+1,merrd[ipt][ir][2][0][0][0] );
        hp->SetBinError(  ipt+1,merrp[ipt][ir][2][0][0][0] );
        hh->SetBinError(  ipt+1,merrh[ipt][ir][2][0][0][0] );

        hdb->SetBinContent(ipt+1,meand[ipt][ir][2][1][0][0]);
        hpb->SetBinContent(ipt+1,meanp[ipt][ir][2][1][0][0]);
        hhb->SetBinContent(ipt+1,meanh[ipt][ir][2][1][0][0]);
        hdb->SetBinError(  ipt+1,merrd[ipt][ir][2][1][0][0]);
        hpb->SetBinError(  ipt+1,merrp[ipt][ir][2][1][0][0]);
        hhb->SetBinError(  ipt+1,merrh[ipt][ir][2][1][0][0]);

        hd3->SetBinContent(ipt+1,meand[ipt][ir][2][0][1][0]);
        hp3->SetBinContent(ipt+1,meanp[ipt][ir][2][0][1][0]);
        hh3->SetBinContent(ipt+1,meanh[ipt][ir][2][0][1][0]);
        hd3->SetBinError(  ipt+1,merrd[ipt][ir][2][0][1][0]);
        hp3->SetBinError(  ipt+1,merrp[ipt][ir][2][0][1][0]);
        hh3->SetBinError(  ipt+1,merrh[ipt][ir][2][0][1][0]);
        
        hdi->SetBinContent(ipt+1,meand[ipt][ir][2][2][0][0]);
        hpi->SetBinContent(ipt+1,meanp[ipt][ir][2][2][0][0]);
        hhi->SetBinContent(ipt+1,meanh[ipt][ir][2][2][0][0]);
        hdi->SetBinError(  ipt+1,merrd[ipt][ir][2][2][0][0]);
        hpi->SetBinError(  ipt+1,merrp[ipt][ir][2][2][0][0]);
        hhi->SetBinError(  ipt+1,merrh[ipt][ir][2][2][0][0]);
      }
    }
    TH1D * hdivide = (TH1D*)hd->Clone();
    hdivide->SetName(Form("hdivide%i%i",ir,usefit));
    hdivide->Divide(hd,hp);
    
    TH1D * hdivideb = (TH1D*)hd->Clone();
    hdivideb->SetName(Form("hdivideb%i%i",ir,usefit));
    hdivideb->Divide(hdb,hpb);
    
    TH1D * hdivide3 = (TH1D*)hd->Clone();
    hdivide3->SetName(Form("hdivide3%i%i",ir,usefit));
    hdivide3->Divide(hd3,hp3);
    
    TH1D * hdividei = (TH1D*)hd->Clone();
    hdividei->SetName(Form("hdividei%i%i",ir,usefit));
    hdividei->Divide(hdi,hpi);
    
    TH1D * hdivideh = (TH1D*)hd->Clone();
    hdivideh->SetName(Form("hdivideh%i%i",ir,usefit));
    hdivideh->Divide(hd,hh);

    TH1D * hcompb = (TH1D*)hd->Clone();
    hcompb->SetName(Form("hcompb%i%i",ir,usefit));
    hcompb->Divide(hdivide,hdivideb);
    
    TH1D * hcomp3 = (TH1D*)hd->Clone();
    hcomp3->SetName(Form("hcomp3%i%i",ir,usefit));
    hcomp3->Divide(hdivide,hdivide3);
    
    TH1D * hcompi = (TH1D*)hd->Clone();
    hcompi->SetName(Form("hcompi%i%i",ir,usefit));
    hcompi->Divide(hdivide,hdividei);
    
    TH1D * hcomph = (TH1D*)hd->Clone();
    hcomph->SetName(Form("hcomph%i%i",ir,usefit));
    hcomph->Divide(hdivide,hdivideh);
    
    TF1 * fline = new TF1(Form("func_%i_%i",ir,usefit),"pol0",10,35);
    hdivide->Fit(fline,"RQIM0");   
    
    TF1 * flineb = new TF1(Form("funcb_%i_%i",ir,usefit),"pol0",10,35);
    hcompb->Fit(flineb,"RQIM0");   
    
    TF1 * fline3 = new TF1(Form("func3_%i_%i",ir,usefit),"pol0",10,35);
    hcomp3->Fit(fline3,"RQIM0"); 
    
    TF1 * flinei = new TF1(Form("funci_%i_%i",ir,usefit),"pol0",10,35);
    hcompi->Fit(flinei,"RQIM0"); 
    
    TF1 * flineh = new TF1(Form("funch_%i_%i",ir,usefit),"pol0",10,35);
    hcomph->Fit(flineh,"RQIM0"); 

    if (usefit) {
      xjf  [ir] = fline->GetParameter(0);
      xjfe [ir] = fline->GetParError(0);
      xjfeb[ir] = abs(1-flineb->GetParameter(0));
      xjfe3[ir] = abs(1-fline3->GetParameter(0));
      xjfei[ir] = abs(1-flinei->GetParameter(0));
      xjfeh[ir] = abs(1-flineh->GetParameter(0));
    }
    else {
      xjm  [ir] = fline->GetParameter(0);
      xjme [ir] = fline->GetParError(0);
      xjmeb[ir] = abs(1-flineb->GetParameter(0));
      xjme3[ir] = abs(1-fline3->GetParameter(0));
      xjmei[ir] = abs(1-flinei->GetParameter(0));
      xjmeh[ir] = abs(1-flineh->GetParameter(0));
    }
  }
  return;
}

void draw_all() {
  gStyle->SetOptStat(0);
  if (!gROOT->IsBatch()) {
    cout << "Run with -b flag or else!!" << endl;
    gROOT->SetBatch(kTRUE);
    //return;
  }
  
  cout << "formatting hists/funcs..." << endl;
  for (int i = 0; i < ana::nPtBins; i++ ) {
    for (int j = 0; j < ana::nJetR; j++ ) {
      for (int k = 0; k < ana::nCalibBins; k++ ) {
        for (int l = 0; l < ana::nIsoBdtBins; l++ ) {
          for (int m = 0; m < ana::n3jetBins; m++ ) {
            for (int n = 0; n < 5; n++ ) {
              //cout << "Working on " << i << " " << j << " " << k << " " << l << " " << m << " " << n << endl;
              float lowcluster = ana::ptBins[i];
              float highcluster = ana::ptBins[i+1];
              float minjet = (k ? ana::jet_calib_pt_cut[j] : ana::jet_pt_cut[j]);
              int lowbin = (int)(minjet/lowcluster/0.08 + 1);
              if (n == 0 || n == 4) {
                fitd[i][j][k][l][m][n] = d.fit(hratiod[i][j][k][l][m][n], lowbin*0.08,2, "RLQI0");
                fitd[i][j][k][l][m][n]->SetParameter(0,fitd[i][j][k][l][m][n]->GetParameter(0)/hratiod[i][j][k][l][m][n]->GetEntries());
                d.format(fitd[i][j][k][l][m][n],0);
              }

              d.format(hratiod[i][j][k][l][m][n],0);
              d.format(hratiop[i][j][k][l][m][n],1);
              dh.format(hratioh[i][j][k][l][m][n],2);
              //d.format(hratioj[i][j][k][l][m][n],2);
              
              if (n == 0 || n == 4) {
                fitp[i][j][k][l][m][n] = d.fit(hratiop[i][j][k][l][m][n], lowbin*0.08,2,"RMQI0");
                fith[i][j][k][l][m][n] = dh.fit(hratioh[i][j][k][l][m][n], lowbin*0.08,2,"RMQI0");
                //fitj[i][j][k][l][m][n] = d.fit(hratioj[i][j][k][l][m][n], (int)(minjet/lowcluster/0.04 + 1)*0.04,2);
                d.format(fitp[i][j][k][l][m][n],1);
                dh.format(fith[i][j][k][l][m][n],1);
              
                hratiod[i][j][k][l][m][n]->GetXaxis()->SetRange(lowbin+1,25);
                hratiop[i][j][k][l][m][n]->GetXaxis()->SetRange(lowbin+1,25);
                hratioh[i][j][k][l][m][n]->GetXaxis()->SetRange(lowbin+1,25);
                meand[i][j][k][l][m][n] = hratiod[i][j][k][l][m][n]->GetMean();
                merrd[i][j][k][l][m][n] = hratiod[i][j][k][l][m][n]->GetMeanError();
                meanp[i][j][k][l][m][n] = hratiop[i][j][k][l][m][n]->GetMean();
                merrp[i][j][k][l][m][n] = hratiop[i][j][k][l][m][n]->GetMeanError();
                meanh[i][j][k][l][m][n] = hratioh[i][j][k][l][m][n]->GetMean();
                merrh[i][j][k][l][m][n] = hratioh[i][j][k][l][m][n]->GetMeanError();
                hratiod[i][j][k][l][m][n]->GetXaxis()->SetRange(1,25);
                hratiop[i][j][k][l][m][n]->GetXaxis()->SetRange(1,25);
                hratioh[i][j][k][l][m][n]->GetXaxis()->SetRange(1,25);
              }
            }
          }
        }
      }
    }
  }

  cout << "Drawing calib plots..." << endl;
  int csize =  700;
  TCanvas * c = new TCanvas("c","",csize*ana::nJetR+200,csize*ana::nPtBins);
  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/abcd_calib.pdf["));
  // 5 for A,B,C,D, and the combined one
  for (int iabcd = 0; iabcd < 5; iabcd++) {
    c->Clear();
    draw_wall(c,iabcd,2);
    c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/abcd_calib.pdf"));
  }
  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/abcd_calib.pdf]"));
  delete c;

  cout << "Drawing single bin..." << endl;
  int i = 5;
  int j = 1;
  TCanvas * cone_calib = new TCanvas("cone_calib","",csize,csize); 
  cone_calib->SaveAs("/home/samson72/sphnx/gammajet/pdfs/oneabcd_calib.pdf[");
  for (int iabcd = 0; iabcd < 5; iabcd++) {
    cone_calib->Clear();
    draw_one(cone_calib,cone_calib->GetName(),"Data","MC Photon", hratiod,hratiop,fitd,fitp,meand,meanp,i,j,2,0,0,iabcd,i,j,2,0,0,iabcd);
    cone_calib->SaveAs("/home/samson72/sphnx/gammajet/pdfs/oneabcd_calib.pdf");
  }
  cone_calib->SaveAs("/home/samson72/sphnx/gammajet/pdfs/oneabcd_calib.pdf]");

  cout << "Drawing many bins..." << endl;
  TCanvas * cmany1 = new TCanvas("cmany1","",csize*4,csize*3);
  TCanvas * cmany2 = new TCanvas("cmany2","",csize*4,csize*3);
  TCanvas * cmany3 = new TCanvas("cmany3","",csize*4,csize*3);
  TCanvas * cmany4 = new TCanvas("cmany4","",csize*4,csize*3);
  TCanvas * cmany5 = new TCanvas("cmany5","",csize*4,csize*3);
  TCanvas * cmany6 = new TCanvas("cmany6","",csize*4,csize*3);
  TCanvas * cmany7 = new TCanvas("cmany7","",csize*4,csize*3);
  TCanvas * cmany8 = new TCanvas("cmany8","",csize*4,csize*3);
  draw_many(
    cmany1, "/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA.pdf", "Data", "MC Photon", "Nominal",
    0, 1, 2, 0, 0, 0,
    1, 1, 2, 0, 0, 0);
  draw_many(
    cmany2, "/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA_narrowBDT.pdf", "Data", "MC Photon", "Narrow BDT",
    0, 1, 2, 1, 0, 0,
    1, 1, 2, 1, 0, 0);
  draw_many(
    cmany3, "/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA_3jetCut.pdf", "Data", "MC Photon", "Third Jet cut",
    0, 1, 2, 0, 1, 0,
    1, 1, 2, 0, 1, 0);
  draw_many(
    cmany4, "/home/samson72/sphnx/gammajet/pdfs/allptbins_MC_R04_regionA_smear.pdf", "MC Photon", "MC smeared", "JER comparison",
    1, 1, 1, 0, 0, 0,
    1, 1, 2, 0, 0, 0);
  draw_many(
    cmany5, "/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA_narrowIso.pdf", "Data", "MC Photon", "Narrow Isolation Energy",
    0, 1, 2, 2, 0, 0,
    1, 1, 2, 2, 0, 0);
  draw_many(
    cmany7, "/home/samson72/sphnx/gammajet/pdfs/allptbins_R08_regionA_Herwig_Pythia.pdf", "Data", "MC Photon", "HERWIG-7.3",
    1, 1, 2, 0, 0, 0,
    3, 1, 2, 0, 0, 0);
  draw_many(
    cmany8, "/home/samson72/sphnx/gammajet/pdfs/allptbins_R08_regionA_Reweight_Pythia.pdf", "MC Photon", "MC Photon Reweighted", "Reweighted",
    1, 1, 2, 0, 0, 0,
    1, 1, 3, 0, 0, 0);

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
  
  

  cout << "Drawing <xj> comparisons..." << endl;
  TCanvas * caxj1comp = new TCanvas("caxj1comp","",1000,600);
  TCanvas * caxj2comp = new TCanvas("caxj2comp","",1000,600);
  TCanvas * caxj3comp = new TCanvas("caxj3comp","",1000,600);
  TCanvas * caxj4comp = new TCanvas("caxj4comp","",1000,600);
  TCanvas * caxj5comp = new TCanvas("caxj5comp","",1000,600);
  TCanvas * caxj6comp = new TCanvas("caxj6comp","",1000,600);
  TCanvas * caxj7comp = new TCanvas("caxj7comp","",1000,600);
  TCanvas * caxj8comp = new TCanvas("caxj8comp","",1000,600);
  TCanvas * caxj9comp = new TCanvas("caxj9comp","",1000,600);
  TCanvas * caxj10comp = new TCanvas("caxj10comp","",1000,600);
  TCanvas * caxj11comp = new TCanvas("caxj11comp","",1000,600);
  comp_comp_axj(caxj1comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_fit_comp.pdf","Nominal", 1,
    0, 2, 0, 0, 0, "Data", kBlue,
    1, 2, 0, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj2comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_narrowBDT_fit_comp.pdf", "Narrow BDT score", 1,
    0, 2, 1, 0, 0, "Data", kBlue,
    1, 2, 1, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj3comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_3jetCut_fit_comp.pdf", "Third Jet cut", 1,
    0, 2, 0, 1, 0, "Data", kBlue,
    1, 2, 0, 1, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj7comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_narrowIso_fit_comp.pdf", "Narrow Iso cut", 1,
    0, 2, 2, 0, 0, "Data", kBlue,
    1, 2, 2, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj9comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_Herwig_fit_comp.pdf", "HERWIG-7.3", 1,
    0, 2, 0, 0, 0, "Data", kBlue,
    3, 2, 0, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj4comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_mean_comp.pdf","Nominal", 0,
    0, 2, 0, 0, 0, "Data", kBlue,
    1, 2, 0, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj5comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_narrowBDT_mean_comp.pdf", "Narrow BDT score", 0,
    0, 2, 1, 0, 0, "Data", kBlue,
    1, 2, 1, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj6comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_3jetCut_mean_comp.pdf", "Third Jet cut", 0,
    0, 2, 0, 1, 0, "Data", kBlue,
    1, 2, 0, 1, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj8comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_narrowIso_mean_comp.pdf", "Narrow Iso cut", 0,
    0, 2, 2, 0, 0, "Data", kBlue,
    1, 2, 2, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj10comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_Herwig_mean_comp.pdf", "HERWIG-7.3", 0,
    0, 2, 0, 0, 0, "Data", kBlue,
    3, 2, 0, 0, 0, "MC Photon", kMagenta+1);
  comp_comp_axj(caxj11comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_clusterSmear_mean_comp.pdf", "Cluster smear", 0,
    0, 2, 0, 0, 0, "Data", kBlue,
    1, 4, 0, 0, 0, "MC Photon", kMagenta+1);

  fillxj(0);
  fillxj(1);
  
  for (int i = 0; i < 2; i++) {
    

  (i == 0 ? cout << endl << endl << "xJ values with mean:" << endl : cout << endl << endl << "xJ values with fit:" << endl << endl);
  cout << "\\begin{center}" << endl;
  cout << "\\begin{tabular}{|| c | c | c | c | c | c ||}" << endl;
  cout << "\\hline" << endl;
  cout << "Jet Radius & Nominal $x_{J\\gamma}$ & Purity & Isolation & Jet topology & Model\\\\ [0.5ex]" << endl;
  cout << "\\hline\\hline" << endl;
  for (int ir = 0; ir < ana::nJetR; ir++) {
    float sys = (i == 0 ? TMath::Sqrt(xjmeb[ir]*xjmeb[ir] + xjme3[ir]*xjme3[ir] + xjmei[ir]*xjmei[ir] + xjmeh[ir]*xjmeh[ir]) : TMath::Sqrt(xjfeb[ir]*xjfeb[ir] + xjfe3[ir]*xjfe3[ir] + xjfei[ir]*xjfei[ir]+xjfeh[ir]*xjfeh[ir]));
    float xj = (i == 0 ? xjm[ir]   : xjf[ir]);
    float es = (i == 0 ? xjme[ir]  : xjfe[ir]);
    float eb = (i == 0 ? xjmeb[ir] : xjfeb[ir]);
    float ei = (i == 0 ? xjmei[ir] : xjfei[ir]);
    float e3 = (i == 0 ? xjme3[ir] : xjfe3[ir]);
    float eh = (i == 0 ? xjmeh[ir] : xjfeh[ir]);
    cout << std::defaultfloat << ana::JetRs[ir] << std::fixed << std::setprecision(4) << " & " << xj << " $\\pm$ " << es << " (stat) $\\pm$ " << sys << " (sys) & " << eb << " & " << ei << " & " << e3 << " & " << eh << " \\\\";
    if (ir == ana::nJetR -1) cout << " [1ex]";
    cout << endl;
  }
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\end{center}" << endl;
  } 
}
