#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
drawer d;
TF1 * fitd[ana::nPtBins][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5];
TF1 * fitp[ana::nPtBins][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5]; 
TF1 * fitj[ana::nPtBins][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5];
const char * histname = "hratio";
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratiod = d.collect_hists(histname,0);
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratiop = d.collect_hists(histname,1);
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hratioj = d.collect_hists(histname,2);


void draw_wall(TCanvas * c, int iabcd, int k) {

  float drawx = 0.80;
  float drawy = 0.92;
  float fontsize = 60;
  int offset = 0;
  string calibstring = (k == 1 ? "JES Calibrated" : "Uncalibrated Jets");
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
      float minjet = (k == 1 ? ana::jet_calib_pt_cut[j] : ana::jet_pt_cut[j]);
      float drawx = 0.15;
      if ((index - 1) % 5 == 0) drawx = 0.35;
      d.drawMany({Form("%.0f GeV < p_{T}^{cluster} < %.0f GeV",lowcluster,highcluster),Form("Jet R=%1.1f",ana::JetRs[j])},drawx,0.85,1,42);
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

void draw_many(TCanvas * c, const char * cname, const char * info1, const char * info2,
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H1,          
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H2,          
    TF1 * f1[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], 
    TF1 * f2[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], 
    int ir1, int icalib1, int ibdt1, int i3jet1, int iabcd1,
    int ir2, int icalib2, int ibdt2, int i3jet2, int iabcd2) {

  float drawx = 0.80;
  float drawy = 0.92;
  float fontsize = 60;
  int offset = 0;
  int ir = 1;
  string calibstring = (icalib1 == 1 ? "JES Calibrated" : "Uncalibrated Jets");
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
      float minjet = (icalib1 == 1 ? ana::jet_calib_pt_cut[ir1] : ana::jet_pt_cut[ir1]);
      float drawx = 0.15;
      if ((index - 1) % 4 == 0) drawx = 0.35;
      d.drawMany({
          Form("#bf{%.0f GeV < p_{T}^{cluster} < %.0f GeV}",lowcluster,highcluster),
          Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info1,  f1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1),f1[ipt][ir1][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1)),
          Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info2,  f2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1),f2[ipt][ir2][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1))
          },drawx,0.88,42,c->GetWh()/3.0);
      TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,1);
      line->SetLineStyle(8);
      line->Draw();
  }
  p->cd();
  d.drawAll({
      info1, 
      info2},
      {
      "Analysis cuts",
      Form("Jet R=%1.1f",ana::JetRs[ir1]),
      calibstring,
      t1[iabcd1].c_str()},
      drawx,drawy,fontsize,c->GetWh());
  TLegend * l2 = new TLegend(drawx,drawy-.3,0.99,drawy-.2);
  l2->SetLineWidth(0);
  l2->AddEntry(H1[0][ir1][icalib1][ibdt1][i3jet1][iabcd1],info1);
  l2->AddEntry(H2[0][ir2][icalib2][ibdt2][i3jet2][iabcd2],info2);
  l2->Draw();
  return;
}

void draw_one(TCanvas * c, const char * cname, const char * info1, const char * info2,
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H1,          
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> H2,          
    TF1 * f1[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], 
    TF1 * f2[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], 
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
  float minjet = (icalib1 == 1 ? ana::jet_calib_pt_cut[ir1] : ana::jet_pt_cut[ir1]);
  TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,h1->GetMaximum());
  line->SetLineStyle(8);
  line->Draw();
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


void comp_axj(TCanvas * c, const char * cname, const char * info,
              TF1 * f1[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], int icalib1, int ibdt1, int i3jet1, int iabcd1, string info1, int color1, 
              TF1 * f2[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], int icalib2, int ibdt2, int i3jet2, int iabcd2, string info2, int color2) {
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
    TH1D * h1 = new TH1D(Form("h_%s_%s",info1.c_str(),f1[0][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetName()),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * h2 = new TH1D(Form("h_%s_%s",info2.c_str(),f1[0][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetName()),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      h1->SetBinContent(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1));
      h2->SetBinContent(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1));
      h1->SetBinError(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1));
      h2->SetBinError(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1));
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
    d.drawAll(
        {
        info1,
        info2
        },
        {
        "Analysis cuts",
        info,
        Form("Jet R=%0.1f",ana::JetRs[ir]), 
        },
        .5, .8, 20, c->GetWh()*0.6);

    TBox *texclude1 = new TBox(10,0,13.0,2);
    texclude1->SetFillColorAlpha(kGray,0.3);
    texclude1->Draw("same");
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
    
    TBox *texclude2 = new TBox(10,0.8,13.0,1.2);
    texclude2->SetFillColorAlpha(kGray,0.3);
    texclude2->Draw("same");

    TF1 * fline = new TF1(Form("func_%s",cname),"pol0",13,35);
    hdivide->Fit(fline,"RQIM0");
    fline->SetLineStyle(9);
    fline->SetLineWidth(3);
    fline->Draw("same");

    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
}

void comp_comp_axj(TCanvas * c, const char * cname, const char * variation,
              TF1 * f1[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], int icalib1, int ibdt1, int i3jet1, int iabcd1, string info1, int color1, 
              TF1 * f2[][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][5], int icalib2, int ibdt2, int i3jet2, int iabcd2, string info2, int color2) {
  c->SaveAs(Form("%s[",cname));
  for (int ir = 0; ir < ana::nJetR; ir++) {
    gPad->SetTicks(1,1);
    TH1D * h1 = new TH1D(Form("hc_%s_%s",info1.c_str(),f1[0][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetName()),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * h2 = new TH1D(Form("hc_%s_%s",info2.c_str(),f1[0][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetName()),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      h1->SetBinContent(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParameter(1));
      h2->SetBinContent(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParameter(1));
      h1->SetBinError(ipt+1,f1[ipt][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetParError(1));
      h2->SetBinError(ipt+1,f2[ipt][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetParError(1));
    }
    TH1D * h1c = new TH1D(Form("hcc_%s_%s",info1.c_str(),f1[0][ir][icalib1][ibdt1][i3jet1][iabcd1]->GetName()),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    TH1D * h2c = new TH1D(Form("hcc_%s_%s",info2.c_str(),f1[0][ir][icalib2][ibdt2][i3jet2][iabcd2]->GetName()),";;<x_{J}>",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      h1c->SetBinContent(ipt+1,fitd[ipt][ir][1][0][0][0]->GetParameter(1));
      h2c->SetBinContent(ipt+1,fitp[ipt][ir][1][0][0][0]->GetParameter(1));
      h1c->SetBinError(ipt+1,fitp[ipt][ir][1][0][0][0]->GetParError(1));
      h2c->SetBinError(ipt+1,fitd[ipt][ir][1][0][0][0]->GetParError(1));
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
        variation,
        },
        .5, .8, 20, c->GetWh());
    
    TBox *texclude = new TBox(10,0.8,13.0,1.2);
    texclude->SetFillColorAlpha(kGray,0.3);
    texclude->Draw("same");

    TF1 * func = new TF1(Form("func_%s",variation),"pol0",13,35);
    hdivide->Fit(func,"RQIM0");
    func->SetLineColor(kSpring+2);
    func->SetLineStyle(9);
    func->SetLineWidth(3);
    func->Draw("same");

    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
}

void draw_all() {
  gStyle->SetOptStat(0);
  if (!gROOT->IsBatch()) {
    cout << "Run with -b flag or else!!" << endl;
    return;
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
              if (n == 0 || n == 4) {
                fitd[i][j][k][l][m][n] = d.fit(hratiod[i][j][k][l][m][n], (int)(minjet/lowcluster/0.08 + 1)*0.08,2, "RLQI0");
                fitd[i][j][k][l][m][n]->SetParameter(0,fitd[i][j][k][l][m][n]->GetParameter(0)/hratiod[i][j][k][l][m][n]->GetEntries());
                d.format(fitd[i][j][k][l][m][n],0);
              }

              d.format(hratiod[i][j][k][l][m][n],0);
              if (k <= 1) d.format(hratiop[i][j][k][l][m][n],1);
              else d.format(hratiop[i][j][k][l][m][n],3);
              //d.format(hratioj[i][j][k][l][m][n],2);
              
              if (n == 0 || n == 4) {
                fitp[i][j][k][l][m][n] = d.fit(hratiop[i][j][k][l][m][n], (int)(minjet/lowcluster/0.08 + 1)*0.08,2,"RMQI0");
                //fitj[i][j][k][l][m][n] = d.fit(hratioj[i][j][k][l][m][n], (int)(minjet/lowcluster/0.04 + 1)*0.04,2);
                if (k <= 1) d.format(fitp[i][j][k][l][m][n],1);
                else d.format(fitp[i][j][k][l][m][n],3);
                 
                //d.format(fitj[i][j][k][l][m][n],2);
                //cout << fitd[i][j][k][l][m][n]->GetParameter(1) << endl;
              
                // JUST FOR TESTING!!!!
                fitd[i][j][k][l][m][n]->SetParameter(1,hratiod[i][j][k][l][m][n]->GetMean());
                fitd[i][j][k][l][m][n]->SetParError(1,hratiod[i][j][k][l][m][n]->GetMeanError());
                fitp[i][j][k][l][m][n]->SetParameter(1,hratiop[i][j][k][l][m][n]->GetMean());
                fitp[i][j][k][l][m][n]->SetParError(1,hratiop[i][j][k][l][m][n]->GetMeanError());
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
    draw_wall(c,iabcd,1);
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
    draw_one(cone_calib,cone_calib->GetName(),"Data","MC Photon", hratiod,hratiop,fitd,fitp,i,j,1,0,0,iabcd,i,j,1,0,0,iabcd);
    cone_calib->SaveAs("/home/samson72/sphnx/gammajet/pdfs/oneabcd_calib.pdf");
  }
  cone_calib->SaveAs("/home/samson72/sphnx/gammajet/pdfs/oneabcd_calib.pdf]");

  cout << "Drawing many bins..." << endl;
  TCanvas * cmany1 = new TCanvas("cmany1","",csize*4,csize*3);
  draw_many(
    cmany1, cmany1->GetName(), "Data", "MC Photon",
    hratiod, hratiop, fitd, fitp,      
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 0);
  cmany1->SaveAs("/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA.pdf");
  
  TCanvas * cmany2 = new TCanvas("cmany2","",csize*4,csize*3);
  draw_many(
    cmany2, cmany2->GetName(), "Data", "MC Photon",
    hratiod, hratiop, fitd, fitp,      
    1, 1, 1, 0, 0,
    1, 1, 1, 0, 0);
  cmany2->SaveAs("/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA_narrowBDT.pdf");
  
  TCanvas * cmany3 = new TCanvas("cmany3","",csize*4,csize*3);
  draw_many(
    cmany3, cmany3->GetName(), "Data", "MC Photon",
    hratiod, hratiop, fitd, fitp,      
    1, 1, 0, 1, 0,
    1, 1, 0, 1, 0);
  cmany3->SaveAs("/home/samson72/sphnx/gammajet/pdfs/allptbins_R04_regionA_3jetCut.pdf");
  
  TCanvas * cmany4 = new TCanvas("cmany4","",csize*4,csize*3);
  draw_many(
    cmany4, cmany4->GetName(), "MC Photon", "MC Photon smeared",
    hratiop, hratiop, fitp, fitp,      
    1, 1, 0, 0, 0,
    1, 2, 0, 0, 0);
  cmany4->SaveAs("/home/samson72/sphnx/gammajet/pdfs/allptbins_MC_R04_regionA_smear.pdf");

  cout << "Drawing <xj>..." << endl;
  TCanvas * caxj1 = new TCanvas("caxj1","",1000,1000);
  comp_axj(caxj1,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA.pdf","Nominal",
    fitd, 1, 0, 0, 0, "Data", kBlue,
    fitp, 1, 0, 0, 0, "MC Photon", kMagenta+1);
  TCanvas * caxj1comp = new TCanvas("caxj1comp","",1000,600);
  comp_comp_axj(caxj1comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_comp.pdf","Nominal",
    fitd, 1, 0, 0, 0, "Data", kBlue,
    fitp, 1, 0, 0, 0, "MC Photon", kMagenta+1);
  
  TCanvas * caxj2 = new TCanvas("caxj2","",1000,1000);
  comp_axj(caxj2,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_narrowBDT.pdf","Narrow BDT score",
    fitd, 1, 1, 0, 0, "Data", kBlue,
    fitp, 1, 1, 0, 0, "MC Photon", kMagenta+1);
  TCanvas * caxj2comp = new TCanvas("caxj2comp","",1000,600);
  comp_comp_axj(caxj2comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_narrowBDT_comp.pdf", "Narrow BDT score",
    fitd, 1, 1, 0, 0, "Data", kBlue,
    fitp, 1, 1, 0, 0, "MC Photon", kMagenta+1);
  
  TCanvas * caxj3 = new TCanvas("caxj3","",1000,1000);
  comp_axj(caxj2,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_3jetCut.pdf", "Third Jet Cut",
    fitd, 1, 0, 1, 0, "Data", kBlue,
    fitp, 1, 0, 1, 0, "MC Photon", kMagenta+1);
  TCanvas * caxj3comp = new TCanvas("caxj3comp","",1000,600);
  comp_comp_axj(caxj3comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_3jetCut_comp.pdf", "Third Jet cut",
    fitd, 1, 0, 1, 0, "Data", kBlue,
    fitp, 1, 0, 1, 0, "MC Photon", kMagenta+1);
  
  //TCanvas * caxj4comp = new TCanvas("caxj4comp","",1000,600);
  //comp_comp_axj(caxj4comp,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_comp.pdf","Nominal",
  //  fitd, 1, 0, 0, 0, "Data", kBlue,
  //  fitp, 1, 0, 0, 0, "MC Photon", kMagenta+1);
  
  
  TCanvas * caxj4 = new TCanvas("caxj4","",1000,1000);
  comp_axj(caxj4,"/home/samson72/sphnx/gammajet/pdfs/axj_regionA_smear.pdf", "JER smeared",
    fitp, 1, 0, 0, 0, "MC Photon", kMagenta+1,
    fitp, 2, 0, 0, 0, "MC Photon smeared", kMagenta+4);


}
