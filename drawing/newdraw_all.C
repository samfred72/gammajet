#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

float xjm[ana::nJetR];
float xjme[ana::nJetR];
float xjmeb[ana::nJetR];
float xjme3[ana::nJetR];
float xjmei[ana::nJetR];
float xjmeh[ana::nJetR];
float xjmeJH[ana::nJetR];
float xjmeJL[ana::nJetR];
float xjf[ana::nJetR];
float xjfe[ana::nJetR];
float xjfeb[ana::nJetR];
float xjfe3[ana::nJetR];
float xjfei[ana::nJetR];
float xjfeh[ana::nJetR];
float xjfeJH[ana::nJetR];
float xjfeJL[ana::nJetR];

struct xjgroup {
  vector<vector<TH1D*>> hists =     vector<vector<TH1D*>>(ana::nPtBins, vector<TH1D*>(ana::nJetR));
  vector<vector<float>> histmeans = vector<vector<float>>(ana::nPtBins, vector<float>(ana::nJetR));
  vector<vector<float>> histerrs =  vector<vector<float>>(ana::nPtBins, vector<float>(ana::nJetR));
  vector<vector<TF1*>> fits =       vector<vector<TF1*>> (ana::nPtBins, vector<TF1*> (ana::nJetR));
  vector<vector<float>> fitmeans =  vector<vector<float>>(ana::nPtBins, vector<float>(ana::nJetR));
  vector<vector<float>> fiterrs =   vector<vector<float>>(ana::nPtBins, vector<float>(ana::nJetR));
  
  void clear() {
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      for (int ir = 0; ir < ana::nJetR; ir++) {
        delete hists[ipt][ir];
        delete fits[ipt][ir];
      }
    }
  }
};


struct xjgroup get_hists(int itype, int icalib, int iib, int i3jet, int iabcd, drawer &d, const char * sim = "pythia") {
  const char * histname = "hratio";
  struct xjgroup group;
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    for (int ir = 0; ir < ana::nJetR; ir++) {
      TH1D * hist = d.get(Form("%s_%i_%i_%i_%i_%i_%i",histname, ipt, ir, icalib, iib, i3jet, iabcd), itype);
      float mean = hist->GetMean();
      float meanerr = hist->GetMeanError();

      TF1 * func = d.fit(hist, 0, 2, "RMQI0");
      float funcmean = func->GetParameter(1);
      float funcerr = func->GetParError(1);
      int iformat = itype; 
      if (strcmp(sim,"herwig") == 0) iformat += 1; 
      d.format(hist,iformat);
      d.format(func,iformat);

      group.hists     [ipt][ir] = hist;
      group.histmeans [ipt][ir] = mean;
      group.histerrs  [ipt][ir] = meanerr;
      group.fits      [ipt][ir] = func;
      group.fitmeans  [ipt][ir] = funcmean;
      group.fiterrs   [ipt][ir] = funcerr;

    }
  }
  return group;
}
struct xjgroup get_short_hists(int itype, const char * histname, int iver, drawer &d, const char * sim = "pythia") {
  struct xjgroup group;
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    for (int ir = 0; ir < ana::nJetR; ir++) {
      TH1D * hist = d.get(Form("%s_%i_%i_%i",histname, ipt, ir, iver),itype);
      float mean = hist->GetMean();
      float meanerr = hist->GetMeanError();

      TF1 * func = d.fit(hist, 0, 2, "RMQI0");
      float funcmean = func->GetParameter(1);
      float funcerr = func->GetParError(1);
              
      if (strcmp(sim,"herwig") == 0) itype += 2; 
      d.format(hist,itype);
      d.format(func,itype);

      group.hists     [ipt][ir] = hist;
      group.histmeans [ipt][ir] = mean;
      group.histerrs  [ipt][ir] = meanerr;
      group.fits      [ipt][ir] = func;
      group.fitmeans  [ipt][ir] = funcmean;
      group.fiterrs   [ipt][ir] = funcerr;

    }
  }
  return group;
}

void draw_many(const char * cname, const char * info1, const char * info2, const char * variation,
    struct xjgroup g1, 
    struct xjgroup g2,
    drawer &d) {
  int csize = 700;
  TCanvas * c = new TCanvas("c","",csize*4,csize*3);
  c->SaveAs(Form("%s[",cname));

  float drawx = 0.77;
  float drawy = 0.92;
  float fontsize = 50;

  int nx = 3;
  int ny = 3;


  for (int ir = 0; ir < ana::nJetR; ir++) {
    int offset = 0;
    TPad * p = new TPad("p","",0,0,1,1);
    p->SetLeftMargin(0.25);
    p->SetBottomMargin(0.25);
    p->Divide(4,3,0,0);
    p->Draw();
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      TH1D * h1 = g1.hists[ipt][ir];
      TH1D * h2 = g2.hists[ipt][ir];
      h1->GetYaxis()->SetRangeUser(0,0.22);

      int index = ipt + 1 + offset;
      if (index % (nx+1) == 0) {index++; offset++;}
      p->cd(index);
      gPad->SetRightMargin(0.02);
      gPad->SetTicks(1,1);
      if ((ipt) < nx*(ny-1)) gPad->SetBottomMargin(0.04);
      gPad->SetTopMargin(0.01);
      if ((ipt) % nx != 0) gPad->SetLeftMargin(0.05);
      if (index == (ny*(nx+1) - 1)) h1->GetXaxis()->SetTitle("x_{J#gamma}"); // bottom right plot
      else h1->GetXaxis()->SetTitle("");
      h1->GetXaxis()->SetTitleSize(0.10);
      h1->GetXaxis()->SetLabelSize(0.08);
      if (index == (nx+1)*(ny-1)+1) h1->GetXaxis()->SetLabelSize(0.07); // bottom left plot
      if (ipt < nx*(ny-1) ) h1->GetXaxis()->SetLabelSize(0.00);
      h1->GetXaxis()->SetNdivisions(406);
      h1->GetYaxis()->SetNdivisions(406);
      //if (index != 1) h1->GetYaxis()->ChangeLabel(-1,-1,0);
      if (index != (ny*(nx+1) - 1)) h1->GetXaxis()->ChangeLabel(-1,-1,0);
      if (index == 1) h1->GetYaxis()->SetTitle("Normalized Counts"); // top left plot
      else h1->GetYaxis()->SetTitle("");
      h1->GetYaxis()->SetTitleSize(0.10);
      h1->GetYaxis()->SetLabelSize(0.08);
      if (ipt % nx != 0 ) h1->GetYaxis()->SetLabelSize(0.00);
      if (ipt == nx*(ny-1)) h1->GetYaxis()->SetLabelSize(0.06); // bottom left plot
      h1->GetYaxis()->SetLabelOffset(0.04);
      h1->GetYaxis()->SetMaxDigits(3);
      h1->GetYaxis()->SetDecimals(2);

      h1->Draw("hist e1 same");
      h2->Draw("hist same");

      float drawx = 0.15;
      float drawy_temp = 0.88;
      if ((index - 1) % 4 == 0) drawx = 0.35;
      if (ipt < 2*nx) drawy_temp = 0.83;
      d.drawMany({
          Form("#bf{%.0f GeV < p_{T}^{#gamma} < %.0f GeV}",ana::ptBins[ipt],ana::ptBins[ipt+1]),
          //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info1,  g1.histmeans[ipt][ir], g1.histerrs[ipt][ir]),
          //Form("#bf{#mu_{%s} = %0.3f #pm %0.3f}",info2,  g2.histmeans[ipt][ir], g2.histerrs[ipt][ir])
          },drawx,drawy_temp,42,c->GetWh()/3.0);
      TLine * mline1 = new TLine(g1.histmeans[ipt][ir],0,g1.histmeans[ipt][ir],0.15);
      mline1->SetLineColor(h1->GetLineColor());
      mline1->SetLineStyle(8);
      TLine * mline2 = new TLine(g2.histmeans[ipt][ir],0,g2.histmeans[ipt][ir],0.15);
      mline2->SetLineColor(h2->GetLineColor());
      mline2->SetLineStyle(7);
      mline1->Draw();
      mline2->Draw();

      if (ir == 2) {
        float m1 = g1.histmeans[ipt][ir];
        float m2 = g2.histmeans[ipt][ir];
        float e1 = g1.histerrs[ipt][ir];
        float e2 = g2.histerrs[ipt][ir];
        float err = (m1/m2) * sqrt( (e1/m1)*(e1/m1) + (e2/m2)*(e2/m2) );

        cout << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << std::fixed << std::setprecision(3) << " \\GeV & " 
          << m1 << " $\\pm$ " << e1 << " \\GeV & " 
          << m2 << " $\\pm$ " << e2 << " \\GeV & " 
          << m1/m2 << " $\\pm$ " << err << " \\\\ " << std::defaultfloat << endl;
      }
    }
    p->cd();
    if (strcmp(variation, "Nominal") == 0) variation = "";
    d.drawAll({
        //info1, 
        //info2,
        },
        {
        //"Analysis cuts",
        Form("Jet R=%1.1f",ana::JetRs[ir]),
        Form("#Delta#phi < %1.0f#pi/%1.0f",ana::oppnum,ana::oppden),
        variation,
        //calibstring,
        //t1[iabcd1].c_str()
        },
        drawx,drawy,fontsize,c->GetWh());
    TLine * mline1 = new TLine(0,0,1,1);
    mline1->SetLineColor(g1.hists[0][0]->GetLineColor());
    mline1->SetLineStyle(8);
    TLine * mline2 = new TLine(0,0,1,1);
    mline2->SetLineColor(g2.hists[0][0]->GetLineColor());
    mline2->SetLineStyle(7);
    TLegend * l2 = new TLegend(drawx,drawy-.3,0.99,drawy-.15);
    l2->SetLineWidth(0);
    l2->AddEntry(g1.hists[0][ir],info1);
    l2->AddEntry(g2.hists[0][ir],info2);
    l2->AddEntry(mline1,Form("Mean %s",info1),"l");
    l2->AddEntry(mline2,Form("Mean %s",info2),"l");
    l2->Draw();
    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
  g1.clear();
  g2.clear();
  delete c;
  return;
}


void draw_one(const char * cname, const char * info1, const char * info2, int ipt, int ir, drawer &d) {
  TCanvas * c = new TCanvas("c", "", 700, 700);

  float drawx = 0.6;
  float drawy = 0.85;
  string calibstring = "JES Calibrated";
  const char * t1 = "#bf{Analysis region A}";
  int colors[2] = {kBlue, kMagenta+1};
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.15);
  
  TH1D * h1 = d.get(Form("hratio_%i_%i_3_0_0_0",ipt,ir),0);
  h1->SetName("h1");
  TH1D * h2 = d.get(Form("hratio_%i_%i_3_0_0_0",ipt,ir),1);
  h2->SetName("h2");

  d.format(h1,0);
  d.format(h2,1);
  h1->GetXaxis()->SetTitle("x_{J#gamma}");
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

  TLine * line1 = new TLine(h1->GetMean(),0,h1->GetMean(),h1->GetMaximum());
  line1->SetLineStyle(9);
  line1->SetLineColor(kBlue);
  line1->Draw();
  TLine * line2 = new TLine(h2->GetMean(),0,h2->GetMean(),h1->GetMaximum());
  line2->SetLineStyle(9);
  line2->SetLineColor(kMagenta+1);
  line2->Draw();
  TLegend * l = new TLegend(.25,.75,.55,.87);
  l->SetLineWidth(0);
  l->SetFillStyle(0);
  l->AddEntry(h1,("data"));
  l->AddEntry(h2,("MC Photon"));
  l->AddEntry(line1,Form("mean %s",info1),"l");
  l->AddEntry(line2,Form("mean %s",info2),"l");
  l->Draw();
  d.drawAll({
      info1,
      info2,
      },
      {
      Form("%0.0f GeV < p_{T}^{#gamma} < %0.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
      "calibrated jet p_{T} > 3 GeV", 
      Form("Jet R=%0.1f",ana::JetRs[ir]), 
      Form("#Delta#phi > %1.0f#pi/%1.0f",ana::oppnum, ana::oppden),
      Form("#eta_{#gamma} < %1.1f",ana::etacut),
      Form("#eta_{Jet} < %1.1f",ana::etacut-ana::JetRs[ir]),
      calibstring,
      t1,
      },
      drawx,drawy,15,c->GetWh());
  c->SaveAs(cname);
  delete h1;
  delete h2;
  delete c;
  return;
}

vector<vector<float>> comp_axj(const char * cname, const char * info1, const char * info2, const char * variation,
    const char * fitfilestem,
    drawer &d) {
  int csize = 1000;
  TCanvas * c = new TCanvas("c","",csize,csize);
  c->SaveAs(Form("%s[",cname));

  vector<vector<float>> xj(8,vector<float>(3));
  
  for (int ir = 0; ir < ana::nJetR; ir++) {
    const char * fitfilename = Form("%s_%s_gammajet.root",fitfilestem,ana::rnames[ir]);
    TFile * infile = TFile::Open(fitfilename,"READ");
    
    TH1D * hMC = new TH1D("hMC",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);
    TH1D * hdata = new TH1D("hdata",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);
    TH1D * hcorrected = new TH1D("hcorrected",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      hMC->SetBinContent(ipt+1, ((TH1D*)infile->Get(Form("hxj_MC_%i",ipt)))->GetMean());
      hdata->SetBinContent(ipt+1, ((TH1D*)infile->Get(Form("hxj_data_%i",ipt)))->GetMean());
      hcorrected->SetBinContent(ipt+1, ((TH1D*)infile->Get(Form("hxj_corrected_%i",ipt)))->GetMean());
      hMC->SetBinError(ipt+1, ((TH1D*)infile->Get(Form("hxj_MC_%i",ipt)))->GetMeanError());
      hdata->SetBinError(ipt+1, ((TH1D*)infile->Get(Form("hxj_data_%i",ipt)))->GetMeanError());
      hcorrected->SetBinError(ipt+1, ((TH1D*)infile->Get(Form("hxj_corrected_%i",ipt)))->GetMeanError());
    }

    TPad * p1 = new TPad("p1","",0,.6,1,1);
    TPad * p2 = new TPad("p2","",0,0,1,0.6);
    p1->Draw();
    p2->Draw();
    p2->SetTopMargin(0);
    p1->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetTicks(1,1);
    gPad->SetBottomMargin(0.03);
    TLegend * lxj = new TLegend(.16,.65,.47,.85);
    lxj->SetLineWidth(0);



    hdata->SetMarkerColor(kBlue); hdata->SetLineColor(kBlue);
    hMC->SetMarkerColor(kMagenta+1); hMC->SetLineColor(kMagenta+1);
    if (strcmp(variation,"Herwig") == 0) {
      hMC->SetMarkerColor(kSpring-1); hMC->SetLineColor(kSpring-1);
    }

    
    hdata->GetYaxis()->SetRangeUser(0,2);
    hdata->GetYaxis()->SetTitleSize(0.09);
    hdata->GetYaxis()->SetLabelSize(0.08);
    hdata->GetYaxis()->SetTitleOffset(0.5);
    hdata->SetMarkerStyle(20);
    hdata->GetXaxis()->SetLabelSize(0);
    hMC->SetMarkerStyle(20);
    hdata->Draw("p");
    hMC->Draw("p same");
    lxj->AddEntry(hdata,info1);
    lxj->AddEntry(hMC,info2);
    lxj->Draw("same");
    //const char * fittext = (usefit ? "Fit method" : "Mean method");
    d.drawAll(
        {
        //info1,
        //info2
        },
        {
        //"Analysis cuts",
        //fittext,
        Form("Jet R=%0.1f",ana::JetRs[ir]), 
        },
        .5, .75, 30, c->GetWh()*0.7);

    TBox *texclude1 = new TBox(10,0,13.0,2);
    texclude1->SetFillColorAlpha(kGray,0.3);
    //texclude1->Draw("same");
    p2->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.02);
    gPad->SetTicks(1,1);
    TH1D * dummy = (TH1D*)hdata->Clone();
    dummy->SetName(Form("d%i",ir));
    dummy->Reset("ICES");

    dummy->GetYaxis()->SetTitle("Data / MC");
    dummy->GetYaxis()->SetTitleSize(0.07);
    dummy->GetYaxis()->SetLabelSize(0.05);
    dummy->GetYaxis()->SetTitleOffset(0.7);
    dummy->GetYaxis()->SetRangeUser(.88,1.08);

    dummy->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    dummy->GetXaxis()->SetLabelSize(0.06);
    dummy->GetXaxis()->SetTitleSize(0.06);
    dummy->GetXaxis()->SetTitleOffset(1);

    dummy->Draw();
    TH1D * hdivide_data = (TH1D*)hdata->Clone(Form("hdivide_data%i",ir));
    hdivide_data->Divide(hdata,hMC);
    hdivide_data->SetMarkerColor(kBlack);
    hdivide_data->SetLineColor(kBlack);
    hdivide_data->SetMarkerStyle(24);
    hdivide_data->SetMarkerSize(1);
    hdivide_data->Draw("same");
    TH1D * hdivide_corrected = (TH1D*)hdata->Clone(Form("hdivide_corrected%i",ir));
    hdivide_corrected->Divide(hcorrected,hMC);
    hdivide_corrected->SetMarkerColor(kBlack);
    hdivide_corrected->SetLineColor(kBlack);
    hdivide_corrected->SetMarkerStyle(20);
    hdivide_corrected->SetMarkerSize(1);
    hdivide_corrected->Draw("same");

    d.drawLine(ana::ptBins[0],1,ana::ptBins[ana::nPtBins],1); 
    
    TF1 * fline = (TF1*)infile->Get("fbest");
    fline->SetLineStyle(9);
    fline->SetLineWidth(3);
    //fline->Draw("same");
    
    TF1 * fHigh = (TF1*)infile->Get("fHigh");
    TF1 * fLow = (TF1*)infile->Get("fLow");
    
    float fbest = fline->GetParameter(0);
    float flow = fLow->Eval(0);
    float fhigh = fHigh->Eval(0);
    TBox *terr = new TBox(10,flow,35,fhigh);
    terr->SetFillColorAlpha(kRed,0.3);
    //terr->Draw("same");

    xj[ir][0] = flow-fbest;
    xj[ir][1] = fbest;
    xj[ir][2] = fhigh-fbest;

    float drawx = 0.20;
    float drawy = 0.82;
    if (strcmp(variation, "Nominal")) d.drawText(Form("Variation: %s",variation),drawx,drawy+0.07,kRed,40);
    d.drawText(Form("#bf{global #it{in situ} a = %2.4f ^{+%2.4f}_{-%2.4f}}",fbest, fhigh-fbest, fbest-flow),drawx,drawy,kRed,40);
    TLegend * ldivide = new TLegend(0.17, 0.24, .6, 0.44);
    ldivide->AddEntry(hdivide_data,"Uncorrected ratio");
    ldivide->AddEntry(hdivide_corrected,"Corrected ratio");
    ldivide->SetLineWidth(0);
    ldivide->Draw();

    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
  delete c;
  return xj;
}

vector<float> comp_comp_axj(const char * cname, const char * info1, const char * info2, const char * variation,
    const char * fitfilestem,
    drawer &d) {
  TCanvas * c = new TCanvas("c","",1000,600);
    
  vector<float> sys_xj(8);
  c->SaveAs(Form("%s[",cname));
  for (int ir = 0; ir < ana::nJetR; ir++) {
    const char * defaultfitfilename = Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_nominal_%s_gammajet.root",ana::rnames[ir]);
    const char * fitfilename = Form("%s_%s_gammajet.root",fitfilestem,ana::rnames[ir]);
    TFile * defaultinfile = TFile::Open(defaultfitfilename,"READ");
    TFile * infile = TFile::Open(fitfilename,"READ");

    TH1D * hdefaultMC = new TH1D("hdefaultMC",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);
    TH1D * hdefaultdata = new TH1D("hdefaultdata",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);
    TH1D * hMC = new TH1D("hMC",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);
    TH1D * hdata = new TH1D("hdata",";<x_{J#gamma}>;",ana::nPtBins,ana::ptBins);

    gPad->SetTicks(1,1);
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      hdefaultMC->SetBinContent(ipt+1,   ((TH1D*)defaultinfile->Get(Form("hxj_MC_%i",ipt)))->GetMean());
      hdefaultdata->SetBinContent(ipt+1, ((TH1D*)defaultinfile->Get(Form("hxj_data_%i",ipt)))->GetMean());
      hdefaultMC->SetBinError(ipt+1,     ((TH1D*)defaultinfile->Get(Form("hxj_MC_%i",ipt)))->GetMeanError());
      hdefaultdata->SetBinError(ipt+1,   ((TH1D*)defaultinfile->Get(Form("hxj_data_%i",ipt)))->GetMeanError());
      
      hMC->SetBinContent(ipt+1,   ((TH1D*)infile->Get(Form("hxj_MC_%i",ipt)))->GetMean());
      hdata->SetBinContent(ipt+1, ((TH1D*)infile->Get(Form("hxj_data_%i",ipt)))->GetMean());
      hMC->SetBinError(ipt+1,     ((TH1D*)infile->Get(Form("hxj_MC_%i",ipt)))->GetMeanError());
      hdata->SetBinError(ipt+1,   ((TH1D*)infile->Get(Form("hxj_data_%i",ipt)))->GetMeanError());
    }
    
    TH1D * dummy = (TH1D*)hdata->Clone(Form("d%i",ir));
    dummy->Reset("ICES");

    gPad->SetBottomMargin(.15);
    dummy->GetYaxis()->SetTitle("Ratio to Data");
    dummy->GetYaxis()->SetTitleSize(0.04);
    dummy->GetYaxis()->SetLabelSize(0.04);
    dummy->GetYaxis()->SetTitleOffset(1);
    dummy->GetYaxis()->SetRangeUser(.8,1.2);

    dummy->GetXaxis()->SetTitle("photon p_{T} [GeV]");
    dummy->GetXaxis()->SetLabelSize(0.04);
    dummy->GetXaxis()->SetTitleSize(0.04);
    dummy->GetXaxis()->SetTitleOffset(1);

    dummy->Draw();
    TH1D * hdivide1 = (TH1D*)hdata->Clone(Form("hdivide1%i",ir));
    TH1D * hdivide2 = (TH1D*)hMC->Clone(Form("hdivide2%i",ir));

    hdivide1->Divide(hdefaultdata,hdefaultMC);
    hdivide2->Divide(hdata,hMC);
    
    TGraphErrors * gdivide1 = new TGraphErrors();
    TGraphErrors * gdivide2 = new TGraphErrors();

    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      gdivide1->AddPoint((ana::ptBins[ipt] + ana::ptBins[ipt+1])/2.0, hdivide1->GetBinContent(ipt+1)); 
      gdivide2->AddPoint((ana::ptBins[ipt] + ana::ptBins[ipt+1])/2.0, hdivide2->GetBinContent(ipt+1)); 
      gdivide1->SetPointError(ipt, (ana::ptBins[ipt+1] - ana::ptBins[ipt])/2.0, hdivide1->GetBinError(ipt+1));
      gdivide2->SetPointError(ipt, (ana::ptBins[ipt+1] - ana::ptBins[ipt])/2.0, hdivide2->GetBinError(ipt+1));
    }
    
    gdivide1->SetLineColor(kBlack);
    gdivide1->SetMarkerColor(kBlack);
    gdivide1->SetMarkerStyle(20);
    
    gdivide2->SetLineColor(kBlue);
    gdivide2->SetMarkerColor(kBlue);
    gdivide2->SetMarkerStyle(20);

    gdivide1->Draw("p same");
    gdivide2->Draw("p same");

    d.drawLine(ana::ptBins[0],1,ana::ptBins[ana::nPtBins],1); 
    
    d.drawAll(
        {
        },
        {
        //"Analysis cuts",
        Form("Jet R=%0.1f",ana::JetRs[ir]), 
        },
        .17, .8, 30, c->GetWh());
    
    TF1 * defaultfline = (TF1*)defaultinfile->Get("fbest");
    TF1 * defaultfHigh = (TF1*)defaultinfile->Get("fHigh");
    TF1 * defaultfLow  = (TF1*)defaultinfile->Get("fLow");
    TF1 * fline = (TF1*)infile->Get("fbest");
    TF1 * fHigh = (TF1*)infile->Get("fHigh");
    TF1 * fLow  = (TF1*)infile->Get("fLow");

    float defaultfbest = defaultfline->GetParameter(0);
    float defaultfhigh = defaultfHigh->Eval(0) - defaultfbest;
    float defaultflow  = defaultfbest - defaultfLow->Eval(0);
    float fbest = fline->GetParameter(0);
    float fhigh = fHigh->Eval(0) - fbest;
    float flow  = fbest - fLow->Eval(0);

    float rbest = fbest / defaultfbest;
    float rhigh = (fbest + fhigh) / (defaultfbest - defaultflow) - rbest;
    float rlow = rbest - (fbest - flow) / (defaultfbest + defaultfhigh);

    if (strcmp(variation, "Herwig") == 0) {
      vector<float> vals = {0.0116,
        0.0127,
        0.0000,
        0.0069,
        0.0163,
        0.0210,
        0.0256};
      cout << ana::JetRs[ir] << fixed << setprecision(4) << " & $" << rbest << "^{+" << rhigh << "}_{-" << rlow << "}$ & " << vals[ir] << "\\\\" << std::defaultfloat << endl;
    } 
    
    TLine * l = new TLine(10,rbest, 35, rbest);
    l->SetLineWidth(3);
    l->SetLineColor(kBlue);
    l->SetLineStyle(9);
    //l->Draw("same");
    TBox *terr = new TBox(10,rlow,35,rhigh);
    terr->SetFillColorAlpha(kBlue,0.3);
    //terr->Draw("same");

    if (strcmp(variation, "Nominal")) d.drawText(Form("Variation: %s",variation),.16,.25,kBlue,40);
    d.drawText(Form("#bf{Relative Global #it{in situ} a = %2.4f ^{+%2.4f}_{-%2.4f}}",rbest, rhigh, rlow),.16,.18,kBlue,40);
    
    TLegend * leg = new TLegend(0.53,0.65,0.88,0.85);
    leg->AddEntry(gdivide1,"Default #it{in situ} ratio");
    leg->AddEntry(gdivide2,"Systematic #it{in situ} ratio");
    leg->SetLineWidth(0);
    leg->Draw();

    c->SaveAs(cname);
    c->Clear();
    sys_xj[ir] = rbest-1;
  }
  c->SaveAs(Form("%s]",cname));
  delete c;
  return sys_xj;
}

void newdraw_all() {
  gStyle->SetOptStat(0);
  if (!gROOT->IsBatch()) {
    cout << "Run with -b flag or else!!" << endl;
    gROOT->SetBatch(kTRUE);
    //return;
  }
  drawer dp = drawer("pythia");
  drawer dh = drawer("herwig");
      
  get_hists(1, 3, 0, 0, 0, dh, "herwig");
  
  int csize = 700;
  cout << "Drawing single bin..." << endl;
  int i = 7;
  int j = 2;
  draw_one(
      "/home/samson72/sphnx/gammajet/pdfs/note/xj_single_bin.pdf","p+p #sqrt{s}=200 GeV","Pythia 8 #gamma+jet", 
      i, j, dp
      );

  cout << "Drawing many bins..." << endl;
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Nominal",
      get_hists(0, 3, 0, 0, 0, dp, "pythia"),
      get_hists(1, 3, 0, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_Herwig.pdf", "p+p #sqrt{s}=200 GeV", "Herwig 3.7", "Herwig",
      get_hists(0, 3, 0, 0, 0, dp, "pythia"),
      get_hists(1, 3, 0, 0, 0, dh, "herwig"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_wideBDT.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Wider BDT Cut",
      get_hists(0, 3, 1, 0, 0, dp, "pythia"),
      get_hists(1, 3, 1, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_narrowISO.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Narrow Isolation Energy Cut",
      get_hists(0, 3, 2, 0, 0, dp, "pythia"),
      get_hists(1, 3, 2, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_3jet.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Third Jet Cut",
      get_hists(0, 3, 0, 1, 0, dp, "pythia"),
      get_hists(1, 3, 0, 1, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_JERhigh.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "High JER Smearing",
      get_hists(0, 5, 0, 0, 0, dp, "pythia"),
      get_hists(1, 5, 0, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_JERlow.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Low JER Smearing",
      get_hists(0, 6, 0, 0, 0, dp, "pythia"),
      get_hists(1, 6, 0, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_scalehigh.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "High EM Scale",
      get_hists(0, 1, 0, 0, 0, dp, "pythia"),
      get_hists(1, 1, 0, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_scalelow.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Low EM Scale",
      get_hists(0, 0, 0, 0, 0, dp, "pythia"),
      get_hists(1, 0, 0, 0, 0, dp, "pythia"),
      dp
      );
  draw_many(
      "/home/samson72/sphnx/gammajet/pdfs/note/allptbins_regionA_Jet.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 Jet", "",
      get_hists(0, 3, 1, 0, 0, dp, "pythia"),
      get_hists(2, 3, 1, 0, 0, dp, "pythia"),
      dp
      );
  
  cout << "Drawing mean xj..." << endl;
  vector<vector<float>> xj_vals;
  vector<vector<float>> xj_vals_herwig;
  xj_vals = comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Nominal",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_nominal",
      dp
      );
  xj_vals_herwig = comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_Herwig.pdf", "p+p #sqrt{s}=200 GeV", "Herwig3.7 #gamma+jet", "Herwig",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_HERWIG",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_wideBDT.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Wider BDT Cut",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_bdt",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_narrowISO.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Narrow Isolation Energy Cut",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_iso",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_3jet.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Third Jet Cut",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_3jet",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_JERhigh.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "High JER Smearing",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_JERhigh",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_JERlow.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Low JER Smearing",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_JERlow",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_scalehigh.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "High EM Scale",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_scalehigh",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_scalelow.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Low EM Scale",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_scalelow",
      dp
      );
  comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/axj_regionA_jet.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 Jet", "",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_nominal_jet",
      dp
      );


  cout << "Drawing <xj> systematics..." << endl;
  vector<float> sys_herwig = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_Herwig.pdf", "p+p #sqrt{s}=200 GeV", "Herwig3.7 #gamma+jet", "Herwig",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_HERWIG",
      dp
      );
  vector<float> sys_bdt = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_wideBDT.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Wider BDT Cut",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_bdt",
      dp
      );
  vector<float> sys_iso = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_narrowISO.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Narrow Isolation Energy Cut",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_iso",
      dp
      );
  vector<float> sys_3jet = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_3jet.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Third Jet Cut",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_3jet",
      dp
      );
  vector<float> sys_JERhigh = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_JERhigh.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "High JER Smearing",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_JERhigh",
      dp
      );
  vector<float> sys_JERlow = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_JERlow.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Low JER Smearing",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_JERlow",
      dp
      );
  vector<float> sys_scalehigh = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_scalehigh.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "High EM Scale",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_scalehigh",
      dp
      );
  vector<float> sys_scalelow = comp_comp_axj(
      "/home/samson72/sphnx/gammajet/pdfs/note/systematic_regionA_scalelow.pdf", "p+p #sqrt{s}=200 GeV", "Pythia 8 #gamma+jet", "Low EM Scale",
      "/home/samson72/sphnx/gammajet/hists/insitu_fit_scalelow",
      dp
      );

  for (int i = 0; i < 1; i++) { 
    (i == 0 ? cout << endl << endl << "xJ values with mean:" << endl : cout << endl << endl << "xJ values with fit:" << endl << endl);
    cout << "Jet Radius & Nominal $x_{J\\gamma}$ & Purity & Isolation & Jet topology & Model & JER & EM Scale\\\\ [0.5ex]" << endl;
    cout << "\\hline\\hline" << endl;
    for (int ir = 0; ir < ana::nJetR; ir++) {
      float herwig_sys = (xj_vals_herwig[ir][1] - xj_vals[ir][1])/2.0;
      
      float sys_down = (i == 0 ? TMath::Sqrt(
            //sys_herwig[ir]*sys_herwig[ir] + 
            herwig_sys*herwig_sys + 
            sys_bdt[ir]*sys_bdt[ir] +
            sys_iso[ir]*sys_iso[ir] +
            sys_3jet[ir]*sys_3jet[ir] +
            (sys_JERhigh[ir] < 0 ? sys_JERhigh[ir]*sys_JERhigh[ir] : 0) +
            (sys_scalelow[ir] < 0 ? sys_scalelow[ir]*sys_scalelow[ir] : 0)
            ) 
          : TMath::Sqrt(0));
      float sys_up = (i == 0 ? TMath::Sqrt(
            herwig_sys*herwig_sys +
            //sys_herwig[ir]*sys_herwig[ir] + 
            sys_bdt[ir]*sys_bdt[ir] +
            sys_iso[ir]*sys_iso[ir] +
            sys_3jet[ir]*sys_3jet[ir] +
            (sys_JERlow[ir] > 0 ? sys_JERlow[ir]*sys_JERlow[ir] : 0) +
            (sys_scalehigh[ir] > 0 ? sys_scalehigh[ir]*sys_scalehigh[ir] : 0)
            ) 
          : TMath::Sqrt(0));
      //float xj =  (i == 0 ? xj_vals[ir][1]  : 0); // THE result
      //float esu = (i == 0 ? fabs(xj_vals[ir][2])  : 0); // upper error
      //float esd = (i == 0 ? fabs(xj_vals[ir][0])  : 0); // lower error
      float xj =  (i == 0 ? (xj_vals[ir][1] + xj_vals_herwig[ir][1])/2.0  : 0); // THE result
      float esu = (i == 0 ? fabs(xj_vals[ir][2])  : 0); // upper error
      float esd = (i == 0 ? fabs(xj_vals[ir][0])  : 0); // lower error
      float eb =  (i == 0 ? fabs(sys_bdt[ir])     : 0);
      float ei =  (i == 0 ? fabs(sys_iso[ir])     : 0);
      float e3 =  (i == 0 ? fabs(sys_3jet[ir])    : 0);
      float eh =  (i == 0 ? fabs(herwig_sys)  : 0);
      //float eh =  (i == 0 ? fabs(sys_herwig[ir])  : 0);
      float eJH = (i == 0 ? (sys_JERhigh[ir] < 0 ? fabs(sys_JERhigh[ir]) : 0)        : 0); //JER is not abs'ded
      float eJL = (i == 0 ? (sys_JERlow[ir] > 0 ? fabs(sys_JERlow[ir]) : 0)          : 0);
      float esH = (i == 0 ? (sys_scalehigh[ir] > 0 ? fabs(sys_scalehigh[ir]) : 0)    : 0); //EM scale is not abs'ded
      float esL = (i == 0 ? (sys_scalelow[ir] < 0 ? fabs(sys_scalelow[ir]) : 0)      : 0);
      
      cout << std::defaultfloat << ana::JetRs[ir] << std::fixed << std::setprecision(4) << " & " << xj << "$^{+" << esu << "}_{-" << esd << "}$ (stat) $^{+" 
        << sys_up << "}_{-" << sys_down << "}$ (sys) & " 
        << eb << " & " << ei << " & " << e3 << " & " << eh 
        << " & $^{+" << eJL << "}_{-" << eJH << "}$"
        << " & $^{+" << esH << "}_{-" << esL << "}$ \\\\";
      if (ir == ana::nJetR -1) cout << " [1ex]";
      cout << endl;
    }
  }
}
