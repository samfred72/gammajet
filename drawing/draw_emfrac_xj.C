#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_emfrac_xj() {
  
  drawer d;
  TCanvas * c = new TCanvas("c","",1000,600);
  const char * cname = "/home/samson72/sphnx/gammajet/pdfs/xj_emfrac.pdf";
  c->SaveAs(Form("%s[",cname));
  int colors[ana::nEmfracBins] = {kRed,kGreen,kBlue};
  for (int ir = 0; ir < ana::nJetR; ir++) {
    gPad->SetLeftMargin(0.15);
    gPad->SetTicks(1,1);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    TLegend * lxj = new TLegend(.16,.65,.47,.85);
    lxj->SetLineWidth(0);
    
    TH1D * hd[3];
    TH1D * hp[3];
    TH1D * hxjd[3];
    TH1D * hxjp[3];
    TH1D * hrat[3];
    TF1 * fline[3];
    for (int iem = 0; iem < ana::nEmfracBins; iem++) {
      hxjd[iem] = new TH1D(Form("xjd%i_%i",iem,ir),";;<x_{J#gamma}>",ana::nPtBins,ana::ptBins);
      hxjp[iem] = new TH1D(Form("xjp%i_%i",iem,ir),";;<x_{J#gamma}>",ana::nPtBins,ana::ptBins);
      hrat[iem] = new TH1D(Form("ratio%i_%i",iem,ir),";p_{T}^{#gamma};<x_{J#gamma,data}>/<x_{J#gamma,MC}>",ana::nPtBins,ana::ptBins);
      fline[iem] = new TF1(Form("fline_%i_%i",iem,ir),"pol0",10,35);
      for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
        hd[iem] = d.get(Form("hratio_emfrac_%i_%i_%i",ipt,ir,iem),0);
        hp[iem] = d.get(Form("hratio_emfrac_%i_%i_%i",ipt,ir,iem),1);
        hxjd[iem]->SetBinContent(ipt + 1,hd[iem]->GetMean());
        hxjd[iem]->SetBinError(ipt + 1,hd[iem]->GetMeanError());
        hxjp[iem]->SetBinContent(ipt + 1,hp[iem]->GetMean());
        hxjp[iem]->SetBinError(ipt + 1,hp[iem]->GetMeanError());
      }
      hrat[iem]->Divide(hxjd[iem],hxjp[iem]);
      hrat[iem]->GetYaxis()->SetRangeUser(0.85,1.3);
      hrat[iem]->SetMarkerColor(colors[iem]);
      hrat[iem]->SetMarkerSize(1);
      hrat[iem]->SetMarkerStyle(20);
      hrat[iem]->SetLineColor(colors[iem]);
      hrat[iem]->Fit(fline[iem],"RIQM0");
      fline[iem]->SetLineColor(colors[iem]);


      float flow = fline[iem]->GetParameter(0) - fline[iem]->GetParError(0);
      float fhigh = fline[iem]->GetParameter(0) + fline[iem]->GetParError(0);
      TBox *terr = new TBox(10,flow,35,fhigh);
      terr->SetFillColorAlpha(colors[iem],0.3);
      
      lxj->AddEntry(hrat[iem],Form("%0.1f < EM frac < %0.1f",ana::emfracBins[iem],ana::emfracBins[iem+1]));
      
      if (iem == 0) hrat[iem]->Draw("p");
      else hrat[iem]->Draw("p same");
      fline[iem]->Draw("same");
      terr->Draw("same");
    }
    lxj->Draw();
    d.drawAll({"p+p #sqrt{s}=200 GeV","Pythia 8 #gamma+jet"},{Form("Jet R=%0.1f",ana::JetRs[ir])},.55,.80,20,600);

    c->SaveAs(cname);
    c->Clear();
  }
  c->SaveAs(Form("%s]",cname));
}

