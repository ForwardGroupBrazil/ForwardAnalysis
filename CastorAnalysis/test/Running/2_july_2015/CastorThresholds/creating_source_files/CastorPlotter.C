/*

Program: Castor Analysis Plotter
Goal: Create histograms.
Project: Diffractive Z Analysis. Castor Studies.

How to Execute:

> root -l CastorPlotter.C

No defense for empyt histograms or files. Private use.

 */

#include <iostream>
#include <fstream>


void CastorPlotter(){
  //StyleNice();
  //MakeMultipleSingle("muon");
  //ThresholdsStudiesChannels("histo_Castor_threshold_p1_unpaired.root","unpaired",0,1);
  //ThresholdsStudiesSectors("histo_Castor_threshold_p2_unpaired.root","unpaired",0,1);
  //ThresholdsStudiesSectors("histo_Castor_threshold_p3_unpaired.root","unpaired",0,1);
  CalculateCorrector("CastorMappingEnergySelection0_step7","ratio_all.root");
  //CalculateCorrector("CastorMappingEnergySelection1_step7","ratio_step1.root");
  //CalculateCorrector("CastorMappingEnergySelection2_step7","ratio_step2.root");
  //CalculateCorrector("CastorMappingEnergySelection3_step7","ratio_step3.root");
  //CalculateCorrector("CastorMappingEnergySelection4_step7","ratio_step4.root");
  //CalculateCorrector("CastorMappingEnergySelection5_step7","ratio_step5.root");
  //CalculateCorrector("CastorMappingEnergySelection6_step7","ratio_step6.root");
  //CalculateCorrectorFlow();

}

//
// Create a list of objects from a rootfile
//
// Output: ListOfHistograms.txt

void RunAll(TString filename, TString logOpt, TString autonorma, TString type, TString th){

  std::ofstream outstring("ListOfHistograms.txt");

  gStyle->SetOptStat(0);
  TFile *f1= new TFile(filename);

  TList* tdf_obs = f1->GetListOfKeys();
  char *hName1 = new char [80];
  char *filesave = new char[80];

  TList *objlist = f1->GetListOfKeys();
  TIterator *itr = objlist->MakeIterator();
  TKey *key;

  TString *currentdir;
  TIter iter(tdf_obs);

  outstring << "" << std::endl;
  outstring << "<< List of Objects >>" << std::endl;
  outstring << " " << std::endl;

  while ((key = static_cast<TKey*>(iter.Next())) != 0) {

    TObject *obj = key->ReadObj();
    TH1F *h = (TH1F*)key->ReadObj();
    sprintf(filesave,"%s.png",h->GetName());

    outstring << obj->GetName() << std::endl;
    cout << "This object has name " << obj->GetName() << endl;

    if (obj->IsA()->InheritsFrom("TDirectory")) {
      TDirectory *d = static_cast<TDirectory*>(obj);
      cout << "Dir: " << d->GetPath() << endl;
    }

    if (h->GetEntries() != 0){ // Defense
      if(obj->IsA()->InheritsFrom("TH1")){
	hName1=h->GetName();
      }
    }

  }

  outstring.close();

}

//
// Castor Phi Plots
//

void SectorPhi(TString name1){

  TCanvas *c1 = new TCanvas(name1,name1);
  c1->cd();

  gStyle->SetOptStat(0);

  TLegend* leg = new TLegend(0.7597956,0.822335,0.9931857,0.9949239,NULL,"brNDC");
  TString legdata, legmc;

  TFile *l1 = TFile::Open("histo_Muon_Reco.root");
  TFile *l2 = TFile::Open("histo_DyToMuMu_Reco.root");
  legdata = "2010 data";
  legmc = "MC";

  TH1F* h_1 = (TH1F*)l1->Get(name1);
  TH1F* h_2 = (TH1F*)l2->Get(name1);

  double weight1, weight2;

  weight1 = 1./h_1->GetEntries();
  weight2 = 1./h_2->GetEntries();

  h_1->SetLineColor(kBlack);
  h_1->SetMarkerStyle(20);
  h_1->SetMarkerSize(0.8);
  h_1->Scale(weight1*100);
  h_1->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_1,legdata,"p");

  h_2->SetLineColor(kRed);
  h_2->SetLineWidth(2);
  h_2->Scale(weight2*100);
  h_2->SetFillColor(kRed-4);
  h_2->SetFillStyle(3020);
  h_2->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_2,legmc,"LFP");

  h_1->Draw("histo");
  h_2->Draw("histosame");
  leg->Draw();

}

//
// Castor: Each sector energy distribution. Each sector energy Vs Total Castor Multiplicity Plots
//

void MakeMultipleSingle(TString type){

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("multiple1","multiple1",2000,1000);
  c1->Divide(4,4);

  TCanvas *c2 = new TCanvas("multiple2","multiple2",2000,1000);
  c2->Divide(4,4);

  TCanvas *c3 = new TCanvas("multiple3","multiple3",2000,1000);
  c3->Divide(4,4);

  TCanvas *c4 = new TCanvas("multiple4","multiple4",500,500);

  TCanvas *c5 = new TCanvas("multiple5","multiple5",2000,250);
  c5->Divide(5,1);

  TCanvas *c6 = new TCanvas("multiple6","multiple6",500,500);

  TCanvas *c7 = new TCanvas("multiple7","multiple7",2000,1000);
  c7->Divide(4,4);

  TCanvas *c8 = new TCanvas("multiple8","multiple8",500,500);

  TCanvas *c9 = new TCanvas("multiple9","multiple9",500,500);

  TCanvas *c10 = new TCanvas("multiple10","multiple10",1000,500);
  c10->Divide(2,1);

  TCanvas *c11 = new TCanvas("multiple11","multiple11",1000,2000);
  c11->Divide(5,5);

  TCanvas *c12 = new TCanvas("multiple12","multiple12",2000,1000);
  c12->Divide(5,5);


  TLegend* leg = new TLegend(0.7597956,0.822335,0.9931857,0.9949239,NULL,"brNDC");

  if (type == "muon" || type == "Muon" || type == "MUON"){
    //TFile *l1 = TFile::Open("histo_Castor_muon.root");
    TFile *l1 = TFile::Open("histo_DyToMuMu_RECO.root");
    TFile *l1mc = TFile::Open("histo_DyToMuMu_RECO.root");
  }
  else if (type == "electron" || type == "Electron" || type == "ELECTRON"){
    TFile *l1 = TFile::Open("histo_Castor_electron.root");
  }
  else {
    std::cout << "Please, put correct option." << std::endl;
    exit(0);
  }

  for (int i=1; i< 17; i++){
    char name[300];
    sprintf(name,"TotalEnergySector%d_step7",i);
    TH1F* h_1 = (TH1F*)l1->Get(name);
    //TH1F* h_2 = h_1->Clone();
    TH1F* h_3 = (TH1F*)l1mc->Get(name);
    //h_2->Reset();
    //h_2->SetMaximum(1.);
    c1->cd(i);
    //gPad->SetLogx(1);
    gPad->SetLogy(1);
    //h_2->Draw();
    //h_1->SetMarkerStyle(8);
    //h_1->SetMarkerSize(1.1);
    h_3->SetLineColor(kRed);
    h_3->SetLineWidth(1.1);
    h_1->GetYaxis()->SetTitleOffset(1.4);
    h_1->GetXaxis()->SetRangeUser(0.,10.);
    h_3->GetXaxis()->SetRangeUser(0.,10.);
    h_1->DrawNormalized("histo");
    h_3->DrawNormalized("histosame");
    h_1->SetMaximum(1.);
  }

  for (int i=1; i< 17; i++){
    gStyle->SetPalette(1);
    char name2[300];
    sprintf(name2,"Sector%d_EnergyVsMultiplicity_step7",i);
    TH2F* h_1b = (TH2F*)l1->Get(name2);
    TH2F* h_2b = h_1b->Clone();
    h_2b->Reset();
    c2->cd(i);
    h_1b->RebinY(3);
    h_2b->Draw("colz");
    h_1b->GetYaxis()->SetTitleOffset(1.4);
    h_1b->Draw("colzsame");
  }

  for (int i=1; i< 17; i++){
    gStyle->SetPalette(1);
    char name3[300];
    sprintf(name3,"AlongZ_Sector%d_EnergyVsModule_step7",i);
    TH2F* h_1c = (TH2F*)l1->Get(name3);
    TH2F* h_2c = h_1c->Clone();
    h_2c->Reset();
    c3->cd(i);
    h_1c->GetYaxis()->SetRangeUser(0.,400.);
    h_2c->GetYaxis()->SetRangeUser(0.,400.);
    h_1c->RebinY(3);
    h_2c->Draw("colz");
    h_1c->GetYaxis()->SetTitleOffset(1.4);
    h_1c->Draw("colzsame");
  }

  for (int i=1; i< 6; i++){
    char name4[300];
    sprintf(name4,"CastorMultiplicityModule%d_step7",i);
    TH1F* h_1d = (TH1F*)l1->Get(name4);
    c4->cd();
    h_1d->GetYaxis()->SetTitleOffset(1.4);
    h_1d->SetMarkerColor(i);
    h_1->SetMarkerSize(1.5);
    h_1d->SetMarkerStyle(8);
    h_1d->GetYaxis()->SetRangeUser(0.,15000.);
    h_1d->SetTitle("Multiplicity");
    char legtitle[300];
    sprintf(legtitle,"Module %d",i);
    leg->AddEntry(h_1d,legtitle,"p");
    if (i==1) {
      h_1d->Draw("p");
    }else{
      h_1d->Draw("psame");
    }
    leg->Draw();

    char name5[300];
    sprintf(name5,"CastorMultiplicityModule%dVsTotal_step7",i);
    TH2F* h_1e = (TH2F*)l1->Get(name5);
    c5->cd(i);
    gStyle->SetPalette(1);
    h_1e->GetYaxis()->SetTitleOffset(1.4);
    h_1e->Draw("colz");
  }

  TH1F* h_1f = (TH1F*)l1->Get("CastorMultiplicity_step7");
  c6->cd();
  h_1f->GetYaxis()->SetTitleOffset(1.4);
  h_1f->Draw();

  for (int i=1; i< 17; i++){
    gStyle->SetPalette(1);
    char name6[300];
    sprintf(name6,"AlongZ_Sector%d_EnergyVsModuleTProf_step7",i);
    TProfile* h_1g = (TProfile*)l1->Get(name6);
    TProfile* h_1gmc = (TProfile*)l1mc->Get(name6);
    c7->cd(i);
    h_1g->SetMarkerStyle(8);
    h_1g->SetMarkerSize(1.1);
    h_1g->GetYaxis()->SetTitleOffset(1.4);
    h_1gmc->SetMarkerStyle(8);
    h_1gmc->SetMarkerColor(kRed);
    h_1gmc->SetMarkerSize(1.1);
    h_1gmc->GetYaxis()->SetTitleOffset(1.4);
    h_1g->Draw("ep");
    h_1gmc->GetYaxis()->SetRangeUser(0.,35.);
    h_1g->GetYaxis()->SetRangeUser(0.,35.);
    h_1gmc->Draw("epsame");
  }

  TProfile* h_1h = (TProfile*)l1->Get("MultiplicityPerModuleTProf_step7");
  c8->cd();
  h_1h->GetYaxis()->SetTitleOffset(1.4);
  h_1h->Draw();

  TH2F* h_1i = (TH2F*)l1->Get("MultiplicityPerModule_step7");
  c9->cd();
  gStyle->SetPalette(1);
  h_1i->GetYaxis()->SetTitleOffset(1.4);
  h_1i->Draw("colz");

  TH2F* h_1j = (TH2F*)l1->Get("CastorMappingMultiplicities_step7");
  c10->cd(1);
  gStyle->SetPalette(1);
  h_1j->GetYaxis()->SetTitle("Sector (#phi)");
  h_1j->GetYaxis()->SetTitleOffset(1.4);
  h_1j->Draw("colz");

  TH2F* h_1k = (TH2F*)l1->Get("CastorMappingEnergy_step7");
  c10->cd(2);
  gStyle->SetPalette(1);
  h_1k->GetYaxis()->SetTitle("Sector (#phi)");
  h_1k->GetYaxis()->SetTitleOffset(1.4);
  h_1k->Draw("colz");


  //TFile *ld = TFile::Open("histo_Castor_muonp2.root");
  TFile *ld = TFile::Open("histo_DyToMuMu_RECO.root");
  int display1;
  for (int i=1; i< 6; i++){
    for (int j=1; j< 6;j++){
      display1++;
      gStyle->SetPalette(1);
      char name7[300];
      sprintf(name7,"CastorMappingMultiplicities_module_%d_snapshot_%d",i,j);
      TH2F* h_1l = (TH2F*)ld->Get(name7);
      c11->cd(display1);
      h_1l->SetMarkerStyle(8);
      h_1l->SetMarkerSize(1.1);
      h_1l->GetYaxis()->SetTitleOffset(1.4);
      h_1l->Draw("colz");
    }
  }

  int display2;
  for (int i=1; i< 6; i++){
    for (int j=1; j< 6;j++){
      display2++;
      gStyle->SetPalette(1);
      char name8[300];
      sprintf(name8,"CastorMappingEnergy_module_%d_snapshot_%d",i,j);
      TH2F* h_1m = (TH2F*)ld->Get(name8);
      c12->cd(display2);
      h_1m->SetMarkerStyle(8);
      h_1m->SetMarkerSize(1.1);
      h_1m->GetYaxis()->SetTitleOffset(1.4);
      h_1m->Draw("colz");
    }
  }

}

//
// Castor: sector thresholds. Thresholds Measurements.
//
// Output: file with each sector threshold.

void ThresholdsStudiesSectors(TString openfile, TString complement, bool printPar, bool printTh){

  gStyle->SetOptStat("em");

  TCanvas *c1 = new TCanvas("sectorenergy","sectorenergy",2000,1000);
  c1->Divide(4,4);
  TCanvas *c2 = new TCanvas("thresholds","thresholds",500,500);
  TCanvas *c3 = new TCanvas("sigma","sigma",500,500);

  TFile *l1 = TFile::Open(openfile);

  TString outtxt = openfile;
  outtxt.ReplaceAll("root","txt");
  TString newname = "FitResults_"+outtxt;
  std::ofstream outstring(newname);
  Double_t x[16], y[16];
  Double_t ex[16], ey[16];
  Double_t ys[16];
  Double_t eys[16];

  for (int i=1; i< 17; i++){

    char name[300];
    char text[300];
    char textth[300];
    char textthf[300];

    sprintf(name,"Sector%d_CastorSumEnergy_with_type_",i);
    TString finalname = name + complement;
    cout << finalname << endl;
    TH1F* h_1 = (TH1F*)l1->Get(finalname);
    c1->cd(i);
    h_1->Fit("gaus","","",-5.,5.);
    h_1->GetXaxis()->SetRangeUser(-5*h_1->GetFunction("gaus")->GetParameter(2),5*h_1->GetFunction("gaus")->GetParameter(2));
    h_1->GetYaxis()->SetRangeUser(0.,1.3*h_1->GetFunction("gaus")->GetMaximum());
    h_1->GetYaxis()->SetTitleOffset(1.4);
    h_1->Draw();

    double limit = 4*h_1->GetFunction("gaus")->GetParameter(2);
    Double_t maxvalue = 0.9*h_1->GetMaximum();

    TLine *line = new TLine(limit,0,limit,0.7*maxvalue);
    line->SetLineColor(8);
    line->SetLineWidth(3);
    line->Draw();

    outstring << "\n<< Sector " << i << " >>" << endl;
    outstring << "Mean: " << h_1->GetFunction("gaus")->GetParameter(1)<< " GeV" << endl;
    outstring << "Sigma: " << h_1->GetFunction("gaus")->GetParameter(2)<< " GeV" << endl;
    outstring << "Threshold (4xsigma): " << 4*h_1->GetFunction("gaus")->GetParameter(2)<< " GeV" << endl;
    outstring << "Error Threshold (4xsigma): " << h_1->GetFunction("gaus")->GetParError(2) << " GeV" << endl;

    y[i-1] = 4*h_1->GetFunction("gaus")->GetParameter(2);
    x[i-1] = i;
    ex[i-1] = 0;
    ey[i-1] = 4*h_1->GetFunction("gaus")->GetParError(2);

    ys[i-1] = h_1->GetFunction("gaus")->GetParameter(2);
    eys[i-1] = h_1->GetFunction("gaus")->GetParError(2);

    TLatex *lt = new TLatex(0.7,160,text);
    sprintf(text,"f(x) = %g*e^{#frac{x-%g}{%g}}",h_1->GetFunction("gaus")->GetParameter(1),h_1->GetFunction("gaus")->GetParameter(2),h_1->GetFunction("gaus")->GetParameter(3));
    lt->SetTextSize(0.1);
    lt->SetTextFont(72);
    lt->SetTextColor(kBlue);
    if (printPar) lt->DrawLatex(2.,1.,text);

    TLatex *ltt = new TLatex(0.7,160,textth);
    sprintf(textth,"Threshold: %g GeV",4*h_1->GetFunction("gaus")->GetParameter(2));
    ltt->SetTextSize(0.08);
    ltt->SetTextFont(72);
    ltt->SetTextColor(kBlack);
    if (printTh) ltt->DrawLatex(-3*h_1->GetFunction("gaus")->GetParameter(2),1.1*h_1->GetFunction("gaus")->GetMaximum(),textth);

  }

  c2->cd();
  gr = new TGraphErrors(16,x,y,ex,ey);
  gr->SetTitle("Castor Threshold per Sector");
  gr->GetYaxis()->SetTitle("Threshold [GeV]");
  gr->GetXaxis()->SetTitle("Sector");
  gr->SetMarkerSize(1.);
  gr->GetYaxis()->SetTitleOffset(1.4);
  gr->SetMarkerStyle(21);
  gr->GetYaxis()->SetRangeUser(0.,1.5*y[0]);
  gr->Draw("AP");

  c3->cd();
  grs = new TGraphErrors(16,x,ys,ex,eys);
  grs->SetTitle("#sigma_{fit} per Sector");
  grs->GetYaxis()->SetTitle("#sigma_{fit} [GeV]");
  grs->GetXaxis()->SetTitle("Sector");
  grs->SetMarkerSize(1.);
  grs->SetMarkerStyle(21);
  grs->GetYaxis()->SetTitleOffset(1.4);
  grs->GetYaxis()->SetRangeUser(0.,1.5*ys[0]);
  grs->Fit("pol0");
  grs->Draw("AP");

  TLatex *lttf = new TLatex(0.7,160,textthf);
  sprintf(textthf,"#sigma_{fit}: %g GeV, Threshold: %g GeV",grs->GetFunction("pol0")->GetParameter(0),4*grs->GetFunction("pol0")->GetParameter(0));
  lttf->SetTextSize(0.03);
  lttf->SetTextFont(72);
  lttf->SetTextColor(kBlack);
  lttf->DrawLatex(1,1.35*ys[0],textthf);

  outstring.close();

}

//
// Castor: channels thresholds. Thresholds Measurements.
//
// Output: file with each channel threshold.

void ThresholdsStudiesChannels(TString openfile, TString complement, bool printPar, bool printTh){

  gStyle->SetOptStat("em");

  TCanvas *c1 = new TCanvas("channelenergy1","channelenergy1",2000,1000);
  c1->Divide(4,4);

  TCanvas *c2 = new TCanvas("channelenergy2","channelenergy2",2000,1000);
  c2->Divide(4,4);

  TCanvas *c3 = new TCanvas("channelenergy3","channelenergy3",2000,1000);
  c3->Divide(4,4);

  TCanvas *c4 = new TCanvas("channelenergy4","channelenergy4",2000,1000);
  c4->Divide(4,4);

  TCanvas *c5 = new TCanvas("channelenergy5","channelenergy5",2000,1000);
  c5->Divide(4,4);

  TCanvas *c6 = new TCanvas("thresholds","thresholds",500,500);
  TCanvas *c7 = new TCanvas("sigma","sigma",500,500);

  TFile *l1 = TFile::Open(openfile);

  TString outtxt = openfile;
  outtxt.ReplaceAll("root","txt");
  TString newname = "FitResults_"+outtxt;
  std::ofstream outstring(newname);
  Double_t x[80], y[80];
  Double_t ex[80], ey[80];
  Double_t ys[80];
  Double_t eys[80];

  for (int i=1; i< 81; i++){

    char name[300];
    char text[300];
    char textth[300];
    char textthf[300];

    sprintf(name,"Channel%d_Energy_with_type_",i);
    TString finalname = name + complement;
    cout << finalname << endl;
    TH1F* h_1 = (TH1F*)l1->Get(finalname);
    if (i>=1 && i<=16) c1->cd(i);
    if (i>=17 && i<=32) c2->cd(i-16);
    if (i>=33 && i<=48) c3->cd(i-32);
    if (i>=49 && i<=64) c4->cd(i-48);
    if (i>=65 && i<=80) c5->cd(i-64);
    h_1->Fit("gaus","","",-5.,5.);
    h_1->GetXaxis()->SetRangeUser(-5*h_1->GetFunction("gaus")->GetParameter(2),5*h_1->GetFunction("gaus")->GetParameter(2));
    h_1->GetYaxis()->SetRangeUser(0.,1.3*h_1->GetFunction("gaus")->GetMaximum());
    h_1->GetYaxis()->SetTitleOffset(1.4);
    h_1->Draw();

    double limit = 4*h_1->GetFunction("gaus")->GetParameter(2);
    Double_t maxvalue = 0.9*h_1->GetMaximum();

    TLine *line = new TLine(limit,0,limit,0.7*maxvalue);
    line->SetLineColor(8);
    line->SetLineWidth(3);
    line->Draw();

    outstring << "\n<< Sector " << i << " >>" << endl;
    outstring << "Mean: " << h_1->GetFunction("gaus")->GetParameter(1)<< " GeV" << endl;
    outstring << "Sigma: " << h_1->GetFunction("gaus")->GetParameter(2)<< " GeV" << endl;
    outstring << "Threshold (4xsigma): " << 4*h_1->GetFunction("gaus")->GetParameter(2)<< " GeV" << endl;
    outstring << "Error Threshold (4xsigma): " << h_1->GetFunction("gaus")->GetParError(2) << " GeV" << endl;

    y[i-1] = 4*h_1->GetFunction("gaus")->GetParameter(2);
    x[i-1] = i;
    ex[i-1] = 0;
    ey[i-1] = 4*h_1->GetFunction("gaus")->GetParError(2);

    ys[i-1] = h_1->GetFunction("gaus")->GetParameter(2);
    eys[i-1] = h_1->GetFunction("gaus")->GetParError(2);

    TLatex *lt = new TLatex(0.7,160,text);
    sprintf(text,"f(x) = %g*e^{#frac{x-%g}{%g}}",h_1->GetFunction("gaus")->GetParameter(1),h_1->GetFunction("gaus")->GetParameter(2),h_1->GetFunction("gaus")->GetParameter(3));
    lt->SetTextSize(0.1);
    lt->SetTextFont(72);
    lt->SetTextColor(kBlue);
    if (printPar) lt->DrawLatex(2.,1.,text);

    TLatex *ltt = new TLatex(0.7,160,textth);
    sprintf(textth,"Threshold: %g GeV",4*h_1->GetFunction("gaus")->GetParameter(2));
    ltt->SetTextSize(0.08);
    ltt->SetTextFont(72);
    ltt->SetTextColor(kBlack);
    if (printTh) ltt->DrawLatex(-3*h_1->GetFunction("gaus")->GetParameter(2),1.1*h_1->GetFunction("gaus")->GetMaximum(),textth);

  }

  c6->cd();
  gr = new TGraphErrors(80,x,y,ex,ey);
  gr->SetTitle("Castor Threshold per Channel");
  gr->GetYaxis()->SetTitle("Threshold [GeV]");
  gr->GetXaxis()->SetTitle("Channel");
  gr->SetMarkerSize(1.);
  gr->GetYaxis()->SetTitleOffset(1.4);
  gr->SetMarkerStyle(21);
  gr->GetYaxis()->SetRangeUser(0.,1.5*y[0]);
  gr->Draw("AP");

  c7->cd();
  grs = new TGraphErrors(80,x,ys,ex,eys);
  grs->SetTitle("#sigma_{fit} per Channel");
  grs->GetYaxis()->SetTitle("#sigma_{fit} [GeV]");
  grs->GetXaxis()->SetTitle("Channel");
  grs->SetMarkerSize(1.);
  grs->SetMarkerStyle(21);
  grs->GetYaxis()->SetTitleOffset(1.4);
  grs->GetYaxis()->SetRangeUser(0.,1.5*ys[0]);
  grs->Fit("pol0");
  grs->Draw("AP");

  TLatex *lttf = new TLatex(0.7,160,textthf);
  sprintf(textthf,"#sigma_{fit}: %g GeV, Threshold: %g GeV",grs->GetFunction("pol0")->GetParameter(0),4*grs->GetFunction("pol0")->GetParameter(0));
  lttf->SetTextSize(0.03);
  lttf->SetTextFont(72);
  lttf->SetTextColor(kBlack);
  lttf->DrawLatex(1,1.35*ys[0],textthf);

  outstring.close();

}

//
// Compare Signal and Noise, same Run Range
//

void CompareSignalAndNoise(TString name,bool logscale, bool swtmsg){

  gStyle->SetOptStat(0);

  // Data
  TFile *l1 = TFile::Open("histo_Castor_threshold_p1_collisions.root");
  TFile *l2 = TFile::Open("histo_Castor_threshold_p1_unpaired.root");
  TFile *l3 = TFile::Open("histo_Castor_threshold_p1_nocollisions.root");

  TCanvas* c1 = new TCanvas(name,name);
  c1->cd()->SetRightMargin(0.04);

  TLegend* leg = new TLegend(0.696, 0.35, 0.94, 0.92);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);

  TString name1 = name + "_collisions";
  TString name2 = name + "_unpaired";
  TString name3 = name + "_no_collisions";

  if(logscale) c1->SetLogy(1);

  TH1F* h_1 = (TH1F*)l1->Get(name1);
  h_1->SetMarkerStyle(20);
  h_1->SetMarkerColor(kBlack);
  h_1->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_1,"Data 2010, Collisions","P");

  TH1F* h_2 = (TH1F*)l2->Get(name2);
  h_2->SetMarkerStyle(20);
  h_2->SetLineColor(kRed);
  h_2->SetLineWidth(1.6);
  h_2->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_2,"Data 2010, Unpaired","L");

  TH1F* h_3 = (TH1F*)l3->Get(name3);
  h_3->SetMarkerStyle(20);
  h_3->SetLineColor(kBlue);
  h_3->SetLineWidth(1.6);
  h_3->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_3,"Data 2010, No Collisions","L");

  Double_t xmin = h_2->GetXaxis()->GetXmin();
  Double_t xmax = h_2->GetXaxis()->GetXmax();

  cout << "X min: " << xmin << " | X max: " << xmax << endl;

  TAxis *axis = h_2->GetXaxis();
  Int_t bin_X_min = axis->FindBin(xmin);
  Int_t bin_X_max = axis->FindBin(xmax);

  cout << "First Bin: " << bin_X_min << " | Last Bin: " << bin_X_max << endl;

  Double_t binCenter_x_min = h_2->GetBinCenter(bin_X_min);
  Double_t binCenter_x_max = h_2->GetBinCenter(bin_X_max);

  int entries=0;
  for (int i = 0; i<= bin_X_max; i++){
    entries += h_2->GetBinContent(i);
  }

  int porcentage = 0.99*entries;
  int count =0;

  cout << "Total Entries: " << entries << endl;
  cout << "Porcentage: " << porcentage << endl;

  Double_t value_limit_95;

  for (int i=0; i<= bin_X_max; i++){

    count += h_2->GetBinContent(i);
    cout << "Count[" << i << "]: " << count << endl;

    if (count > porcentage){
      cout << "Passed Count[" << i << "]: " << count << endl;
      value_limit_95 = axis->GetBinCenter(i);
      break;
    }
  }

  Double_t maxvalue = h_2->GetMaximum();
  c1->Update();

  TLine *line = new TLine(value_limit_95,0,value_limit_95,maxvalue);
  line->SetLineColor(8);
  line->SetLineWidth(3);

  h_1->Scale(100./h_1->GetEntries());
  h_2->Scale(100./h_2->GetEntries());
  h_3->Scale(100./h_3->GetEntries());

  h_2->Draw("histo");
  h_1->Draw("psame");
  h_3->Draw("histosame");
  leg->Draw("same");
  gPad->RedrawAxis();

  if (swtmsg){
    text1 = new TLatex(5.0,10.0,"Castor Pedestals");
    text1->SetNDC();
    text1->SetTextAlign(13);
    text1->SetX(0.5);
    text1->SetY(0.5);
    text1->SetTextFont(42);
    text1->SetTextSizePixels(24);// dflt=28
    text1->Draw();
  }

  c1->SaveAs(name+TString(".png"));

}

//
// Draw Castor Noise with Different Run Ranges
//

void NoiseDifferentPeriods(TString name,bool logscale){

  gStyle->SetOptStat(0);

  // Data
  TFile *l1 = TFile::Open("histo_Castor_threshold_p1_unpaired.root");
  TFile *l2 = TFile::Open("histo_Castor_threshold_p2_unpaired.root");
  TFile *l3 = TFile::Open("histo_Castor_threshold_p3_unpaired.root");

  TCanvas* c1 = new TCanvas(name,name);
  c1->cd()->SetRightMargin(0.04);

  TLegend* leg = new TLegend(0.696, 0.35, 0.94, 0.92);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);

  TString name1 = name + "_unpaired";
  TString name2 = name + "_unpaired";
  TString name3 = name + "_unpaired";

  if(logscale) c1->SetLogy(1);

  TH1F* h_1 = (TH1F*)l1->Get(name1);
  h_1->SetMarkerStyle(20);
  h_1->SetMarkerColor(kBlack);
  h_1->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_1,"Run 136035-138559","L");

  TH1F* h_2 = (TH1F*)l2->Get(name2);
  h_2->SetMarkerStyle(20);
  h_2->SetLineColor(kRed);
  h_2->SetLineWidth(1.6);
  h_2->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_2,"Run 138560-141955","L");

  TH1F* h_3 = (TH1F*)l3->Get(name3);
  h_3->SetMarkerStyle(20);
  h_3->SetLineColor(kBlue);
  h_3->SetLineWidth(1.6);
  h_3->GetYaxis()->SetTitle("Fraction of Events %");
  leg->AddEntry(h_3,"Run 141956-149291","L");

  Double_t xmin = h_2->GetXaxis()->GetXmin();
  Double_t xmax = h_2->GetXaxis()->GetXmax();

  cout << "X min: " << xmin << " | X max: " << xmax << endl;

  TAxis *axis = h_2->GetXaxis();
  Int_t bin_X_min = axis->FindBin(xmin);
  Int_t bin_X_max = axis->FindBin(xmax);

  cout << "First Bin: " << bin_X_min << " | Last Bin: " << bin_X_max << endl;

  Double_t binCenter_x_min = h_2->GetBinCenter(bin_X_min);
  Double_t binCenter_x_max = h_2->GetBinCenter(bin_X_max);

  int entries=0;
  for (int i = 0; i<= bin_X_max; i++){
    entries += h_2->GetBinContent(i);
  }

  int porcentage = 0.99*entries;
  int count =0;

  cout << "Total Entries: " << entries << endl;
  cout << "Porcentage: " << porcentage << endl;

  Double_t value_limit_95;

  for (int i=0; i<= bin_X_max; i++){

    count += h_2->GetBinContent(i);
    cout << "Count[" << i << "]: " << count << endl;

    if (count > porcentage){
      cout << "Passed Count[" << i << "]: " << count << endl;
      value_limit_95 = axis->GetBinCenter(i);
      break;
    }
  }

  Double_t maxvalue = h_2->GetMaximum();
  c1->Update();

  TLine *line = new TLine(value_limit_95,0,value_limit_95,maxvalue);
  line->SetLineColor(8);
  line->SetLineWidth(3);

  h_1->Scale(100./h_1->GetEntries());
  h_2->Scale(100./h_2->GetEntries());
  h_3->Scale(100./h_3->GetEntries());

  h_2->Draw("histo");
  h_1->Draw("histosame");
  h_3->Draw("histosame");
  leg->Draw("same");
  gPad->RedrawAxis();
  line->Draw();

  c1->SaveAs(name+TString(".png"));

}

//
// Make Correction Evolution Plots
//
void CalculateCorrectorFlow(){

  TCanvas *c1 = new TCanvas("module1","module1",1200,1200);
  c1->Divide(4,4);

  TCanvas *c2 = new TCanvas("module2","module2",1200,1200);
  c2->Divide(4,4);

  TCanvas *c3 = new TCanvas("module3","module3",1200,1200);
  c3->Divide(4,4);

  TCanvas *c4 = new TCanvas("module4","module4",1200,1200);
  c4->Divide(4,4);

  TCanvas *c5 = new TCanvas("module5","module5",1200,1200);
  c5->Divide(4,4);

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile *l1 = TFile::Open("ratio_step1.root");
  TFile *l2 = TFile::Open("ratio_step2.root");
  TFile *l3 = TFile::Open("ratio_step3.root");
  TFile *l4 = TFile::Open("ratio_step4.root");
  TFile *l5 = TFile::Open("ratio_step5.root");
  TFile *l6 = TFile::Open("ratio_step6.root");

  TH2F* h_1 = (TH2F*)l1->Get("channelcorrector");
  TH2F* h_2 = (TH2F*)l2->Get("channelcorrector");
  TH2F* h_3 = (TH2F*)l3->Get("channelcorrector");
  TH2F* h_4 = (TH2F*)l4->Get("channelcorrector");
  TH2F* h_5 = (TH2F*)l5->Get("channelcorrector");
  TH2F* h_6 = (TH2F*)l6->Get("channelcorrector");

  for(int j=1;j<6;j++){
    for (int i=1;i<17;i++){

      char name[300];
      char title[300];
      char textfit[300];
      double acorr = 0;

      if(j==1){
	sprintf(name,"Channel%d",i);
	sprintf(title,"Channel %d;#sum E_{CASTOR}; Coefficient Factor",i);
	c1->cd(i);
      }

      if(j==2){
	sprintf(name,"Channel%d",i+16);
	sprintf(title,"Channel %d;#sum E_{CASTOR}; Coefficient Factor",i+16);
	c2->cd(i);
      }

      if(j==3){
	sprintf(name,"Channel%d",i+32);
	sprintf(title,"Channel %d;#sum E_{CASTOR}; Coefficient Factor",i+32);
	c3->cd(i);
      }

      if(j==4){
	sprintf(name,"Channel%d",i+48);
	sprintf(title,"Channel %d;#sum E_{CASTOR}; Coefficient Factor",i+48);
	c4->cd(i);
      }

      if(j==5){
	sprintf(name,"Channel%d",i+64);
	sprintf(title,"Channel %d;#sum E_{CASTOR}; Coefficient Factor",i+64);
	c5->cd(i);
      }

      TH1F* h=new TH1F(name,title, 6, 0, 300);
      h->SetBinContent(1,h_1->GetBinContent(j,i));
      h->SetBinContent(2,h_2->GetBinContent(j,i));
      h->SetBinContent(3,h_3->GetBinContent(j,i));
      h->SetBinContent(4,h_4->GetBinContent(j,i));
      h->SetBinContent(5,h_5->GetBinContent(j,i));
      h->SetBinContent(6,h_6->GetBinContent(j,i));
      acorr = (h_1->GetBinContent(j,i) + h_2->GetBinContent(j,i) + h_3->GetBinContent(j,i) + h_4->GetBinContent(j,i) + h_5->GetBinContent(j,i) + h_6->GetBinContent(j,i))/6.;
      /* h->SetBinError(1,h_1->GetBinError(j,i));
	 h->SetBinError(2,h_2->GetBinError(j,i));
	 h->SetBinError(3,h_3->GetBinError(j,i));
	 h->SetBinError(4,h_4->GetBinError(j,i));
	 h->SetBinError(5,h_5->GetBinError(j,i));
	 h->SetBinError(6,h_6->GetBinError(j,i));*/
      h->GetYaxis()->SetRangeUser(0.,3.);
      h->SetMarkerSize(1);
      h->SetMarkerStyle(20);
      h->Draw("P");
      //h->Fit("pol0");

      sprintf(textfit,"%g",acorr);
      TLatex *hfit = new TLatex(0.7,160,textfit);
      hfit->SetTextSize(0.1);
      hfit->SetTextFont(72);
      hfit->SetTextColor(kBlack);
      hfit->DrawLatex(175,2.5,textfit);

      TLine *line = new TLine(0,acorr,300,acorr);
      line->SetLineColor(8);
      line->SetLineWidth(3);
      line->Draw();

    }
  }

}


//
// Histograms: create correction mapping
//
void CalculateCorrector(TString name,TString output){

  TCanvas *c1 = new TCanvas("factors","factors",1400,400);
  c1->Divide(3,1);

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile *l1 = TFile::Open("histo_Castor_muon.root");
  TFile *l2 = TFile::Open("histo_DyToMuMu_RECO_nocor.root");

  TH2F* h_1 = (TH2F*)l1->Get(name);
  TH2F* h_2 = (TH2F*)l2->Get(name);

  int nBinsX = h_1->GetXaxis()->GetNbins();
  int nBinsY = h_2->GetYaxis()->GetNbins();

  // Integral Bin by Bin
  double factor1d = h_1->Integral(1,1,0,nBinsY);
  double factor2d = h_1->Integral(2,2,0,nBinsY);
  double factor3d = h_1->Integral(3,3,0,nBinsY);
  double factor4d = h_1->Integral(4,4,0,nBinsY);
  double factor5d = h_1->Integral(5,5,0,nBinsY);

  double factor1m = h_2->Integral(1,1,0,nBinsY);
  double factor2m = h_2->Integral(2,2,0,nBinsY);
  double factor3m = h_2->Integral(3,3,0,nBinsY);
  double factor4m = h_2->Integral(4,4,0,nBinsY);
  double factor5m = h_2->Integral(5,5,0,nBinsY);

  for(int i=1;i<=nBinsY; i++){
    h_1->SetBinContent(1,i,h_1->GetBinContent(1,i)/factor1d);
    h_1->SetBinContent(2,i,h_1->GetBinContent(2,i)/factor2d);
    h_1->SetBinContent(3,i,h_1->GetBinContent(3,i)/factor3d);
    h_1->SetBinContent(4,i,h_1->GetBinContent(4,i)/factor4d);
    h_1->SetBinContent(5,i,h_1->GetBinContent(5,i)/factor5d);

    h_2->SetBinContent(1,i,h_2->GetBinContent(1,i)/factor1m);
    h_2->SetBinContent(2,i,h_2->GetBinContent(2,i)/factor2m);
    h_2->SetBinContent(3,i,h_2->GetBinContent(3,i)/factor3m);
    h_2->SetBinContent(4,i,h_2->GetBinContent(4,i)/factor4m);
    h_2->SetBinContent(5,i,h_2->GetBinContent(5,i)/factor5m);
  }

  double factor1ds = h_1->Integral(1,1,0,nBinsY);
  double factor2ds = h_1->Integral(2,2,0,nBinsY);
  double factor3ds = h_1->Integral(3,3,0,nBinsY);
  double factor4ds = h_1->Integral(4,4,0,nBinsY);
  double factor5ds = h_1->Integral(5,5,0,nBinsY);

  double factor1ms = h_2->Integral(1,1,0,nBinsY);
  double factor2ms = h_2->Integral(2,2,0,nBinsY);
  double factor3ms = h_2->Integral(3,3,0,nBinsY);
  double factor4ms = h_2->Integral(4,4,0,nBinsY);
  double factor5ms = h_2->Integral(5,5,0,nBinsY);

  TH2F* ratio = (TH2F*)h_1->Clone();
  ratio->SetTitle("Correction Channels Factors");
  ratio->Divide(h_2);

  c1->cd(1);
  h_1->SetTitle("Pattern Data Energy, normalized");
  h_1->Draw();
  c1->cd(2);
  h_2->SetTitle("Pattern MC Energy, normalized");
  h_2->Draw();
  c1->cd(3);
  ratio->Draw();

  TFile f(output,"recreate");
  ratio->Write("channelcorrector");

}


//
// Histograms: configure colors
//

void TH1Draw(TH1F* h_1, TH1F* h_2, TH1F* h_3){

  h_1->SetLineColor(kBlack);
  h_1->SetMarkerStyle(20);
  h_1->SetMarkerSize(0.8);
  h_1->GetYaxis()->SetTitle("N Events");

  h_2->SetLineColor(kRed);
  h_2->SetLineWidth(2);
  h_2->SetFillColor(kRed-4);
  h_2->SetFillStyle(3020);
  h_2->GetYaxis()->SetTitle("N Events");

  h_3->SetLineColor(kBlue);
  h_3->SetLineWidth(2);
  h_3->SetFillColor(kBlue-4);
  h_3->SetFillStyle(3020);
  h_3->GetYaxis()->SetTitle("N Events");

}
