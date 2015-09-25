/*

   >> P R O G R A M  F I T T E R
   -----------------------------
GOAL: Create fit object to be used in Diffractive Z/W analysis,
in order to correct castor energy distribution

PROJECT: Diffractive Z/W boson analysis

 */

#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TAxis.h"
#include <vector>

std::string mc_file_bkg;
std::string mc_file_signal;

double xmin, xmax;
int rebin;

void GenerateCastorFitFile(){

  data = "";
  mc_file_bkg = "histo_castor_DyToEE_Reco_pu.root";
  mc_file_signal = "histo_castor_PomwigMinus_electron_Reco.root";

  // X axis range. If xmin == xmax, no changes.
  xmin = 0.;
  xmax = 1500.;

  // Rebin. If rebin <= 0, no changes.
  rebin = 20;

  MyResolution("sumECastorMinus_zeropileup","gensumECastorMinus_zeropileup");
  //MyResolution("sumECastorMinus_zeropileup","gensumECastorPionsMinus_zeropileup");

}

// C M S   O F F I C I A L   S T Y L E
//------------------------------------

void loadFiles(TString hname1, TString hname2){

  TFile *l_file_signal  = TFile::Open(mc_file_signal.c_str());
  TFile *l_file_bkg  = TFile::Open(mc_file_bkg.c_str());

  TH1F* h_reco_signal = (TH1F*)l_file_signal->Get(hname1);
  TH1F* h_gen_signal = (TH1F*)l_file_signal->Get(hname2);

  TH1F* h_reco_bkg = (TH1F*)l_file_bkg->Get(hname1);
  TH1F* h_gen_bkg = (TH1F*)l_file_bkg->Get(hname2);

  //h_reco_signal->Scale(0.3*0.001329954557739);
  //h_reco_bkg->Scale(0.020137974016854);

  //h_gen_signal->Scale(0.3*0.001329954557739);
  //h_gen_bkg->Scale(0.020137974016854);

  TList *listreco = new TList;
  listreco->Add(h_reco_bkg);
  listreco->Add(h_reco_signal);

  TList *listgen = new TList;
  listgen->Add(h_gen_bkg);
  listgen->Add(h_gen_signal);

  TH1F *reco = (TH1F*)h_reco_bkg->Clone("reco");
  reco->Reset();
  reco->Merge(listreco);

  TH1F *gen = (TH1F*)h_gen_bkg->Clone("gen");
  gen->Reset();
  gen->Merge(listgen);

}

void MyResolution(TString hname1, TString hname2){

  loadFiles(hname1, hname2);

  if(xmin != xmax) reco->GetXaxis()->SetRangeUser(xmin,xmax);
  if(rebin>0) reco->Rebin(rebin);
  reco->Sumw2();

  if(xmin != xmax) gen->GetXaxis()->SetRangeUser(xmin,xmax);
  gen->Sumw2();
  if(rebin>0) gen->Rebin(rebin);

  TH1F *ratio=(TH1F*)reco->Clone();
  ratio->Divide(gen);
  ratio->SetStats(0);

  g1 = new TF1("linearity","pol9",0.,1500.);
  ratio->Fit(g1,"R");

  TFile *f1 = new TFile("ToFit.root","recreate"); 
  f1->cd();
  g1->SetNpx(500);
  g1->Write();
  ratio->Write();
  f1->Close();

}
