
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include <math.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <string>
#include <iostream>
#include <cmath>
#include "TLegend.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TPaveText.h"
#include "tdrstyle.C"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"

void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.2);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.1);
  //pad1->SetGrid();
  pad1->SetGrid(10,10);
}

void plotEfficiency(){
  TString fileName = "MiniAOD_effi_80x_DYtoLL.root"; 
  TString treePath = "cutBased/Ntuple";
  int bins = 20;
  int low  = 0;
  int high = 100;

  //Plotting Variables
  TString variable = "tauPt";
  TString GenCut= "genTauPt>20&&abs(genTauEta)<2.4";
  TString RecoCut= "tauPt>18&&dmf>-1&&goodReco_>0&&" + GenCut;

  //Style
  TString xaxis = "Tau Pt [GeV]";
  int markerstyle = 20;
  Color_t color = kRed;
  TString outFileName = "plotEfficiency";

  setTDRStyle();
  //tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");

  TFile *tauFile    = new TFile("MiniAOD_effi_80x_DYtoLL.root");

  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600); Tcan->cd();  Tcan->SetFillColor(0);
  TPad* pad1 = new TPad("pad1","The pad",0,0,0.95,0.95);

  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  TTree* tauTree = (TTree*)tauFile->Get(treePath);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  TH1F* Denom;
  Denom = new TH1F("Denom","Denom",bins,low,high);
  Denom->Sumw2();
  tauTree->Draw(variable+">>+Denom",GenCut);

  TH1F* Num;
  Num = new TH1F("Num","Num",bins,low,high);
  tauTree->Draw(variable+">>+Num",RecoCut);

  Num->Divide(Denom);
  gStyle->SetErrorX(0.5);

  Num->GetXaxis()->SetTitle(xaxis);
  Num->GetYaxis()->SetTitle("Efficiency ");
  Num->GetYaxis()->SetTitleOffset(0.8);
  Num->SetMarkerStyle(markerstyle);
  Num->SetMarkerColor(color+2);
  Num->Draw("P");

  Num->SetFillStyle(1001);
  Num->SetLineWidth((short)1.5);
  Num->SetMaximum(1.1);
  Num->SetMinimum(0);

  Tcan->cd();
  Tcan->SaveAs(outFileName+".pdf");
}
