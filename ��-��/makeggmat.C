//Make prompt g-g coincidence matrix using RADWARE approach.
//Usage: root -l makeggmat.C
//The output file will be analyzed by gg.C
#include <iostream>
#include <sstream>
#include "TString.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"

TString rootfile="../decay46_123_all.root"; // Name of input file
TString ggfile="46_123gg600s.root"; // Name of output file

void makeggmat()
{
  int nebin=2500;//number of bin for TH2I
  double emax=2500;//maximum energy .
  double dtmin=300;//min. of decaytime
  double dtmax=800;//max. of decaytime.
  double ctmax=400;//coincidence time width
  bool   backsub=false;//true;
  cout<<"Generate gg matix to ["<<ggfile<<"] from ["<<rootfile<<"]"<<endl;
  TFile *fin=TFile::Open(rootfile);
  TTree *tree=(TTree*)fin->Get("tree");
  //make gg 2D histogram
  //----------------------------------------------------------
  // TString atcut="caxt>0&&caxt<1000 && cayt>0&&cayt<1000 && abs(caxt-cayt)<200";
  TString sctcut=Form("abs(caxt-cayt)<%f",ctmax);
  TString sgtree=Form("caxe:caye>>gg(%i,0,%f,%i,0,%f)",nebin,emax,nebin,emax);
  TString sgcut= Form("%s && decaytime>%f&&decaytime<%f",sctcut.Data(),dtmin,dtmax);
  tree->Draw(sgtree,sgcut);
  auto hgg=(TH2I*)gROOT->FindObject("gg"); //gg matrix
  TString sbtree=Form("caxe:caye>>bg0(%i,0,%f,%i,0,%f)",nebin,emax,nebin,emax);
  TString sbcut= Form("%s && decaytime>%f && decaytime<%f",sctcut.Data(),4000.,4000+dtmax-dtmin);  
  tree->Draw(sbtree,sbcut);  
  auto hbg0=(TH2I*)gROOT->FindObject("bg0"); //bgg matrix
  if(backsub) hgg->Add(hgg,hbg0,1,-1);//background subtruction
  //make ae histogram.
  //----------------------------------------------------------
  //TString satcut="at>0 && at<1000";
  TString satcut="1";
  TString sgtree1=Form("ae>>ae(%i,0,%f)",nebin*2,emax);//for precise peak value;
  TString sgcut1= Form("%s && decaytime>%f && decaytime<%f",satcut.Data(),dtmin,dtmax);
  TString sbtree1=Form("ae>>aeb0(%i,0,%f)",nebin,emax);
  TString sbcut1= Form("%s && decaytime>%f && decaytime<%f",satcut.Data(),4000.,4000+dtmax-dtmin);  
  tree->Draw(sgtree1,sgcut1);
  auto hae=(TH1D*)gROOT->FindObject("ae"); //gg matrix
  tree->Draw(sbtree1,sbcut);  
  auto haeb0=(TH1D*)gROOT->FindObject("aeb0"); //bgg matrix
  if(backsub) hae->Add(hae,haeb0,1,-1);//background subtruction.

  //----------------------------------------------------------
  auto hx=(TH1D*)hgg->ProjectionX("Tpj"); //total projection spectrum(Tpj)
  TSpectrum *sx= new TSpectrum(500);
  Int_t nfoundx=sx->Search(hx,20,"",0.1);
  auto hbx=sx->Background(hx,9,"same"); //background of Tpj
  hbx->SetName("TpjBg");hbx->SetTitle("TpjBg");
  auto hpeakx=(TH2I*)hbx->Clone("TpjPeak");
  hpeakx->Add(hx,hbx,1,-1); //net peaks of Tpj
  hpeakx->Draw("same");

  TH2I* hggb=new TH2I("ggbmat","bgmat for gg",nebin,0,emax,nebin,0,emax);
  hggb->Reset();
  Double_t T,Pi,Pj,pi,pj,Bij,x,y;
  T=hx->Integral();
  for(int i=0;i<hgg->GetNbinsX();i++) {
    for(int j=0;j<hgg->GetNbinsY();j++) {
      Pi=hx->GetBinContent(i+1);
      Pj=hx->GetBinContent(j+1);
      pi=hpeakx->GetBinContent(i+1);
      pj=hpeakx->GetBinContent(j+1);
      Bij=(Pi*Pj-pi*pj)/T;
      x=hx->GetBinCenter(i+1);
      y=hx->GetBinCenter(j+1);
      hggb->Fill(x,y,Bij);// backgound gg matix
    }
  }
  TH2I* hggmat=new TH2I("ggmat","gg with backsub",nebin,0,emax,nebin,0,emax);
  hggmat->Add(hgg,hggb,1,-1);// background subtracted gg matix.
  TFile *fout=new TFile(ggfile,"RECREATE");
  hgg->Write();//"gg"- gg matix
  hae->Write();
  hx->Write();//"Tpj" -total projection spectrum
  hbx->Write();//"TpjBg" - background of Tpj
  hpeakx->Write();//"TpjPeak" -net peaks of Tpj
  hggb->Write();//"ggbmat" - background gg matix
  hggmat->Write();//"ggmat" -background subtracted gg matix 
  fout->Write();
  fout->Close(); 
 
}

