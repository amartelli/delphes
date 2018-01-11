#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"





void drawHistos(){

  gROOT->Reset();     
  gROOT->Macro("~/public/setStyle.C");

  TFile *_file0;
  _file0 =  TFile::Open("outputAnalysis/results.root");  


  int iColors[4] = {kBlack, kBlue, kRed, kGreen+2};

  std::vector<std::string> histoNames;
  histoNames.push_back("htrkPtRel_bFlav");
  histoNames.push_back("htrkPtRel_cFlav");
  histoNames.push_back("htrkPtRel_lightFlav");

  TH1F* htrkPtRel[3];
  for(int ij=0; ij<3; ++ij){
    htrkPtRel[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkPtRel[ij]->SetName(histoNames.at(ij).c_str());
    htrkPtRel[ij]->SetLineColor(iColors[ij+1]);
    htrkPtRel[ij]->SetMarkerColor(iColors[ij+1]);
    htrkPtRel[ij]->SetMarkerStyle(7);
    htrkPtRel[ij]->Rebin(1);
  }

  histoNames.clear();
  histoNames.push_back("htrkIP2Dxy_bFlav");
  histoNames.push_back("htrkIP2Dxy_cFlav");
  histoNames.push_back("htrkIP2Dxy_lightFlav");

  TH1F* htrkD0[3];
  for(int ij=0; ij<3; ++ij){
    htrkD0[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkD0[ij]->SetName(histoNames.at(ij).c_str());
    htrkD0[ij]->SetLineColor(iColors[ij+1]);
    htrkD0[ij]->SetMarkerColor(iColors[ij+1]);
    htrkD0[ij]->SetMarkerStyle(7);
    htrkD0[ij]->Rebin(4);
  }

  histoNames.clear();
  histoNames.push_back("htrkEtaRel_bFlav");
  histoNames.push_back("htrkEtaRel_cFlav");
  histoNames.push_back("htrkEtaRel_lightFlav");

  TH1F* htrkEtaRel[3];
  for(int ij=0; ij<3; ++ij){
    htrkEtaRel[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkEtaRel[ij]->SetName(histoNames.at(ij).c_str());
    htrkEtaRel[ij]->SetLineColor(iColors[ij+1]);
    htrkEtaRel[ij]->SetMarkerColor(iColors[ij+1]);
    htrkEtaRel[ij]->SetMarkerStyle(7);
    htrkEtaRel[ij]->Rebin(4);
  }


  histoNames.clear();
  histoNames.push_back("htrkPPar_bFlav");
  histoNames.push_back("htrkPPar_cFlav");
  histoNames.push_back("htrkPPar_lightFlav");

  TH1F* htrkPPar[3];
  for(int ij=0; ij<3; ++ij){
    htrkPPar[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkPPar[ij]->SetName(histoNames.at(ij).c_str());
    htrkPPar[ij]->SetLineColor(iColors[ij+1]);
    htrkPPar[ij]->SetMarkerColor(iColors[ij+1]);
    htrkPPar[ij]->SetMarkerStyle(7);
    htrkPPar[ij]->Rebin(2);
  }

  histoNames.clear();
  histoNames.push_back("htrkPPar_bFlav");
  histoNames.push_back("htrkPPar_cFlav");
  histoNames.push_back("htrkPPar_lightFlav");

  TH1F* hGenPtrkPPar[3];
  for(int ij=0; ij<3; ++ij){
    hGenPtrkPPar[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    hGenPtrkPPar[ij]->SetName(histoNames.at(ij).c_str());
    hGenPtrkPPar[ij]->SetLineColor(iColors[ij+1]);
    hGenPtrkPPar[ij]->SetMarkerColor(iColors[ij+1]);
    hGenPtrkPPar[ij]->SetMarkerStyle(7);
    hGenPtrkPPar[ij]->Rebin(2);
  }


  //pT > 80

  histoNames.clear();
  histoNames.push_back("htrkPtRel_bFlav_80");
  histoNames.push_back("htrkPtRel_cFlav_80");
  histoNames.push_back("htrkPtRel_lightFlav_80");

  TH1F* htrkPtRel_80[3];
  for(int ij=0; ij<3; ++ij){
    htrkPtRel_80[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkPtRel_80[ij]->SetName(histoNames.at(ij).c_str());
    htrkPtRel_80[ij]->SetLineColor(iColors[ij+1]);
    htrkPtRel_80[ij]->SetMarkerColor(iColors[ij+1]);
    htrkPtRel_80[ij]->SetMarkerStyle(7);
    htrkPtRel[ij]->Rebin(2);
  }

  histoNames.clear();
  histoNames.push_back("htrkIP2Dxy_bFlav_80");
  histoNames.push_back("htrkIP2Dxy_cFlav_80");
  histoNames.push_back("htrkIP2Dxy_lightFlav_80");

  TH1F* htrkD0_80[3];
  for(int ij=0; ij<3; ++ij){
    htrkD0_80[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkD0_80[ij]->SetName(histoNames.at(ij).c_str());
    htrkD0_80[ij]->SetLineColor(iColors[ij+1]);
    htrkD0_80[ij]->SetMarkerColor(iColors[ij+1]);
    htrkD0_80[ij]->SetMarkerStyle(7);
    htrkD0_80[ij]->Rebin(4);
  }

  // gentrack
  histoNames.clear();
  histoNames.push_back("hGenPtrkPtRel_bFlav");
  histoNames.push_back("hGenPtrkPtRel_cFlav");
  histoNames.push_back("hGenPtrkPtRel_lightFlav");

  TH1F* hGenPtrkPtRel[3];
  for(int ij=0; ij<3; ++ij){
    hGenPtrkPtRel[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    hGenPtrkPtRel[ij]->SetName(histoNames.at(ij).c_str());
    hGenPtrkPtRel[ij]->SetLineColor(iColors[ij+1]);
    hGenPtrkPtRel[ij]->SetMarkerColor(iColors[ij+1]);
    hGenPtrkPtRel[ij]->SetMarkerStyle(7);
    hGenPtrkPtRel[ij]->Rebin(2);
  }

  histoNames.clear();
  histoNames.push_back("hGenPtrkIP2Dxy_bFlav");
  histoNames.push_back("hGenPtrkIP2Dxy_cFlav");
  histoNames.push_back("hGenPtrkIP2Dxy_lightFlav");

  TH1F* hGenPtrkD0[3];
  for(int ij=0; ij<3; ++ij){
    hGenPtrkD0[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    hGenPtrkD0[ij]->SetName(histoNames.at(ij).c_str());
    hGenPtrkD0[ij]->SetLineColor(iColors[ij+1]);
    hGenPtrkD0[ij]->SetMarkerColor(iColors[ij+1]);
    hGenPtrkD0[ij]->SetMarkerStyle(7);
    hGenPtrkD0[ij]->Rebin(4);
  }

  //computed 
  histoNames.clear();
  histoNames.push_back("htrkIP2Dxy_bFlav_m");
  histoNames.push_back("htrkIP2Dxy_cFlav_m");
  histoNames.push_back("htrkIP2Dxy_lightFlav_m");

  TH1F* htrkSip2D[3];
  for(int ij=0; ij<3; ++ij){
    htrkSip2D[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkSip2D[ij]->SetName(histoNames.at(ij).c_str());
    htrkSip2D[ij]->SetLineColor(iColors[ij+1]);
    htrkSip2D[ij]->SetMarkerColor(iColors[ij+1]);
    htrkSip2D[ij]->SetMarkerStyle(7);
    htrkSip2D[ij]->Rebin(4);
  }

  histoNames.clear();
  histoNames.push_back("htrkIP3D_bFlav_m");
  histoNames.push_back("htrkIP3D_cFlav_m");
  histoNames.push_back("htrkIP3D_lightFlav_m");

  TH1F* htrkSip3D[3];
  for(int ij=0; ij<3; ++ij){
    htrkSip3D[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkSip3D[ij]->SetName(histoNames.at(ij).c_str());
    htrkSip3D[ij]->SetLineColor(iColors[ij+1]);
    htrkSip3D[ij]->SetMarkerColor(iColors[ij+1]);
    htrkSip3D[ij]->SetMarkerStyle(7);
    htrkSip3D[ij]->Rebin(4);
  }

  //htrackSumJetEtRatio
  histoNames.clear();
  histoNames.push_back("htrackSumJetEtRatio_bFlav");
  histoNames.push_back("htrackSumJetEtRatio_cFlav");
  histoNames.push_back("htrackSumJetEtRatio_lightFlav");

  TH1F* htrkSumJetEtRatio[3];
  for(int ij=0; ij<3; ++ij){
    htrkSumJetEtRatio[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkSumJetEtRatio[ij]->SetName(histoNames.at(ij).c_str());
    htrkSumJetEtRatio[ij]->SetLineColor(iColors[ij+1]);
    htrkSumJetEtRatio[ij]->SetMarkerColor(iColors[ij+1]);
    htrkSumJetEtRatio[ij]->SetMarkerStyle(7);
    htrkSumJetEtRatio[ij]->Rebin(4);
  }

  //htrackSumJetDeltaR
  histoNames.clear();
  histoNames.push_back("htrackSumJetDeltaR_bFlav");
  histoNames.push_back("htrackSumJetDeltaR_cFlav");
  histoNames.push_back("htrackSumJetDeltaR_lightFlav");

  TH1F* htrkSumJetDeltaR[3];
  for(int ij=0; ij<3; ++ij){
    htrkSumJetDeltaR[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    htrkSumJetDeltaR[ij]->SetName(histoNames.at(ij).c_str());
    htrkSumJetDeltaR[ij]->SetLineColor(iColors[ij+1]);
    htrkSumJetDeltaR[ij]->SetMarkerColor(iColors[ij+1]);
    htrkSumJetDeltaR[ij]->SetMarkerStyle(7);
    htrkSumJetDeltaR[ij]->Rebin(4);
  }

  //gen
  histoNames.clear();
  histoNames.push_back("hGenPtrackSumJetEtRatio_bFlav");
  histoNames.push_back("hGenPtrackSumJetEtRatio_cFlav");
  histoNames.push_back("hGenPtrackSumJetEtRatio_lightFlav");

  TH1F* hGenPtrkSumJetEtRatio[3];
  for(int ij=0; ij<3; ++ij){
    hGenPtrkSumJetEtRatio[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    hGenPtrkSumJetEtRatio[ij]->SetName(histoNames.at(ij).c_str());
    hGenPtrkSumJetEtRatio[ij]->SetLineColor(iColors[ij+1]);
    hGenPtrkSumJetEtRatio[ij]->SetMarkerColor(iColors[ij+1]);
    hGenPtrkSumJetEtRatio[ij]->SetMarkerStyle(7);
    hGenPtrkSumJetEtRatio[ij]->Rebin(4);
  }

  histoNames.clear();
  histoNames.push_back("hGenPtrackSumJetDeltaR_bFlav");
  histoNames.push_back("hGenPtrackSumJetDeltaR_cFlav");
  histoNames.push_back("hGenPtrackSumJetDeltaR_lightFlav");

  TH1F* hGenPtrkSumJetDeltaR[3];
  for(int ij=0; ij<3; ++ij){
    hGenPtrkSumJetDeltaR[ij] = (TH1F*)_file0->Get(histoNames.at(ij).c_str());
    hGenPtrkSumJetDeltaR[ij]->SetName(histoNames.at(ij).c_str());
    hGenPtrkSumJetDeltaR[ij]->SetLineColor(iColors[ij+1]);
    hGenPtrkSumJetDeltaR[ij]->SetMarkerColor(iColors[ij+1]);
    hGenPtrkSumJetDeltaR[ij]->SetMarkerStyle(7);
    hGenPtrkSumJetDeltaR[ij]->Rebin(4);
  }



  TH1F* hFlavor = (TH1F*)_file0->Get("hFlavor");
  TH1F* hJetPt = (TH1F*)_file0->Get("hJetPt");
  TH1F* hJetEta = (TH1F*)_file0->Get("hJetEta");
  TH1F* hNNeutrals = (TH1F*)_file0->Get("hNNeutrals");
  TH1F* hNCharged  = (TH1F*)_file0->Get("hNCharged");
  //  TH1F* htrkEtaRel = (TH1F*)_file0->Get("htrkEtaRel");
  // TH1F* htrkSumJetEtRatio = (TH1F*)_file0->Get("htrackSumJetEtRatio");
  // TH1F* htrkSumJetDeltaR = (TH1F*)_file0->Get("htrackSumJetDeltaR");


  TLegend *legAll = new TLegend(0.7,0.75,0.85,0.90,NULL,"brNDC");
  legAll->SetTextFont(42);
  legAll->SetTextSize(0.04);
  legAll->SetFillColor(kWhite);
  legAll->SetLineColor(kWhite);
  legAll->SetShadowColor(kWhite);
  legAll->AddEntry(htrkPtRel[0], "b jet", "l");
  legAll->AddEntry(htrkPtRel[1], "c jet", "l");
  legAll->AddEntry(htrkPtRel[2], "udsg ", "l");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* c1 = new TCanvas();
  c1->cd();
  htrkPtRel[0]->GetXaxis()->SetTitle("PtRel");
  htrkPtRel[0]->GetYaxis()->SetRangeUser(0., 0.05);
  htrkPtRel[0]->Draw("h");
  htrkPtRel[1]->Draw("h, same");
  htrkPtRel[2]->Draw("h, same");
  legAll->Draw("same");
  c1->Print("outputAnalysis/plotsCompare_noSmear/htrkPtRel.png", "png");
  c1->Print("outputAnalysis/plotsCompare_noSmear/htrkPtRel.pdf", "pdf");
  c1->Print("outputAnalysis/plotsCompare_noSmear/htrkPtRel.root", "root");


  TCanvas* c1g = new TCanvas();
  c1g->cd();
  hGenPtrkPtRel[0]->GetXaxis()->SetTitle("PtRel");
  hGenPtrkPtRel[0]->GetYaxis()->SetRangeUser(0., 0.05);
  hGenPtrkPtRel[0]->Draw("h");
  hGenPtrkPtRel[1]->Draw("h, same");
  hGenPtrkPtRel[2]->Draw("h, same");
  legAll->SetHeader("genTrack");
  legAll->Draw("same");
  c1g->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkPtRel.png", "png");
  c1g->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkPtRel.pdf", "pdf");
  c1g->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkPtRel.root", "root");


  TCanvas* c2 = new TCanvas();
  c2->cd();
  htrkPtRel_80[0]->GetXaxis()->SetTitle("PtRel");
  //  htrkPtRel[0]->GetYaxis()->SetTitle("events");
  htrkPtRel_80[0]->Draw("h");
  htrkPtRel_80[1]->Draw("h, same");
  htrkPtRel_80[2]->Draw("h, same");
  legAll->SetHeader("jet pT > 80GeV");
  legAll->Draw("same");
  c2->Print("outputAnalysis/plotsCompare_noSmear/htrkPtRel_80.png", "png");
  c2->Print("outputAnalysis/plotsCompare_noSmear/htrkPtRel_80.pdf", "pdf");
  c2->Print("outputAnalysis/plotsCompare_noSmear/htrkPtRel_80.root", "root");




  TCanvas* c3 = new TCanvas();
  gPad->SetLogy();
  c3->cd();
  htrkD0[0]->GetXaxis()->SetTitle("d0 (mm)");
  //  htrkD0[0]->GetYaxis()->SetTitle("events");
  htrkD0[0]->Draw("h");
  htrkD0[1]->Draw("h, same");
  htrkD0[2]->Draw("h, same");
  legAll->SetHeader("");
  legAll->Draw("same");
  c3->Print("outputAnalysis/plotsCompare_noSmear/htrkD0.png", "png");
  c3->Print("outputAnalysis/plotsCompare_noSmear/htrkD0.pdf", "pdf");
  c3->Print("outputAnalysis/plotsCompare_noSmear/htrkD0.root", "root");

  TCanvas* c3g = new TCanvas();
  gPad->SetLogy();
  c3g->cd();
  hGenPtrkD0[0]->GetXaxis()->SetTitle("d0 (mm)");
  //  htrkD0[0]->GetYaxis()->SetTitle("events");
  hGenPtrkD0[0]->Draw("h");
  hGenPtrkD0[1]->Draw("h, same");
  hGenPtrkD0[2]->Draw("h, same");
  legAll->SetHeader("genTrack");
  legAll->Draw("same");
  c3g->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkD0.png", "png");
  c3g->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkD0.pdf", "pdf");
  c3g->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkD0.root", "root");

  TCanvas* c2d = new TCanvas();
  gPad->SetLogy();
  c2d->cd();
  htrkSip2D[0]->GetXaxis()->SetTitle("sip 2D");
  //  htrkD0[0]->GetYaxis()->SetTitle("events");
  htrkSip2D[0]->Draw("h");
  htrkSip2D[1]->Draw("h, same");
  htrkSip2D[2]->Draw("h, same");
  legAll->SetHeader("computed");
  legAll->Draw("same");
  c2d->Print("outputAnalysis/plotsCompare_noSmear/htrkSip2D.png", "png");
  c2d->Print("outputAnalysis/plotsCompare_noSmear/htrkSip2D.pdf", "pdf");
  c2d->Print("outputAnalysis/plotsCompare_noSmear/htrkSip2D.root", "root");

  TCanvas* c3d = new TCanvas();
  gPad->SetLogy();
  c3d->cd();
  htrkSip3D[0]->GetXaxis()->SetTitle("sip 3D");
  //  htrkD0[0]->GetYaxis()->SetTitle("events");
  htrkSip3D[0]->Draw("h");
  htrkSip3D[1]->Draw("h, same");
  htrkSip3D[2]->Draw("h, same");
  legAll->SetHeader("computed");
  legAll->Draw("same");
  c3d->Print("outputAnalysis/plotsCompare_noSmear/htrkSip3D.png", "png");
  c3d->Print("outputAnalysis/plotsCompare_noSmear/htrkSip3D.pdf", "pdf");
  c3d->Print("outputAnalysis/plotsCompare_noSmear/htrkSip3D.root", "root");


  TCanvas* c4 = new TCanvas();
  gPad->SetLogy();
  c4->cd();
  htrkD0_80[0]->GetXaxis()->SetTitle("d0 (mm)");
  //  htrkD0[0]->GetYaxis()->SetTitle("events");
  htrkD0_80[0]->Draw("h");
  htrkD0_80[1]->Draw("h, same");
  htrkD0_80[2]->Draw("h, same");
  legAll->SetHeader("jet pT > 80GeV");
  legAll->Draw("same");
  c4->Print("outputAnalysis/plotsCompare_noSmear/htrkD0_80.png", "png");
  c4->Print("outputAnalysis/plotsCompare_noSmear/htrkD0_80.pdf", "pdf");
  c4->Print("outputAnalysis/plotsCompare_noSmear/htrkD0_80.root", "root");


  TCanvas* c5 = new TCanvas();
  c5->cd();
  hFlavor->GetXaxis()->SetTitle("jet pdgID");
  hFlavor->Draw();
  c5->Print("outputAnalysis/plotsCompare_noSmear/hFlavour.png", "png");
  c5->Print("outputAnalysis/plotsCompare_noSmear/hFlavour.pdf", "pdf");
  c5->Print("outputAnalysis/plotsCompare_noSmear/hFlavour.root", "root");


  TCanvas* c6 = new TCanvas();
  gPad->SetLogy();
  c6->cd();
  hJetPt->GetXaxis()->SetTitle("jet pT");
  hJetPt->Draw();
  c6->Print("outputAnalysis/plotsCompare_noSmear/hJetPt.png", "png");
  c6->Print("outputAnalysis/plotsCompare_noSmear/hJetPt.pdf", "pdf");
  c6->Print("outputAnalysis/plotsCompare_noSmear/hJetPt.root", "root");

  TCanvas* c7 = new TCanvas();
  c7->cd();
  hJetEta->GetXaxis()->SetTitle("jet #eta");
  hJetEta->Draw();
  c7->Print("outputAnalysis/plotsCompare_noSmear/hJetEta.png", "png");
  c7->Print("outputAnalysis/plotsCompare_noSmear/hJetEta.pdf", "pdf");
  c7->Print("outputAnalysis/plotsCompare_noSmear/hJetEta.root", "root");


  TCanvas* c8 = new TCanvas();
  c8->cd();
  hNNeutrals->GetXaxis()->SetTitle("#n. neutrals in jet");
  hNNeutrals->Draw();
  c8->Print("outputAnalysis/plotsCompare_noSmear/hNNeutrals.png", "png");
  c8->Print("outputAnalysis/plotsCompare_noSmear/hNNeutrals.pdf", "pdf");
  c8->Print("outputAnalysis/plotsCompare_noSmear/hNNeutrals.root", "root");

  TCanvas* c9 = new TCanvas();
  c9->cd();
  hNCharged->GetXaxis()->SetTitle("#n. charged in jet");
  hNCharged->Draw();
  c9->Print("outputAnalysis/plotsCompare_noSmear/hNCharged.png", "png");
  c9->Print("outputAnalysis/plotsCompare_noSmear/hNCharged.pdf", "pdf");
  c9->Print("outputAnalysis/plotsCompare_noSmear/hNCharged.root", "root");

  ////
  TCanvas* c10 = new TCanvas();
  c10->cd();
  htrkEtaRel[0]->GetXaxis()->SetTitle("EtaRel");
  htrkEtaRel[0]->Draw("h");
  htrkEtaRel[1]->Draw("h, same");
  htrkEtaRel[2]->Draw("h, same");
  legAll->SetHeader("");
  legAll->Draw("same");
  c10->Print("outputAnalysis/plotsCompare_noSmear/htrkEtaRel.png", "png");
  c10->Print("outputAnalysis/plotsCompare_noSmear/htrkEtaRel.pdf", "pdf");
  c10->Print("outputAnalysis/plotsCompare_noSmear/htrkEtaRel.root", "root");


  TCanvas* c11 = new TCanvas();
  c11->cd();
  htrkSumJetEtRatio[0]->GetXaxis()->SetTitle("#Sigma(trk)_{ET})/ Jet_{ET}");
  htrkSumJetEtRatio[0]->Draw("h");
  htrkSumJetEtRatio[1]->Draw("h, same");
  htrkSumJetEtRatio[2]->Draw("h, same");
  legAll->SetHeader("");
  legAll->Draw("same");
  c11->Print("outputAnalysis/plotsCompare_noSmear/htrkSumJetEtRatio.png", "png");
  c11->Print("outputAnalysis/plotsCompare_noSmear/htrkSumJetEtRatio.pdf", "pdf");
  c11->Print("outputAnalysis/plotsCompare_noSmear/htrkSumJetEtRatio.root", "root");

  TCanvas* c12 = new TCanvas();
  gPad->SetLogy();
  c12->cd();
  htrkSumJetDeltaR[0]->GetXaxis()->SetTitle("DeltaR(#Sigma(trk), Jet)");
  htrkSumJetDeltaR[0]->Draw("h");
  htrkSumJetDeltaR[1]->Draw("h, same");
  htrkSumJetDeltaR[2]->Draw("h, same");
  legAll->SetHeader("");
  legAll->Draw("same");
  c12->Print("outputAnalysis/plotsCompare_noSmear/htrkSumJetDeltaR.png", "png");
  c12->Print("outputAnalysis/plotsCompare_noSmear/htrkSumJetDeltaR.pdf", "pdf");
  c12->Print("outputAnalysis/plotsCompare_noSmear/htrkSumJetDeltaR.root", "root");


  TCanvas* c13 = new TCanvas();
  c13->cd();
  htrkPPar[0]->GetXaxis()->SetTitle("PPar");
  htrkPPar[0]->Draw("h");
  htrkPPar[1]->Draw("h, same");
  htrkPPar[2]->Draw("h, same");
  legAll->SetHeader("");
  legAll->Draw("same");
  c13->Print("outputAnalysis/plotsCompare_noSmear/htrkPPar.png", "png");
  c13->Print("outputAnalysis/plotsCompare_noSmear/htrkPPar.pdf", "pdf");
  c13->Print("outputAnalysis/plotsCompare_noSmear/htrkPPar.root", "root");

  TCanvas* c14 = new TCanvas();
  c14->cd();
  hGenPtrkPPar[0]->GetXaxis()->SetTitle("PPar");
  hGenPtrkPPar[0]->Draw("h");
  hGenPtrkPPar[1]->Draw("h, same");
  hGenPtrkPPar[2]->Draw("h, same");
  legAll->SetHeader("genTrack");
  legAll->Draw("same");
  c14->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkPPar.png", "png");
  c14->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkPPar.pdf", "pdf");
  c14->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkPPar.root", "root");


  TCanvas* c15 = new TCanvas();
  c15->cd();
  hGenPtrkSumJetEtRatio[0]->GetXaxis()->SetTitle("#Sigma(trk)_{ET})/ Jet_{ET}");
  hGenPtrkSumJetEtRatio[0]->Draw("h");
  hGenPtrkSumJetEtRatio[1]->Draw("h, same");
  hGenPtrkSumJetEtRatio[2]->Draw("h, same");
  legAll->SetHeader("genTrack");
  legAll->Draw("same");
  c15->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkSumJetEtRatio.png", "png");
  c15->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkSumJetEtRatio.pdf", "pdf");
  c15->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkSumJetEtRatio.root", "root");

  TCanvas* c16 = new TCanvas();
  gPad->SetLogy();
  c16->cd();
  hGenPtrkSumJetDeltaR[0]->GetXaxis()->SetTitle("DeltaR(#Sigma(trk), Jet)");
  hGenPtrkSumJetDeltaR[0]->Draw("h");
  hGenPtrkSumJetDeltaR[1]->Draw("h, same");
  hGenPtrkSumJetDeltaR[2]->Draw("h, same");
  legAll->SetHeader("genTrack");
  legAll->Draw("same");
  c16->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkSumJetDeltaR.png", "png");
  c16->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkSumJetDeltaR.pdf", "pdf");
  c16->Print("outputAnalysis/plotsCompare_noSmear/hGenPtrkSumJetDeltaR.root", "root");


}

//  LocalWords:  cd
