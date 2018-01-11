/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/myTest.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <iostream>
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"

//////////
#include "ExRootAnalysis/ExRootResult.h"

#include "ExRootAnalysis/ExRootUtilities.h"

#include "TROOT.h"
#include "TFile.h"
#include "TClass.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TH2.h"
#include "THStack.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TFolder.h"
///////////////////////////


#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1* hJetPt;
  TH1* hJetEta;

  TH1* hNCharged;
  TH1* hNNeutrals;
  TH1* hFlavor;

  // TH1* htrkEtaRel;
  // TH1* htrkPtRel;

  TH1* htrackSumJetEtRatio;
  TH1* htrackSumJetDeltaR;


  TH1* htrackSumJetEtRatio_bFlav;
  TH1* htrackSumJetEtRatio_cFlav;
  TH1* htrackSumJetEtRatio_lightFlav;

  TH1* htrackSumJetDeltaR_bFlav;
  TH1* htrackSumJetDeltaR_cFlav;
  TH1* htrackSumJetDeltaR_lightFlav;


  TH1* htrkIP2Dxy_bFlav;
  TH1* htrkPtRel_bFlav;
  TH1* htrkEtaRel_bFlav;
  TH1* htrkPPar_bFlav;

  TH1* htrkIP2Dxy_cFlav;
  TH1* htrkPtRel_cFlav;
  TH1* htrkEtaRel_cFlav;
  TH1* htrkPPar_cFlav;

  TH1* htrkIP2Dxy_lightFlav;
  TH1* htrkPtRel_lightFlav;
  TH1* htrkEtaRel_lightFlav;
  TH1* htrkPPar_lightFlav;


  TH1* htrkIP2Dxy_bFlav_no0;
  TH1* htrkIP2Dxy_cFlav_no0;
  TH1* htrkIP2Dxy_lightFlav_no0;

  TH1* htrkIP2Dxy_bFlav_m;
  TH1* htrkIP2Dxy_cFlav_m;
  TH1* htrkIP2Dxy_lightFlav_m;

  TH1* htrkIP3D_bFlav_m;
  TH1* htrkIP3D_cFlav_m;
  TH1* htrkIP3D_lightFlav_m;

  TH1* htrkIP2Dxy_bFlav_80;
  TH1* htrkPtRel_bFlav_80;

  TH1* htrkIP2Dxy_cFlav_80;
  TH1* htrkPtRel_cFlav_80;

  TH1* htrkIP2Dxy_lightFlav_80;
  TH1* htrkPtRel_lightFlav_80;


  //////
  TH1* hGenPtrkPtRel_bFlav;
  TH1* hGenPtrkPtRel_cFlav;
  TH1* hGenPtrkPtRel_lightFlav;
  
  TH1* hGenPtrkIP2Dxy_bFlav;
  TH1* hGenPtrkIP2Dxy_cFlav;
  TH1* hGenPtrkIP2Dxy_lightFlav;

  TH1* hGenPtrackSumJetEtRatio_bFlav;
  TH1* hGenPtrackSumJetEtRatio_cFlav;
  TH1* hGenPtrackSumJetEtRatio_lightFlav;

  TH1* hGenPtrackSumJetDeltaR_bFlav;
  TH1* hGenPtrackSumJetDeltaR_cFlav;
  TH1* hGenPtrackSumJetDeltaR_lightFlav;

  TH1* hGenPtrkPPar_bFlav;
  TH1* hGenPtrkPPar_cFlav;
  TH1* hGenPtrkPPar_lightFlav;
};



//float etaRel(TLorentzVector* dir, TLorentzVector track)
float etaRel(const TVector3 &dir, const TVector3 &track)
{

  TDatabasePDG *pdg;
  pdg = TDatabasePDG::Instance();
  TParticlePDG* pion = pdg->GetParticle(Int_t(211));

  double momPar = dir.Dot(track);
  double energy = std::sqrt(track.Mag2() + pion->Mass()*pion->Mass() );
  
  return 0.5 * std::log((energy + momPar) / (energy - momPar));
}


//float ptRel(TLorentzVector* dir, TLorentzVector track)
float ptRel(const TVector3 &dir, const TVector3 &track)
{

  float magnitute = track.Perp(dir.Unit());
  return magnitute;

  //return track.Cross(dir.Unit()).Mag();

  //  float magnitude = (dir.Unit()).Dot(track);
  //TVector3 projection = dir.Unit() * magnitude;
  //return (track - projection).Mag();
}

float pPar(const TVector3 &dir, const TVector3 &track)
{
  float magnitude = (dir.Unit()).Dot(track);
  return magnitude;
}


float deltaR(float eta1, float eta2, float phi1, float phi2){

  float dEta = std::abs(eta1 - eta2);
  float dPhi = (std::abs(phi1 - phi2) > TMath::Pi()) ? (std::abs(phi1 - phi2) - TMath::Pi()) : std::abs(phi1 - phi2);

  return sqrt(dEta*dEta + dPhi*dPhi);
}

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  /*
  float minD0value = -1.;
  float maxD0value = 1.;
  int D0bins = 600;
  */
  float minD0value = -0.3;
  float maxD0value = 0.3;
  int D0bins = 600;

  float scaleD0value = 30.;
  int D0binsScale = 1000;
  /////////////////////////////////////////////////////////////////
  plots->hJetPt = result->AddHist1D(
   "hJetPt","hJetPt","p_{T} [GeV]", "N_{jet}", 100,0.,100.);

  plots->hJetEta = result->AddHist1D(
   "hJetEta","hJetEta","#eta", "N_{jet}", 250,-5.,5.);

  plots->hNCharged = result->AddHist1D(
   "hNCharged","hNCharged"," NCharged ", "N_{jet}", 100,0.,100.);

  plots->hNNeutrals = result->AddHist1D(
   "hNNeutrals","hNNeutrals"," NNeutral ", "N_{jet}", 100,0.,100.);

  plots->hFlavor = result->AddHist1D(
   "hFlavor","hFlavor"," flavor ", "N_{jet}", 30,0.,30.);

  ///
  // plots->htrkEtaRel = result->AddHist1D(
  //  "htrkEtaRel","htrkEtaRel","trk etaRel ", "N_{trk in jet}", 250,-10.,10.);

  // plots->htrkPtRel = result->AddHist1D(
  //  "htrkPtRel","htrkPtRel","trk PtRel ", "N_{trk in jet}", 100,0.,10.);

  plots->htrackSumJetEtRatio = result->AddHist1D(
   "htrackSumJetEtRatio", "htrackSumJetEtRatio", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);
  plots->htrackSumJetDeltaR = result->AddHist1D(
   "htrackSumJetDeltaR", "htrackSumJetDeltaR", "trackSumJetDeltaR", "N_{trk in jet}", 200, 0., 10.);

  ///by flavor
  plots->htrackSumJetEtRatio_bFlav = result->AddHist1D(
   "htrackSumJetEtRatio_bFlav", "htrackSumJetEtRatio_bFlav", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);
  plots->htrackSumJetEtRatio_cFlav = result->AddHist1D(
   "htrackSumJetEtRatio_cFlav", "htrackSumJetEtRatio_cFlav", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);
  plots->htrackSumJetEtRatio_lightFlav = result->AddHist1D(
   "htrackSumJetEtRatio_lightFlav", "htrackSumJetEtRatio_lightFlav", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);

  plots->htrackSumJetDeltaR_bFlav = result->AddHist1D(
   "htrackSumJetDeltaR_bFlav", "htrackSumJetDeltaR_bFlav", "trackSumJetDeltaR_bFlav", "N_{trk in jet}", 200, 0., 10.);
  plots->htrackSumJetDeltaR_cFlav = result->AddHist1D(
   "htrackSumJetDeltaR_cFlav", "htrackSumJetDeltaR_cFlav", "trackSumJetDeltaR_bFlav", "N_{trk in jet}", 200, 0., 10.);
  plots->htrackSumJetDeltaR_lightFlav = result->AddHist1D(
   "htrackSumJetDeltaR_lightFlav", "htrackSumJetDeltaR_lightFlav", "trackSumJetDeltaR_bFlav", "N_{trk in jet}", 200, 0., 10.);


  // ptRel
  plots->htrkPtRel_bFlav = result->AddHist1D(
   "htrkPtRel_bFlav","htrkPtRel_bFlav","trk PtRel flavour b", "N_{trk in jet}", 100,0.,10.);

  plots->htrkPtRel_cFlav = result->AddHist1D(
   "htrkPtRel_cFlav","htrkPtRel_cFlav","trk PtRel flavour c", "N_{trk in jet}", 100,0.,10.);

  plots->htrkPtRel_lightFlav = result->AddHist1D(
   "htrkPtRel_lightFlav","htrkPtRel_lightFlav","trk PtRel flavour light", "N_{trk in jet}", 100,0.,10.);
  ///////// D0
  plots->htrkIP2Dxy_bFlav = result->AddHist1D(
   "htrkIP2Dxy_bFlav","htrkIP2Dxy_bFlav"," trk IP2D flavour b ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->htrkIP2Dxy_cFlav = result->AddHist1D(
  "htrkIP2Dxy_cFlav","htrkIP2Dxy_cFlav"," trk IP2D flavour c ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->htrkIP2Dxy_lightFlav = result->AddHist1D(
   "htrkIP2Dxy_lightFlav","htrkIP2Dxy_lightFlav"," trk IP2D flavour light ", "N_{trk in jet}", D0bins,minD0value,maxD0value);
  //////////etarel
  plots->htrkEtaRel_bFlav = result->AddHist1D(
   "htrkEtaRel_bFlav","htrkEtaRel_bFlav","trk EtaRel flavour b", "N_{trk in jet}", 200,-5.,5.);

  plots->htrkEtaRel_cFlav = result->AddHist1D(
   "htrkEtaRel_cFlav","htrkEtaRel_cFlav","trk EtaRel flavour c", "N_{trk in jet}", 200,-5.,5.);

  plots->htrkEtaRel_lightFlav = result->AddHist1D(
   "htrkEtaRel_lightFlav","htrkEtaRel_lightFlav","trk EtaRel flavour light", "N_{trk in jet}", 200,-5.,5.);
  //ppar
  plots->htrkPPar_bFlav = result->AddHist1D(
   "htrkPPar_bFlav","htrkPPar_bFlav","trk PPar flavour b", "N_{trk in jet}", 300,-30.,30.);

  plots->htrkPPar_cFlav = result->AddHist1D(
   "htrkPPar_cFlav","htrkPPar_cFlav","trk PPar flavour c", "N_{trk in jet}", 300,-30.,30.);

  plots->htrkPPar_lightFlav = result->AddHist1D(
   "htrkPPar_lightFlav","htrkPtRel_lightFlav","trk PPar flavour light", "N_{trk in jet}", 300,-30.,30.);


  /////////
  plots->htrkIP2Dxy_bFlav_no0 = result->AddHist1D(
   "htrkIP2Dxy_bFlav_no0","htrkIP2Dxy_bFlav_no0"," trk IP2D flavour b ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->htrkIP2Dxy_cFlav_no0 = result->AddHist1D(
  "htrkIP2Dxy_cFlav_no0","htrkIP2Dxy_cFlav_no0"," trk IP2D flavour c ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->htrkIP2Dxy_lightFlav_no0 = result->AddHist1D(
   "htrkIP2Dxy_lightFlav_no0","htrkIP2Dxy_lightFlav_no0"," trk IP2D flavour light ", "N_{trk in jet}", D0bins,minD0value,maxD0value);
  ///computed by me
  plots->htrkIP2Dxy_bFlav_m = result->AddHist1D(
   "htrkIP2Dxy_bFlav_m","htrkIP2Dxy_bFlav_m"," trk IP2D flavour b ", "N_{trk in jet}", D0binsScale,scaleD0value*minD0value,scaleD0value*maxD0value);

  plots->htrkIP2Dxy_cFlav_m = result->AddHist1D(
  "htrkIP2Dxy_cFlav_m","htrkIP2Dxy_cFlav_m"," trk IP2D flavour c ", "N_{trk in jet}", D0binsScale,scaleD0value*minD0value,scaleD0value*maxD0value);

  plots->htrkIP2Dxy_lightFlav_m = result->AddHist1D(
   "htrkIP2Dxy_lightFlav_m","htrkIP2Dxy_lightFlav_m"," trk IP2D flavour light ", "N_{trk in jet}",D0binsScale,scaleD0value*minD0value,scaleD0value*maxD0value);
  ///// //////// ///// ///// 
  plots->htrkIP3D_bFlav_m = result->AddHist1D(
   "htrkIP3D_bFlav_m","htrkIP3D_bFlav_m"," trk IP2D flavour b ", "N_{trk in jet}",D0binsScale,scaleD0value*minD0value,scaleD0value*maxD0value);

  plots->htrkIP3D_cFlav_m = result->AddHist1D(
  "htrkIP3D_cFlav_m","htrkIP3D_cFlav_m"," trk IP2D flavour c ", "N_{trk in jet}",D0binsScale,scaleD0value*minD0value,scaleD0value*maxD0value);

  plots->htrkIP3D_lightFlav_m = result->AddHist1D(
   "htrkIP3D_lightFlav_m","htrkIP3D_lightFlav_m"," trk IP2D flavour light ", "N_{trk in jet}",D0binsScale,scaleD0value*minD0value,scaleD0value*maxD0value);


  //
  plots->htrkPtRel_bFlav_80 = result->AddHist1D(
   "htrkPtRel_bFlav_80","htrkPtRel_bFlav_80","trk PtRel flavour b", "N_{trk in jet}", 100,0.,10.);

  plots->htrkPtRel_cFlav_80 = result->AddHist1D(
   "htrkPtRel_cFlav_80","htrkPtRel_cFlav_80","trk PtRel flavour c", "N_{trk in jet}", 100,0.,10.);

  plots->htrkPtRel_lightFlav_80 = result->AddHist1D(
   "htrkPtRel_lightFlav_80","htrkPtRel_lightFlav_80","trk PtRel flavour light", "N_{trk in jet}", 100,0.,10.);
  /////////
  plots->htrkIP2Dxy_bFlav_80 = result->AddHist1D(
   "htrkIP2Dxy_bFlav_80","htrkIP2Dxy_bFlav_80"," trk IP2D flavour b ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->htrkIP2Dxy_cFlav_80 = result->AddHist1D(
   "htrkIP2Dxy_cFlav_80","htrkIP2Dxy_cFlav_80"," trk IP2D flavour c ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->htrkIP2Dxy_lightFlav_80 = result->AddHist1D(
   "htrkIP2Dxy_lightFlav_80","htrkIP2Dxy_lightFlav_80"," trk IP2D flavour light ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  ////////////////////////////////////////////////////////////////////////

  plots->hGenPtrkPtRel_bFlav = result->AddHist1D(
   "hGenPtrkPtRel_bFlav","hGenPtrkPtRel_bFlav","trk PtRel flavour b", "N_{trk in jet}", 100,0.,10.);

  plots->hGenPtrkPtRel_cFlav = result->AddHist1D(
   "hGenPtrkPtRel_cFlav","hGenPtrkPtRel_cFlav","trk PtRel flavour c", "N_{trk in jet}", 100,0.,10.);

  plots->hGenPtrkPtRel_lightFlav = result->AddHist1D(
   "hGenPtrkPtRel_lightFlav","hGenPtrkPtRel_lightFlav","trk PtRel flavour light", "N_{trk in jet}", 100,0.,10.);
  /////////
  
  plots->hGenPtrkIP2Dxy_bFlav = result->AddHist1D(
   "hGenPtrkIP2Dxy_bFlav","hGenPtrkIP2Dxy_bFlav"," trk IP2D flavour b ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->hGenPtrkIP2Dxy_cFlav = result->AddHist1D(
   "hGenPtrkIP2Dxy_cFlav","hGenPtrkIP2Dxy_cFlav"," trk IP2D flavour c ", "N_{trk in jet}", D0bins,minD0value,maxD0value);

  plots->hGenPtrkIP2Dxy_lightFlav = result->AddHist1D(
   "hGenPtrkIP2Dxy_lightFlav","hGenPtrkIP2Dxy_lightFlav"," trk IP2D flavour light ", "N_{trk in jet}", D0bins,minD0value,maxD0value);
  ///
  plots->hGenPtrackSumJetEtRatio_bFlav = result->AddHist1D(
   "hGenPtrackSumJetEtRatio_bFlav", "hGenPtrackSumJetEtRatio_bFlav", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);
  plots->hGenPtrackSumJetEtRatio_cFlav = result->AddHist1D(
   "hGenPtrackSumJetEtRatio_cFlav", "hGenPtrackSumJetEtRatio_cFlav", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);
  plots->hGenPtrackSumJetEtRatio_lightFlav = result->AddHist1D(
   "hGenPtrackSumJetEtRatio_lightFlav", "hGenPtrackSumJetEtRatio_lightFlav", "trackSumJetEtRatio", "N_{trk in jet}", 200, 0., 3.);

  plots->hGenPtrackSumJetDeltaR_bFlav = result->AddHist1D(
   "hGenPtrackSumJetDeltaR_bFlav", "hGenPtrackSumJetDeltaR_bFlav", "trackSumJetDeltaR_bFlav", "N_{trk in jet}", 200, 0., 10.);
  plots->hGenPtrackSumJetDeltaR_cFlav = result->AddHist1D(
   "hGenPtrackSumJetDeltaR_cFlav", "hGenPtrackSumJetDeltaR_cFlav", "trackSumJetDeltaR_bFlav", "N_{trk in jet}", 200, 0., 10.);
  plots->hGenPtrackSumJetDeltaR_lightFlav = result->AddHist1D(
   "hGenPtrackSumJetDeltaR_lightFlav", "hGenPtrackSumJetDeltaR_lightFlav", "trackSumJetDeltaR_bFlav", "N_{trk in jet}", 200, 0., 10.);
  
  //ppar
  plots->hGenPtrkPPar_bFlav = result->AddHist1D(
   "hGenPtrkPPar_bFlav","hGenPtrkPPar_bFlav","trk PPar flavour b", "N_{trk in jet}", 300,-30.,30.);

  plots->hGenPtrkPPar_cFlav = result->AddHist1D(
   "hGenPtrkPPar_cFlav","hGenPtrkPPar_cFlav","trk PPar flavour c", "N_{trk in jet}", 300,-30.,30.);

  plots->hGenPtrkPPar_lightFlav = result->AddHist1D(
   "hGenPtrkPPar_lightFlav","hGenPtrkPPar_lightFlav","trk PPar flavour light", "N_{trk in jet}", 300,-30.,30.);
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  // TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  // TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  // TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  GenParticle *gentrack;
  Electron *electron;
  Photon *photon;
  Muon *muon;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  Candidate* candidate;
  TLorentzVector momentum;
  TLorentzVector genmomentum;
  //  Float_t DzAll;
  //  Float_t trkSumPt;
  Float_t trkSumEt;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry){
  //for(entry = 0; entry < 10; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // cout << "--  New event -- " << endl;

    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i){
      jet = (Jet*) branchJet->At(i);

      plots->hJetPt->Fill(jet->PT);
      plots->hJetEta->Fill(jet->Eta);

      plots->hFlavor->Fill(jet->Flavor);
      plots->hNCharged->Fill(jet->NCharged);
      plots->hNNeutrals->Fill(jet->NNeutrals);

      momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
      genmomentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
      trkSumEt = 0.;

      float jpx = ((Candidate*)jet)->Momentum.Px();
      float jpy = ((Candidate*)jet)->Momentum.Py();
      float jpz = ((Candidate*)jet)->Momentum.Pz();

      // cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;

      // Loop over all jet's constituents
      Track* highestPtConst;
      float highestPtConstpT = 0.;
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);

        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle:: Class()){
	  //	  std::cout << " const e' Gen "<< std::endl;
	}
        else if(object->IsA() == Track::Class())
        {
	  //	  std::cout << " const e' Track "<< std::endl;
          track = (Track*) object;
	  gentrack = (GenParticle*) track->Particle.GetObject();
	  //cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << " D0 = " << track->D0 << "  D0error = " << track->ErrorD0 << endl;
          momentum += track->P4();
          genmomentum += gentrack->P4();
	  

	  float trackPt = track->PT;
	  //	  std::cout << " trackPt = " << trackPt << std::endl;

	  if(trackPt >= 1. && std::abs(track->Eta) < 2.5){

	  float d0 = track->D0;
	  float gend0 = gentrack->D0;

	  // if(d0 == 0.) std::cout << " jet->Flavor = " << jet->Flavor << " track->PT = " << track->PT << " track->Eta = " << track->Eta 
	  // 			 << " jet->PT = " << jet->PT << " jet->Eta = " << jet->Eta << std::endl;

	  if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21){
	    plots->htrkIP2Dxy_lightFlav->Fill(d0);
	    plots->hGenPtrkIP2Dxy_lightFlav->Fill(gend0);
	  }
	  if(jet->Flavor == 4){
	    plots->htrkIP2Dxy_cFlav->Fill(d0);
	    plots->hGenPtrkIP2Dxy_cFlav->Fill(gend0);
	  }
	  if(jet->Flavor == 5){
	    plots->htrkIP2Dxy_bFlav->Fill(d0);
	    plots->hGenPtrkIP2Dxy_bFlav->Fill(gend0);
	  }
	  
	  if(d0 != 0.){
	    if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21)    plots->htrkIP2Dxy_lightFlav_no0->Fill(d0);
	    if(jet->Flavor == 4)    plots->htrkIP2Dxy_cFlav_no0->Fill(d0);
	    if(jet->Flavor == 5)    plots->htrkIP2Dxy_bFlav_no0->Fill(d0);
	  }


	  //opzione calcolalo da te
	  if(trackPt >= 1. && std::abs(track->Eta) < 2.5){
	    float md0 = TMath::Abs(track->D0);
	    float xd = track->Xd;
	    float yd = track->Yd;
	    float zd = track->Zd;
	    float dd0 = TMath::Abs(track->ErrorD0);
	    float dz = TMath::Abs(track->DZ);
	    float ddz = TMath::Abs(track->ErrorDZ);

	    float sign3D = (jpx*xd + jpy*yd + jpz*zd > 0.0) ? 1 : -1;
	    //add transverse and longitudinal significances in quadrature
	    float sip3D = sign3D * TMath::Sqrt( TMath::Power(d0 / dd0, 2) + TMath::Power(dz / ddz, 2) );

	    float sign2D = (jpx*xd + jpy*yd > 0.0) ? 1 : -1;
	    float sip2D = sign2D * d0 / TMath::Abs(dd0);
	    //float sip2D = sign2D * d0;

	    //std::cout << " sip2D = " << sip2D << " sip3D = " << sip3D << " jet->Flavor = " << jet->Flavor << std::endl;

	    if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21){
	      plots->htrkIP2Dxy_lightFlav_m->Fill(sip2D);
	      plots->htrkIP3D_lightFlav_m->Fill(sip3D);
	    }
	    if(jet->Flavor == 4){
	      plots->htrkIP2Dxy_cFlav_m->Fill(sip2D);
	      plots->htrkIP3D_cFlav_m->Fill(sip3D);
	    }
	    if(jet->Flavor == 5){
	      plots->htrkIP2Dxy_bFlav_m->Fill(sip2D);
	      plots->htrkIP3D_bFlav_m->Fill(sip3D);
	    }
	  }
	  }

	  //opzione track and genTrack
	  // plots->htrkEtaRel->Fill(etaRel( (((Candidate*)jet)->Momentum).Vect(),   track->P4().Vect() ) );
	  // plots->htrkPtRel->Fill(ptRel( (((Candidate*)jet)->Momentum.Vect()),     track->P4().Vect() ) );
	  if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21){
	    plots->htrkPtRel_lightFlav->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(), track->P4().Vect() ) );
	    plots->htrkPPar_lightFlav->Fill(pPar( (((Candidate*)jet)->Momentum).Vect(), track->P4().Vect() ) );
	    plots->htrkEtaRel_lightFlav->Fill(etaRel( (((Candidate*)jet)->Momentum).Vect(), track->P4().Vect() ) );

	    plots->hGenPtrkPtRel_lightFlav->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(), gentrack->P4().Vect() ) );
	    plots->hGenPtrkPPar_lightFlav->Fill(pPar( (((Candidate*)jet)->Momentum).Vect(), gentrack->P4().Vect() ) );
	  }
	  if(jet->Flavor == 4){
	    plots->htrkPtRel_cFlav->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );
	    plots->htrkPPar_cFlav->Fill(pPar( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );
	    plots->htrkEtaRel_cFlav->Fill(etaRel( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );

	    plots->hGenPtrkPtRel_cFlav->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(),     gentrack->P4().Vect() ) );
	    plots->hGenPtrkPPar_cFlav->Fill(pPar( (((Candidate*)jet)->Momentum).Vect(),     gentrack->P4().Vect() ) );
	  }
	  if(jet->Flavor == 5){
	    plots->htrkPtRel_bFlav->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );
	    plots->htrkPPar_bFlav->Fill(pPar( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );
	    plots->htrkEtaRel_bFlav->Fill(etaRel( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );

	    plots->hGenPtrkPtRel_bFlav->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(),     gentrack->P4().Vect() ) );
    	    plots->hGenPtrkPPar_bFlav->Fill(pPar( (((Candidate*)jet)->Momentum).Vect(),     gentrack->P4().Vect() ) );
	  }


	  trkSumEt += track->P4().Et();

	  if(jet->PT > 80 && track->PT >= 1. && std::abs(track->Eta) < 2.5){
	    if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21)       plots->htrkIP2Dxy_lightFlav_80->Fill(track->D0);
	    if(jet->Flavor == 4)    plots->htrkIP2Dxy_cFlav_80->Fill(track->D0);
	    if(jet->Flavor == 5)    plots->htrkIP2Dxy_bFlav_80->Fill(track->D0);

	    //opzione track

	    if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21) 
	      plots->htrkPtRel_lightFlav_80->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(), track->P4().Vect() ) );
	    if(jet->Flavor == 4)            plots->htrkPtRel_cFlav_80->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );
	    if(jet->Flavor == 5)            plots->htrkPtRel_bFlav_80->Fill(ptRel( (((Candidate*)jet)->Momentum).Vect(),     track->P4().Vect() ) );

	  }

        }
        else if(object->IsA() == Tower::Class())
        {
	  //	  std::cout << " const e' Tower "<< std::endl;
          tower = (Tower*) object;
          // cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          //momentum += tower->P4();
        }
      }
  
      plots->htrackSumJetEtRatio->Fill(momentum.Et()/jet->P4().Et());
      plots->htrackSumJetDeltaR->Fill(deltaR(jet->Eta, momentum.Eta(), jet->Phi, momentum.Phi()) );

      if(jet->Flavor == 1 || jet->Flavor == 2 || jet->Flavor == 3 || jet->Flavor == 21){
	plots->htrackSumJetEtRatio_lightFlav->Fill(momentum.Et()/jet->P4().Et());
	plots->htrackSumJetDeltaR_lightFlav->Fill(deltaR(jet->Eta, momentum.Eta(), jet->Phi, momentum.Phi()) );

	plots->hGenPtrackSumJetEtRatio_lightFlav->Fill(genmomentum.Et()/jet->P4().Et());
	plots->hGenPtrackSumJetDeltaR_lightFlav->Fill(deltaR(jet->Eta, genmomentum.Eta(), jet->Phi, genmomentum.Phi()) );
      }
      if(jet->Flavor == 4){
	plots->htrackSumJetEtRatio_cFlav->Fill(momentum.Et()/jet->P4().Et());
	plots->htrackSumJetDeltaR_cFlav->Fill(deltaR(jet->Eta, momentum.Eta(), jet->Phi, momentum.Phi()) );

	plots->hGenPtrackSumJetEtRatio_cFlav->Fill(genmomentum.Et()/jet->P4().Et());
	plots->hGenPtrackSumJetDeltaR_cFlav->Fill(deltaR(jet->Eta, genmomentum.Eta(), jet->Phi, genmomentum.Phi()) );
      }
      if(jet->Flavor == 5){
	plots->htrackSumJetEtRatio_bFlav->Fill(momentum.Et()/jet->P4().Et());
	plots->htrackSumJetDeltaR_bFlav->Fill(deltaR(jet->Eta, momentum.Eta(), jet->Phi, momentum.Phi()) );

	plots->hGenPtrackSumJetEtRatio_bFlav->Fill(genmomentum.Et()/jet->P4().Et());
	plots->hGenPtrackSumJetDeltaR_bFlav->Fill(deltaR(jet->Eta, genmomentum.Eta(), jet->Phi, genmomentum.Phi()) );
      }

    }
  }


  plots->htrkIP2Dxy_bFlav->Scale(1./plots->htrkIP2Dxy_bFlav->Integral());
  plots->htrkIP2Dxy_cFlav->Scale(1./plots->htrkIP2Dxy_cFlav->Integral());
  plots->htrkIP2Dxy_lightFlav->Scale(1./plots->htrkIP2Dxy_lightFlav->Integral());
   //
  plots->htrkIP2Dxy_bFlav_80->Scale(1./plots->htrkIP2Dxy_bFlav_80->Integral());
  plots->htrkIP2Dxy_cFlav_80->Scale(1./plots->htrkIP2Dxy_cFlav_80->Integral());
  plots->htrkIP2Dxy_lightFlav_80->Scale(1./plots->htrkIP2Dxy_lightFlav_80->Integral());

  plots->hGenPtrkIP2Dxy_bFlav->Scale(1./plots->hGenPtrkIP2Dxy_bFlav->Integral());
  plots->hGenPtrkIP2Dxy_cFlav->Scale(1./plots->hGenPtrkIP2Dxy_cFlav->Integral());
  plots->hGenPtrkIP2Dxy_lightFlav->Scale(1./plots->hGenPtrkIP2Dxy_lightFlav->Integral());
   //
  plots->htrkIP2Dxy_bFlav_m->Scale(1./plots->htrkIP2Dxy_bFlav_m->Integral());
  plots->htrkIP2Dxy_cFlav_m->Scale(1./plots->htrkIP2Dxy_cFlav_m->Integral());
  plots->htrkIP2Dxy_lightFlav_m->Scale(1./plots->htrkIP2Dxy_lightFlav_m->Integral());

  //mancano da scalare
  // _no0  == eliminando gli zeri  

  plots->htrkPtRel_bFlav->Scale(1./plots->htrkPtRel_bFlav->Integral());
  plots->htrkPtRel_cFlav->Scale(1./plots->htrkPtRel_cFlav->Integral());
  plots->htrkPtRel_lightFlav->Scale(1./plots->htrkPtRel_lightFlav->Integral());

  plots->htrkEtaRel_bFlav->Scale(1./plots->htrkEtaRel_bFlav->Integral());
  plots->htrkEtaRel_cFlav->Scale(1./plots->htrkEtaRel_cFlav->Integral());
  plots->htrkEtaRel_lightFlav->Scale(1./plots->htrkEtaRel_lightFlav->Integral());

  plots->htrkPPar_bFlav->Scale(1./plots->htrkPPar_bFlav->Integral());
  plots->htrkPPar_cFlav->Scale(1./plots->htrkPPar_cFlav->Integral());
  plots->htrkPPar_lightFlav->Scale(1./plots->htrkPPar_lightFlav->Integral());

  plots->hGenPtrkPPar_bFlav->Scale(1./plots->hGenPtrkPPar_bFlav->Integral());
  plots->hGenPtrkPPar_cFlav->Scale(1./plots->hGenPtrkPPar_cFlav->Integral());
  plots->hGenPtrkPPar_lightFlav->Scale(1./plots->hGenPtrkPPar_lightFlav->Integral());


  plots->htrkPtRel_bFlav_80->Scale(1./plots->htrkPtRel_bFlav_80->Integral());
  plots->htrkPtRel_cFlav_80->Scale(1./plots->htrkPtRel_cFlav_80->Integral());
  plots->htrkPtRel_lightFlav_80->Scale(1./plots->htrkPtRel_lightFlav_80->Integral());

  ////
  plots->hGenPtrkPtRel_bFlav->Scale(1./plots->hGenPtrkPtRel_bFlav->Integral());
  plots->hGenPtrkPtRel_cFlav->Scale(1./plots->hGenPtrkPtRel_cFlav->Integral());
  plots->hGenPtrkPtRel_lightFlav->Scale(1./plots->hGenPtrkPtRel_lightFlav->Integral());

  /////////
  plots->htrackSumJetEtRatio_bFlav->Scale(1./plots->htrackSumJetEtRatio_bFlav->Integral());
  plots->htrackSumJetEtRatio_cFlav->Scale(1./plots->htrackSumJetEtRatio_cFlav->Integral());
  plots->htrackSumJetEtRatio_lightFlav->Scale(1./plots->htrackSumJetEtRatio_lightFlav->Integral());

  plots->htrackSumJetDeltaR_bFlav->Scale(1./plots->htrackSumJetDeltaR_bFlav->Integral());
  plots->htrackSumJetDeltaR_cFlav->Scale(1./plots->htrackSumJetDeltaR_cFlav->Integral());
  plots->htrackSumJetDeltaR_lightFlav->Scale(1./plots->htrackSumJetDeltaR_lightFlav->Integral());

  plots->hGenPtrackSumJetEtRatio_bFlav->Scale(1./plots->hGenPtrackSumJetEtRatio_bFlav->Integral());
  plots->hGenPtrackSumJetEtRatio_cFlav->Scale(1./plots->hGenPtrackSumJetEtRatio_cFlav->Integral());
  plots->hGenPtrackSumJetEtRatio_lightFlav->Scale(1./plots->hGenPtrackSumJetEtRatio_lightFlav->Integral());

  plots->hGenPtrackSumJetDeltaR_bFlav->Scale(1./plots->hGenPtrackSumJetDeltaR_bFlav->Integral());
  plots->hGenPtrackSumJetDeltaR_cFlav->Scale(1./plots->hGenPtrackSumJetDeltaR_cFlav->Integral());
  plots->hGenPtrackSumJetDeltaR_lightFlav->Scale(1./plots->hGenPtrackSumJetDeltaR_lightFlav->Integral());

}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("outputAnalysis/");
}

//------------------------------------------------------------------------------

void myTest(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("outputAnalysis/results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
