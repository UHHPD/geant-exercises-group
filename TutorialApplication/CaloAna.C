#include "include/TutorialApplication.hh"
#include "TGeoManager.h"
#include "TVirtualGeoTrack.h"
#include "TDatabasePDG.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <iostream>
#include <string>

Float_t XofFirstSecondary()
{
  TObjArray* tracks = gGeoManager->GetListOfTracks();
  
  if(tracks->GetEntriesFast() < 2) return 999;
  
  //get first secondary
  TVirtualGeoTrack* track = (TVirtualGeoTrack*)tracks->At(1);
  
  //get the first point of this track
  Double_t x,y,z,t;
  track->GetPoint(0,x,y,z,t);
  
  return x;
}



Int_t CountChargedinScint()
{
  Int_t ncharged = 0;

  TObjArray* tracks = gGeoManager->GetListOfTracks();

  for(Int_t i=0,l=tracks->GetEntriesFast();i<l;++i) {
    TVirtualGeoTrack* track = (TVirtualGeoTrack*)tracks->At(i);
    //track->Print();
    
    Int_t pdg = track->GetPDG();
    Double_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    if( charge == 0 ) continue;
    
    Double_t x,y,z,t;
    TGeoNode* lastnode = NULL;
    for(Int_t j=0,lp=track->GetNpoints();j<lp;++j) {
      track->GetPoint(j,x,y,z,t);
      TGeoNode* node = gGeoManager->FindNode(x,y,z);
      if(! node ) continue;
      //node->Print();
      
      if( lastnode == node ) continue;
      lastnode = node;
      //is scintillator ?
      //std::cout << node->GetMedium()->GetMaterial()->GetZ() << std::endl;
      
      if(node->GetMedium()->GetMaterial()->GetZ() == 1)
	++ncharged;
    }
    //std::cout << "charge:" << ncharged << std::endl;
  }
  //std::cout << ncharged << std::endl;
  return ncharged;
}


void CaloAna()
{
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();
// initialize geometry: volumes and materials of a Sampling Calorimeter   
  Double_t AbsWid=.3;         //Absorber width
  Double_t SciWid=1.;         //Scintillator width, 
  Double_t SizeFact=2;      //size of the calorimeter in interaction lengths labmda_I, 4.
  Int_t IMat=1;               //material 1:Pb 2:Fe 
  TString geom("geometry/SamplingCalorimeter(");
  geom+=AbsWid; geom.Append(",");
  geom+=SciWid; geom.Append(",");
  geom+=SizeFact; geom.Append(",");
  geom+=IMat; geom.Append(")");
  app->InitMC(geom);
  // xxx app->InitMC("geometry/SamplingCalorimeter(2.,1.,2.,1)"); 

// Book histogram(s)
 // for shower analyis
  TH1F* hx  = new TH1F("hx","starting point of shower",100, 0.,25.);
  hx->SetYTitle("# of entries");  
  hx->SetXTitle("x of first vertex [cm]");
  TH1F* hwidth  = new TH1F("hwidth","width of the shower",50,0.,25.);
  hwidth->SetYTitle("# of entries");
  hwidth->SetXTitle("width of the shower [cm]");
  TH1F* hlength = new TH1F("hlength","length of the shower",50,0.,25.);
  hlength->SetYTitle("# of entries");
  hlength->SetXTitle("length of the shower [cm]");
  // for hit counting
  TProfile* hcounts = new TProfile("hcounts","Counts vs particle energy",
				   20,0.,10.,"s");
  TH2F* energyresolution= new TH2F("energyresolution","Energy Resolution",100,0.,10.,100,-1.5,1.5);
  energyresolution->SetXTitle("Energy in GeV");
  energyresolution->SetYTitle("/Delta E/E");
  // option "s": show sigma(i) instead of sigma(i)/sqrt(n_i)
  hcounts->SetXTitle("energy [GeV]");
  hcounts->SetYTitle("mean number of counts");
  TH2D* hresponse = new TH2D("hresponse","measured energy/particle energy vs particle energy; energy [GeV]; response",
			     20,0.,10.,50,0,2);
  TH1F* e_pi= new TH1F("e/pi","e/pi",100,0,4);
  e_pi->SetYTitle("counts");
  e_pi->SetXTitle("e/pi");
  std::vector<float> e;
  std::vector<float> pi;

//simulate events at fixed momentum
  TH1F* hhelp; // for analysis of internal histograms
  Double_t xp[1]={0.90},xq[1];

  unsigned int nevt = 100;
  double       p = 10;//GeV

  app->SetPrimaryPDG(-11); 
  /* PDG codes     22: Photon    +/-11: e+/-  +-13: muon   
               +/-211: pion    +/-2212: proton              */
  app->SetPrimaryMomentum(p);
  for(unsigned int i = 0 ; i < nevt ; ++i) {
    app->RunMC(1,!(i%100)); 
    // fill starting point of shower (pos. of first secondary)
    hx->Fill(XofFirstSecondary());
    // access GEANT internal histograms
    hhelp = (TH1F*) gROOT->Get("hEdepTrans"); assert(hhelp);
    // and evaluate quantiles x-wise (-> radius)
    hhelp->GetQuantiles(1,xq,xp);
    // fill the width of the event as two times the max. radius
    hwidth->Fill(2.*xq[0]);
    hhelp = (TH1F*) gROOT->Get("hEdepLong"); assert(hhelp);
    hhelp->GetQuantiles(1,xq,xp);
    hlength->Fill(xq[0]);
    e.push_back(CountChargedinScint());
    
    // reset internal histograms
    app->FinishRun();
  }
  
  app->SetPrimaryPDG(-211); 
  /* PDG codes     22: Photon    +/-11: e+/-  +-13: muon   
               +/-211: pion    +/-2212: proton              */
  app->SetPrimaryMomentum(p);
  for(unsigned int i = 0 ; i < nevt ; ++i) {
    app->RunMC(1,!(i%100)); 
    pi.push_back(CountChargedinScint());
    
    // reset internal histograms
    app->FinishRun();
  }
  
  for(int i=0;i<e.size();i++){
  e_pi->Fill(e[i]/pi[i]);
    
  }
  app->SetPrimaryPDG(-11); 
  // events at different momenta
  nevt = 10000; p = 0.1;
  double stepping = 9.9 / nevt;
  // generate a large number of events
  for(unsigned int i=0;i<nevt;++i) {
    app->SetPrimaryMomentum(p);
    app->RunMC(1,!(i%1000));
    hcounts->Fill(p,CountChargedinScint());
    float E_res=(p-(CountChargedinScint()+2.31)/115.19)/p; // for electron
    //float E_res=((p-(CountChargedinScint()+18.39)/60.88))/p; //for pion
    energyresolution->Fill(p,E_res);
    p += stepping;
    
    // reset internal histograms
    app->FinishRun();
  }
  TObjArray aSlices;
  energyresolution->FitSlicesY();
  TH2F *energyresolution_0 = (TH2F*)gDirectory->Get("energyresolution_1");
  // display results  
  TCanvas* c = new TCanvas(); c->Divide(2,2);
  c->cd(1);  hx->Draw();
  c->cd(2);  hwidth->Draw();
  c->cd(3);  hlength->Draw();
  c->cd(4);  hcounts->Draw();
  TCanvas* d = new TCanvas(); d->Divide(1,3);
  d->cd(1);  energyresolution->Draw();
  d->cd(2); energyresolution_0->Draw();
  d->cd(3); e_pi->Draw();
}
