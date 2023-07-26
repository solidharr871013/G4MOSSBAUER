//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file BMSEventAction.cc
/// \brief Implementation of the BMSEventAction class

#include "BMSEventAction.hh"
#include "BMSRunAction.hh"
#include "BMSSteppingAction.hh"
#include "BMSAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4THitsMap.hh"

#include "BMSLayer1Hit.hh"
#include "BMSLayer1SD.hh"

#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception("BMSEventAction::EndOfEventAction()",
                  "BMSCode001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception("BMSEventAction::EndOfEventAction()",
                "BMSCode001", JustWarning, msg);
  }
  return hc;
}
}

BMSEventAction::BMSEventAction(BMSRunAction* runAction)
: G4UserEventAction(),
  fMessenger(nullptr),
  fRunAction(runAction),
  fEdep(0.),
  fTotalNumber(0),
  fTransmitNumber(0),
  fBMS(0),
  fTMS(0),
  fLayer1ID(-1),
  fLayer2ID(-1),
  fLayer3ID(-1),
  f2TotalE_low(0.6330),
  f2TotalE_high(0.693),
  f2ScatterE_low(0.020),
  f2ScatterE_high(0.16)
{
    DefineCommands();
     G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSEventAction::~BMSEventAction()
{ delete fMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BMSEventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  fTotalNumber = 0;
  fTransmitNumber = 0;
  fBMS = 0;
  fTMS = 0;

  if(fLayer1ID == -1 ){
      auto sdManager = G4SDManager::GetSDMpointer();

      fLayer1ID = sdManager->GetCollectionID("Layer1/LayerColl");
      fLayer2ID = sdManager->GetCollectionID("Layer2/LayerColl");
      fLayer3ID = sdManager->GetCollectionID("Layer3/LayerColl");
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BMSEventAction::EndOfEventAction(const G4Event* event)
{   
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();

    G4double eThreshold = 1*keV;

    auto hcLayer1 = GetHC(event, fLayer1ID);
    if(!hcLayer1) return;

    auto hcLayer2 = GetHC(event, fLayer2ID);
    if(!hcLayer2) return;

/*****************************************************************/
    //G4int test = BMSEventAction::GetTotalNumber();
    //G4cout << "now the total number is " << test << G4endl;


  // accumulate statistics in run action
   // if(fEdep>eThreshold){
     //  fRunAction->AddEdep(fEdep);
       //analysis->FillH1(0,fEdep);
   // }

    if(hcLayer1){

        G4double edepInLayer1 = 0.;

        for (unsigned long itr=0; itr != hcLayer1->GetSize(); ++itr) {

            auto hit = static_cast<BMSLayer1Hit*>(hcLayer1->GetHit(itr));
            edepInLayer1 += hit->GetEdep();
        }

        G4double EdepBack = AddFluction(edepInLayer1);

        if(edepInLayer1>eThreshold){
            analysis->FillH1(0,EdepBack);
            ++fTotalNumber;
        }
        fRunAction->AddTotalNumber(fTotalNumber);

        if(EdepBack<(17.3)*keV && EdepBack>(11.4)*keV){
            ++fBMS;
        }
        if(EdepBack<8.4*keV && EdepBack>4.5*keV){
            ++fBMS;
        }
        fRunAction->AddBMS(fBMS);
    }

    if(hcLayer2){

        G4double edepInLayer2 = 0.;

        for (unsigned long itr=0; itr != hcLayer2->GetSize(); ++itr) {

            auto hit = static_cast<BMSLayer1Hit*>(hcLayer2->GetHit(itr));
            edepInLayer2 += hit->GetEdep();
        }
        
        G4double EdepTrans = AddFluctionTrans(edepInLayer2);

        if(edepInLayer2>eThreshold){
            analysis->FillH1(1,EdepTrans);
            ++fTransmitNumber;
        }

        fRunAction->AddTransmitNumber(fTransmitNumber);

        if(EdepTrans<(14.4*(1+0.06))*keV && EdepTrans>14.4*(1-0.06)*keV){
            ++fTMS;
        }
        fRunAction->AddTMS(fTMS);
    }

    



}

G4double BMSEventAction::AddFluction(G4double val){

    //calculate the coefficient of the resolution, based on 662keV
    G4double FWHM = 2065, Energy = 5900, Resolution = FWHM/Energy,
             coefficient = Resolution*std::sqrt(Energy/1000000);

    G4double valVar = coefficient*(std::sqrt(val))/(2*std::sqrt(2*std::log(2)));
    //G4double valVar = 0.08*val/(2*std::sqrt(2*std::log(2)));


    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<G4double> distribution(val,valVar);
    G4double var_new = distribution(generator);

    return var_new;
}

G4double BMSEventAction::AddFluctionTrans(G4double val){

    //calculate the coefficient of the resolution, based on 662keV
    G4double FWHM = 5760, Energy = 144000, Resolution = FWHM/Energy,
             coefficient = Resolution*std::sqrt(Energy/1000000);

    G4double valVar = coefficient*(std::sqrt(val))/(2*std::sqrt(2*std::log(2)));
    //G4double valVar = 0.08*val/(2*std::sqrt(2*std::log(2)));


    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<G4double> distribution(val,valVar);
    G4double var_new = distribution(generator);

    return var_new;
}

void BMSEventAction::Set2TotalELow(G4double val){f2TotalE_low = val;}
void BMSEventAction::Set2TotalEHigh(G4double val){f2TotalE_high = val;}
void BMSEventAction::Set2ScatterELow(G4double val){f2ScatterE_low = val;}
void BMSEventAction::Set2ScatterEHigh(G4double val){f2ScatterE_high = val;}



void BMSEventAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/BMSCamera/event/",
                                        "Event Control");


    auto& ThresTotalEHighCmd
      = fMessenger->DeclareMethodWithUnit("ThresTotalEHigh","MeV",
                                  &BMSEventAction::Set2TotalEHigh,
                                  "Set high threshold of total edep in a 2coincident event.");
    ThresTotalEHighCmd.SetParameterName("Energy", true);
    ThresTotalEHighCmd.SetRange("Energy>=0. && Energy<100.");
    ThresTotalEHighCmd.SetDefaultValue("0.710");

    auto& ThresTotalELowCmd
      = fMessenger->DeclareMethodWithUnit("ThresTotalELow","MeV",
                                  &BMSEventAction::Set2TotalELow,
                                  "Set low threshold of total edep in a 2coincident event.");
    ThresTotalELowCmd.SetParameterName("Energy", true);
    ThresTotalELowCmd.SetRange("Energy>=0. && Energy<100.");
    ThresTotalELowCmd.SetDefaultValue("0.615");

    auto& ThresScatterEHighCmd
      = fMessenger->DeclareMethodWithUnit("ThresScatterEHigh","MeV",
                                  &BMSEventAction::Set2ScatterEHigh,
                                  "Set high threshold of Scatter layer edep in a 2coincident event.");
    ThresScatterEHighCmd.SetParameterName("Energy", true);
    ThresScatterEHighCmd.SetRange("Energy>=0. && Energy<100.");
    ThresScatterEHighCmd.SetDefaultValue("0.165");

    auto& ThresScatterELowCmd
      = fMessenger->DeclareMethodWithUnit("ThresScatterELow","MeV",
                                  &BMSEventAction::Set2ScatterELow,
                                  "Set Low threshold of Scatter layer edep in a 2coincident event.");
    ThresScatterELowCmd.SetParameterName("Energy", true);
    ThresScatterELowCmd.SetRange("Energy>=0. && Energy<100.");
    ThresScatterELowCmd.SetDefaultValue("0.01");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
