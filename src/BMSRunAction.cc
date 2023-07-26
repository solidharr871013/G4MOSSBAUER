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
/// \file BMSRunAction.cc
/// \brief Implementation of the BMSRunAction class

#include "BMSRunAction.hh"
#include "BMSPrimaryGeneratorAction.hh"
#include "BMSDetectorConstruction.hh"
#include "BMSEventAction.hh"
#include "BMSAnalysis.hh"
// #include "BMSRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSRunAction::BMSRunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  f3LayerCoincidence(0),
  fTotalNumber(0),
  fTransmitNumber(0),
  fBMS(0),
  fTMS(0),
  f2LayerCoincidence(0)
{ 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 
  accumulableManager->RegisterAccumulable(f3LayerCoincidence);
  accumulableManager->RegisterAccumulable(fTotalNumber);
  accumulableManager->RegisterAccumulable(fTransmitNumber);
  accumulableManager->RegisterAccumulable(fBMS);
  accumulableManager->RegisterAccumulable(fTMS);
  accumulableManager->RegisterAccumulable(f2LayerCoincidence);

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetVerboseLevel(1);
  analysis->SetFileName("BMSHisto");
  analysis->CreateH1("EdepInBackscattering","E",1500,0,150*keV);//ID=0
  analysis->CreateH1("EdepInTransmitting","E",1500,0,150*keV);//ID=1

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSRunAction::~BMSRunAction()
{delete G4AnalysisManager::Instance();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BMSRunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BMSRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4double edepTest = fEdep.GetValue();
  G4cout << "the test energy value is " << edepTest/keV << G4endl;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  // cast int to double type
  G4double NO3LayerCoincidence = f3LayerCoincidence.GetValue();
  G4double NOTotalEvent = fTotalNumber.GetValue();
  G4double NoTransmit = fTransmitNumber.GetValue();
  G4double NoBMS = fBMS.GetValue();
  G4double NoTMS = fTMS.GetValue();
  G4double No2LayerCoincidence = f2LayerCoincidence.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const BMSDetectorConstruction* detectorConstruction
   = static_cast<const BMSDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  G4double BMS3EventRatio = NO3LayerCoincidence/NOTotalEvent,
           BMS2EventRatio = No2LayerCoincidence/NOTotalEvent;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const BMSPrimaryGeneratorAction* generatorAction
   = static_cast<const BMSPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }



  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << "The number of total event is " << NOTotalEvent << G4endl
     << "The number of transmitting event is " << NoTransmit << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  
  analysis->Write();
  analysis->CloseFile();


if(IsMaster()){

/****************Start txt format output***************/
    std::ofstream SingleBMSEfficiency;


    SingleBMSEfficiency.open(fOutput + "BMS_TMS.txt", std::ios_base::app);
    if (SingleBMSEfficiency.is_open()){

      SingleBMSEfficiency << NoBMS << "; " << NoTMS << ";" << std::endl;
      }
    else SingleBMSEfficiency << "Unable to open file" << std::endl;

/****************End txt format output*********************/

}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BMSRunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void BMSRunAction::Add2LayerCoincidence(G4int number){
    f2LayerCoincidence += number;
}

void BMSRunAction::Add3LayerCoincidence(G4int number){
    f3LayerCoincidence += number;
}

void BMSRunAction::AddBMS(G4int number){
    fBMS += number;
}

void BMSRunAction::AddTMS(G4int number){
    fTMS += number;
}

void BMSRunAction::AddTotalNumber(G4int totalnumber){
    fTotalNumber += totalnumber;
}

void BMSRunAction::AddTransmitNumber(G4int transmitNumber){
  fTransmitNumber += transmitNumber;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

