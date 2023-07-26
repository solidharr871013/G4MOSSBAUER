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
/// \file BMSEventAction.hh
/// \brief Definition of the BMSEventAction class

#ifndef BMSEventAction_h
#define BMSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class BMSRunAction;
class G4GenericMessenger;


/// Event action class
///

class BMSEventAction : public G4UserEventAction
{
  public:
    BMSEventAction(BMSRunAction* runAction);
    virtual ~BMSEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) { fEdep += edep; }
    
    void AddTotalNumber(G4int totalnumber) {fTotalNumber += totalnumber;}
    void AddTransmitNumber(G4int transmitNumber){fTransmitNumber += transmitNumber;}

    void AddBMS(G4int BMSnum){fBMS += BMSnum;}
    void AddTMS(G4int TMSnum){fTMS += TMSnum;}

    void Set2TotalEHigh(G4double val);
    void Set2TotalELow(G4double val);
    void Set2ScatterEHigh(G4double val);
    void Set2ScatterELow(G4double val);

    G4double AddFluction(G4double val);
    G4double AddFluctionTrans(G4double val);

    G4int GetTotalNumber() { return fTotalNumber; }

  private:
    void DefineCommands();
    G4GenericMessenger *fMessenger;

    BMSRunAction* fRunAction;
    G4double     fEdep;

    G4int        fTotalNumber;
    G4int fTransmitNumber;
    G4int fBMS, fTMS;

    G4int fLayer1ID;
    G4int fLayer2ID;
    G4int fLayer3ID;

    G4double f2TotalE_low, f2TotalE_high, f2ScatterE_low, f2ScatterE_high;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
