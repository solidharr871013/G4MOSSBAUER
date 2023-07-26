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
/// \file BMSPrimaryGeneratorAction.cc
/// \brief Implementation of the BMSPrimaryGeneratorAction class

#include "BMSPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "random"
#include "G4GenericMessenger.hh"
#include "cmath"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSPrimaryGeneratorAction::BMSPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fMessenger(nullptr),
  fParticleGun(nullptr),
  velocity(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  DefineCommands();

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(3000*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSPrimaryGeneratorAction::~BMSPrimaryGeneratorAction()
{
    delete  fMessenger;
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BMSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.






    //the following code is for generating two seperated source, the source is aparted in
    // the distance of source_divergence in x-axis.


    //G4int randNumber = rand()%2;
    G4double randNumber = G4UniformRand();

/**********ring************
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);

    G4double cellNo = 20, cellSide = 2.*mm;
    G4double source_x = 0*cm, source_y = 0*cm, source_z = -50*cm;
    G4double detectDim = 0.5*cellNo*cellSide;

    G4double targetX = detectDim*(2*G4UniformRand()-1),
             targetY = detectDim*(2*G4UniformRand()-1),
             targetZ = 25*mm;

    G4double radius = 10.5*cm;
    source_x = radius*std::cos(G4UniformRand()*2*pi);
    source_y = radius*std::sin(G4UniformRand()*2*pi);
*********ring end***********/

//    G4double cellNo = 16, cellSide = 2.*mm;

//    G4double detectDim = 0.5*cellNo*cellSide;

//    G4double targetX = detectDim*(4*G4UniformRand()-2),
//             targetY = detectDim*(4*G4UniformRand()-2),
//             targetZ = -27*mm;

//    G4double source_x = 20*cm, source_y = 0*cm, source_z = -20*cm;
//    G4double phi = pi+0.01111*pi*(2*G4UniformRand()-1),//45 degree
//             cos_phi = std::cos(phi),
//             sin_phi = std::sin(phi);
//    G4double theta = 0.7853981634+0.1745329252*(2*G4UniformRand()-1),
//             cos_theta = cos(theta),
//             sin_theta = sin(theta);

//             cos_theta = 0.8*(1-0.25*G4UniformRand()),
//             sin_theta = std::sqrt(1-cos_theta*cos_theta);

//     G4double source_x = -11.547*cm, source_y = 0*cm, source_z = -20*cm;
//     G4double phi = 0.0872664626*(2*G4UniformRand()-1),//30 degree
//              cos_phi = std::cos(phi),
//              sin_phi = std::sin(phi);
//     G4double theta = 0.5235987756+0.1745329252*(2*G4UniformRand()-1),
//              cos_theta = std::cos(theta),
//              sin_theta = std::sin(theta);

//              cos_theta = 0.9*(1-0.233333*G4UniformRand()),
//              sin_theta = std::sqrt(1-cos_theta*cos_theta);

//    G4double source_x = -34.641*cm, source_y = 0*cm, source_z = -20*cm;
//    G4double phi = 1/180*pi*(2*G4UniformRand()-1),//60 degree
//             cos_phi = std::cos(phi),
//             sin_phi = std::sin(phi);
//    G4double theta = 1.047197551+0.1745329252*(2*G4UniformRand()-1),
//             cos_theta = cos(theta),
//             sin_theta = sin(theta);

//             cos_theta = 0.6*(1-0.3333*G4UniformRand()),
//             sin_theta = std::sqrt(1-cos_theta*cos_theta);

//    G4double source_x = 0*cm, source_y = 0*cm, source_z = -100*cm;
    G4double phi = 2*pi*G4UniformRand(),//0 degree
             cos_phi = cos(phi),
             sin_phi = sin(phi);
    G4double theta = 0.5*pi*G4UniformRand(),
             cos_theta = cos(theta),
             sin_theta = sin(theta);

    //velocity = 299792458.0;

//    G4double cos_theta = 1-0.001*G4UniformRand(),
//             sin_theta = std::sqrt(1-cos_theta*cos_theta);
/**************multi-source*************************/

    G4double x_direction = sin_theta*cos_phi,
             y_direction = sin_theta*sin_phi,
             z_direction = cos_theta;

    G4double holderinsideThickness = 0.0023*mm, holderinsideRmax = 4*mm, sToOL = 9.65*mm;

    G4double phi2 = 2*pi*G4UniformRand(),//0 degree
             cos_phi2 = cos(phi2),
             sin_phi2 = sin(phi2);

    G4double sourceR = holderinsideRmax*G4UniformRand();



    G4double source_x = sourceR*cos_phi2, 
             source_y = sourceR*sin_phi2, 
             source_z = sToOL+holderinsideThickness*G4UniformRand();

//////////////////////////////////////////////////////
/////////  iron energy  //////////////////////////////
///////////////////////////////////////////////////////

    G4double E122 = 85.6, E136 = 10.68, E14 = 9.16;
    //G4double E122 = 0, E136 = 0, E14 = 100;
    G4double P_122 = E122/(E122+E136+E14), P_136 = E136/(E122+E136+E14);
    G4double velocityC = 299792458.0;

    if(randNumber<P_122){

        fParticleGun->SetParticleEnergy(122.06065*(1+(velocity)*0.001/velocityC)*cos_theta*keV);
    }

    else if(randNumber>=P_122 && randNumber<(P_122+P_136)){


        fParticleGun->SetParticleEnergy(136.47356*(1+(velocity)*0.001/velocityC)*cos_theta*keV);
    }

    else if(randNumber>=(P_122+P_136) && randNumber<=1){

       G4double randEmit = G4UniformRand(), sourceRecoilless = 0.8;

       if(randEmit < sourceRecoilless){
           fParticleGun->SetParticleEnergy(14.4129*(1+(velocity)*0.001/velocityC)*cos_theta*keV);
       }
       else{fParticleGun->SetParticleEnergy(14.4*keV);}
       
    }

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

/////////////////////////////////////////////////////
///////////// tin energy ////////////////////////////
///////////////////////////////////////////////////// 
    // G4double velocityC = 299792458.0;

    // fParticleGun->SetParticleEnergy(23.871*(1+(velocity)*0.001/velocityC)*keV);

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));


    // G4double randEmit = G4UniformRand(), sourceRecoilless = 0.75;

    // if(randEmit < sourceRecoilless){
    //     fParticleGun->SetParticleEnergy(14.4129*keV);
    // }
    // else{fParticleGun->SetParticleEnergy(14.4*keV);}
    
    // fParticleGun->SetParticleEnergy(122.06065*keV);
    // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));

    fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,-1*source_z));

    fParticleGun->GeneratePrimaryVertex(anEvent);



/*************multi-source end**********************/

/********* single source ************



    G4double x_direction = sin_theta*cos_phi,
             y_direction = sin_theta*sin_phi,
             z_direction = cos_theta;
//    G4double x_direction = targetX-source_x,
//             y_direction = targetY-source_y,
//             z_direction = targetZ-source_z;

    //fParticleGun->SetParticleEnergy(300*keV);

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));

    fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

    fParticleGun->GeneratePrimaryVertex(anEvent);

**********single source end******************/

}


void BMSPrimaryGeneratorAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/BMSCamera/setVelocity/",
                                        "gun Control");


    auto& DeltaECmd
      = fMessenger->DeclareMethod("velocity",
                                  &BMSPrimaryGeneratorAction::setDeltaE,
                                  "Set mossbauer Energy.");
    DeltaECmd.SetParameterName("velocity", true);
    //DeltaECmd.SetRange("velocity>=-20. && velocity<=20.");
    DeltaECmd.SetDefaultValue("0");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

