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
/// \file BMSDetectorConstruction.cc
/// \brief Implementation of the BMSDetectorConstruction class

#include "BMSDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4GenericMessenger.hh"

#include "cmath"

#include "BMSLayer1SD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSDetectorConstruction::BMSDetectorConstruction()
: G4VUserDetectorConstruction(),
  fMessenger(nullptr),
  fLogicalWorld(nullptr),
  BMSdetlog(nullptr)
{
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BMSDetectorConstruction::~BMSDetectorConstruction()
{
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* BMSDetectorConstruction::Construct()
{
   //
   // Construct the material to be used
   //
   ConstructMaterial();
   G4Material* world_mat = G4Material::GetMaterial("G4_AIR");
   G4Material* lead_mat = G4Material::GetMaterial("lead");
   //G4Material* sample_mat = G4Material::GetMaterial("G4_STAINLESS-STEEL");
   G4Material* sourceholder_mat = G4Material::GetMaterial("Ti_alloy");
   G4Material* source_mat = G4Material::GetMaterial("sourceMatrix");
   //G4Material* Siboard_mat = G4Material::GetMaterial("Siboard");
   //G4Material* cell_mat = G4Material::GetMaterial("Gd3Ca5Al5O12");
   //G4Material* cell_mat = G4Material::GetMaterial("Lu2SiO5");
   //G4Material* layer_mat = G4Material::GetMaterial("BaSO4");
   //G4Material* sample_mat = G4Material::GetMaterial("G4_LANTHANUM_OXYSULFIDE");
   //G4Material* sample_mat = G4Material::GetMaterial("G4_FERRIC_OXIDE");
   //G4Material* sample_mat = G4Material::GetMaterial("SnSb");
   
   //G4Material* sample_mat = G4Material::GetMaterial("G4_Fe");
   //G4Material* sample_mat = G4Material::GetMaterial("iron");
   G4Material* sample_mat = G4Material::GetMaterial("test_mat");

  //
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  /////////////////////////////////////////////////////////////////////
  G4double world_sizeXY = 10*m;
  G4double world_sizeZ  = 10*m;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  fLogicalWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

                      
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(nullptr,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fLogicalWorld,            //its logical volume
                      "World",               //its name
                      nullptr,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

                      

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  G4double container_side = 10*mm;
  G4ThreeVector container_pos = G4ThreeVector(0,0,-202*mm);

  G4RotationMatrix* layer_rotation = new G4RotationMatrix();
  layer_rotation->rotateY(0.*deg);

  G4Box* containerSolid = new G4Box("containerSolid",0.5*container_side,0.5*container_side,1*mm);
  G4LogicalVolume* containerLog = new G4LogicalVolume(containerSolid,lead_mat,"containerLog");
  new G4PVPlacement(layer_rotation,
                    container_pos,
                    containerLog,
                    "containerPhy",
                    fLogicalWorld,
                    false,
                    0,
                    true);


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

  G4double sourceholderThickness = 13*mm, sourceholderRmin = 0, sourceholderRmax = 6*mm;
  G4double holderinsideThickness = 0.006*mm, holderinsideRmin = 0, holderinsideRmax = 4*mm;
  G4double sToOL = 9.65*mm;

  G4Tubs* motherHolder = new G4Tubs("motherHolder", sourceholderRmin,
                                    sourceholderRmax, 0.5*sourceholderThickness, 0, 2*pi);
  G4Tubs* subHolder = new G4Tubs("subHolder", holderinsideRmin, holderinsideRmax, 0.5*holderinsideThickness,
                                  0, 2*pi);

  G4LogicalVolume* sourceMatrixLog = new G4LogicalVolume(subHolder, source_mat, "sourceMatrixLog");
  G4ThreeVector sourceMatrixPos = G4ThreeVector(0,0,-1*sToOL-0.5*holderinsideThickness);
  new G4PVPlacement(nullptr,
                    sourceMatrixPos,
                    sourceMatrixLog,
                    "sourceMatrixPhy",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

  G4ThreeVector ztrans = G4ThreeVector(0,0,0.5*sourceholderThickness-0.5*holderinsideThickness);
  G4RotationMatrix* rots = new G4RotationMatrix;
  rots->rotateX(0);
  G4SubtractionSolid* sourceholderSolid = new G4SubtractionSolid("sourceholderSolid",
                                                             motherHolder, subHolder, rots,ztrans);

  G4LogicalVolume* sourceholderLog = new G4LogicalVolume(sourceholderSolid, sourceholder_mat, "sourceholderLog");

  G4ThreeVector sourceholderPos = G4ThreeVector(0,0,-1*sToOL-0.5*sourceholderThickness);

  new G4PVPlacement(nullptr,
                    sourceholderPos,
                    sourceholderLog,
                    "sourceholderPhy",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

  //G4Material* BMSdet_mat = G4Material::GetMaterial("Silicon");
  G4Material* BMSdet_mat = G4Material::GetMaterial("LaBr");
  

  G4double BMSdetThickness = 1*mm, BMSdetringRmin = 6*mm, BMSdetringRmax = 13.2*mm;
  //G4Tubs* BMSdetsolid = new G4Tubs("BMSdetsolid", BMSdetringRmin, BMSdetringRmax, 0.5*BMSdetThickness, 0, 2*pi);

  //G4ThreeVector BMSdetposition = G4ThreeVector(0,0,0.5*BMSdetThickness);
  
///////////////////////////////////////////////////
//////// the covering mode BMS detector ///////////
///////////////////////////////////////////////////


   G4double BMStrdxymin = 12*mm+2*BMSdetThickness, BMStrdxymax = 40*mm+2*BMSdetThickness, BMStrdy = 10*mm;
   G4double BMStrdxymin2 = BMStrdxymin-2*BMSdetThickness, BMStrdxymax2 = BMStrdxymax-2*BMSdetThickness;
  
   G4Trd* BMSdetMothersolid = new G4Trd("BMSdetMothersolid",
                                        0.5*BMStrdxymin,0.5*BMStrdxymax,
                                        0.5*BMStrdxymin,0.5*BMStrdxymax,0.5*BMStrdy);
                                    
   G4Trd* BMSdetDausolid = new G4Trd("BMSdetMothersolid",
                                        0.5*BMStrdxymin2,0.5*BMStrdxymax2,
                                        0.5*BMStrdxymin2,0.5*BMStrdxymax2,0.5*BMStrdy);

   G4SubtractionSolid* BMSdetsolid = new G4SubtractionSolid("BMSdetsolid",
                                                              BMSdetMothersolid, BMSdetDausolid, rots,G4ThreeVector(0,0,0));

   G4ThreeVector BMSdetposition = G4ThreeVector(0,0,0.5*BMSdetThickness+0.5*BMStrdy);
/////////////////////////////////////////////////////////////


  BMSdetlog = new G4LogicalVolume(BMSdetsolid, BMSdet_mat, "BMSdetlog");

  
  new G4PVPlacement(nullptr,
                    BMSdetposition,
                    BMSdetlog,
                    "BMSdetphy",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);

//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

  //G4double sample_xy = 100*mm, sample_z = 1*mm;
  //G4Box* sampleSolid = new G4Box("sampleSolid", 0.5*sample_xy, 0.5*sample_xy, 0.5*sample_z);

  G4double sample_r = 10*mm, sample_z = 0.0021*mm;

  G4double sample_dis = 10*mm;

  G4Tubs* sampleSolid = new G4Tubs("sampleSolid",0,sample_r,0.5*sample_z,0,2*pi);
  G4ThreeVector samplePos = G4ThreeVector(0, 0, 0.5*BMSdetThickness+0.5*sample_z+sample_dis);

  G4LogicalVolume* sampleLog = new G4LogicalVolume(sampleSolid, sample_mat, "sampleLog");
  new G4PVPlacement(nullptr,
                    samplePos,
                    sampleLog,
                    "samplePhy",
                    fLogicalWorld,
                    false,
                    0,
                    checkOverlaps);



/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

  G4Material* TMSdet_mat = G4Material::GetMaterial("LaBr");

  G4double tranDetThickness = 2*mm, tranDetXY = 2*cm;
  G4Box* tranDetSolid = new G4Box("tranDetSolid", 0.5*tranDetXY, 0.5*tranDetXY, 0.5*tranDetThickness);
  tranDetlog = new G4LogicalVolume(tranDetSolid, TMSdet_mat, "tranDetlog");

  G4ThreeVector tranDetPosition = G4ThreeVector(0, 0, 0.5*BMSdetThickness+0.5*sample_z+sample_dis+0.5*tranDetThickness+2*cm);
  new G4PVPlacement(
    nullptr,
    tranDetPosition,
    tranDetlog,
    "transDetPhy",
    fLogicalWorld,
    false,
    0,
    checkOverlaps
  );


/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

G4double slitL = 0.5*mm, shd1Thickness = 2*mm,
        shd1Rmax = sourceholderRmax+shd1Thickness+slitL, shd1Rmin = sourceholderRmax+slitL, 
        shd1Height = sToOL+0.5*sourceholderThickness;

G4Tubs* shd1Solid = new G4Tubs("shd1Solid", shd1Rmin, shd1Rmax, 0.5*shd1Height, 0, 2*pi);
G4LogicalVolume* shd1Log = new G4LogicalVolume(shd1Solid, lead_mat, "shd1Log");

G4ThreeVector shd1Pos = G4ThreeVector(0,0,-0.5*shd1Height);
new G4PVPlacement(nullptr,
                  shd1Pos,
                  shd1Log,
                  "shd1Phy",
                  fLogicalWorld,
                  false,
                  0,
                  checkOverlaps);
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
G4double adjL = 0.5*tan(22.5/180*pi), 
         shd2Rmax = shd1Rmin, shd2Rmin = holderinsideRmax-adjL, shd2Height = sToOL-slitL;

G4Tubs* shd2Solid = new G4Tubs("shd2Solid", shd2Rmin, shd2Rmax, 0.5*shd2Height, 0, 2*pi);
G4LogicalVolume* shd2Log = new G4LogicalVolume(shd2Solid, lead_mat, "shd2Log");

G4ThreeVector shd2Pos = G4ThreeVector(0,0,-0.5*shd2Height);
new G4PVPlacement(nullptr,
                  shd2Pos,
                  shd2Log,
                  "shd2Phy",
                  fLogicalWorld,
                  false,
                  0,
                  checkOverlaps);
///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

G4double shd3Rmin = shd2Rmin, shd3Rmax = sourceholderRmax, shd3Height = 1*mm;

G4Tubs* shd3Solid = new G4Tubs("shd3Solid", shd3Rmin, shd3Rmax, 0.5*shd3Height, 0, 2*pi);
G4LogicalVolume* shd3Log = new G4LogicalVolume(shd3Solid, lead_mat, "shd3Log");

G4ThreeVector shd3Pos = G4ThreeVector(0,0,0.5*shd3Height);
new G4PVPlacement(nullptr,
                  shd3Pos,
                  shd3Log,
                  "shd3Phy",
                  fLogicalWorld,
                  false,
                  0,
                  checkOverlaps);


//   fLayer1Log->SetVisAttributes(G4VisAttributes::GetInvisible());


  //
  //always return the physical World
  //
  return physWorld;
}


void BMSDetectorConstruction::ConstructSDandField(){

    auto sdManger = G4SDManager::GetSDMpointer();
    G4String SDname;

    auto Layer1 = new BMSLayer1SD(SDname = "/Layer1");
    sdManger->AddNewDetector(Layer1);
    BMSdetlog->SetSensitiveDetector(Layer1);

    auto layer2 = new BMSLayer1SD(SDname = "/Layer2");
    sdManger->AddNewDetector(layer2);
    tranDetlog->SetSensitiveDetector(layer2);


}

void BMSDetectorConstruction::ConstructMaterial(){
    // Get NistManager
    G4NistManager* man = G4NistManager::Instance();

    man->FindOrBuildMaterial("G4_AIR");

    // GAGG material
    G4Element*  O = man->FindOrBuildElement( "O");
    G4Element*  S = man->FindOrBuildElement( "S");
    G4Element* Si = man->FindOrBuildElement("Si");
    G4Element* Lu = man->FindOrBuildElement("Lu");
    G4Element* Gd = man->FindOrBuildElement("Gd");
    G4Element* Ca = man->FindOrBuildElement("Ca");
    G4Element* Fe = man->FindOrBuildElement("Fe");
    G4Element* Al = man->FindOrBuildElement("Al");
    G4Element* Br = man->FindOrBuildElement("Br");
    G4Element* La = man->FindOrBuildElement("La");
    G4Element* Ba = man->FindOrBuildElement("Ba");
    G4Element* Rh = man->FindOrBuildElement("Rh");
    G4Element* Pb = man->FindOrBuildElement("Pb");
    G4Element* Y  = man->FindOrBuildElement("Y");
    G4Element* Sn  = man->FindOrBuildElement("Sn");
    G4Element* Sb  = man->FindOrBuildElement("Sb");

    G4Material* LaBr = new G4Material("LaBr",5.06*g/cm3,2);
    LaBr->AddElement(La,1);
    LaBr->AddElement(Br,3);

    G4Material* LSO = new G4Material("Lu2SiO5",7.1*g/cm3,3);
    LSO->AddElement(Lu,2);
    LSO->AddElement(Si,1);
    LSO->AddElement(O,5);

    G4Material* YSO = new G4Material("Y2SiO5",4.5*g/cm3,3);
    YSO->AddElement(Y,2);
    YSO->AddElement(Si,1);
    YSO->AddElement(O,5);

    G4Material* GAGG = new G4Material("Gd3CaAl4O12", 6.63*g/cm3,4);
    GAGG->AddElement(Gd,3);
    GAGG->AddElement(Ca,1);
    GAGG->AddElement(Al,4);
    GAGG->AddElement(O,12);

    G4Material* BaSO4 = new G4Material("BaSO4", 4.49*g/cm3,3);
    BaSO4->AddElement(Ba,1);
    BaSO4->AddElement(S,1);
    BaSO4->AddElement(O,4);

    G4Material* Silicon = new G4Material("Silicon", 2.33*g/cm3,1);
    Silicon->AddElement(Si,1);

    G4Material* iron = new G4Material("iron", 7.874*g/cm3,1);
    iron->AddElement(Fe,1);

    G4Material* SnSb = new G4Material("SnSb", 4.82*g/cm3,2);
    SnSb->AddElement(Sn,1);
    SnSb->AddElement(Sb,1);

    man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    man->FindOrBuildMaterial("G4_LANTHANUM_OXYSULFIDE");
    man->FindOrBuildMaterial("G4_Fe");
    man->FindOrBuildMaterial("G4_FERROUS_SULFATE");
    man->FindOrBuildMaterial("G4_FERRIC_OXIDE");

    G4Material* lead = new G4Material("lead", 11.35*g/cm3,1);
    lead->AddElement(Pb,1);

    G4Material* sourceMatrix = new G4Material("sourceMatrix", 12.41*g/cm3,1);
    sourceMatrix->AddElement(Rh,1);

    G4Material* Aluminium = man->FindOrBuildMaterial("G4_Al"),
              * Zrinium = man->FindOrBuildMaterial("G4_Zr"),
              * Monium = man->FindOrBuildMaterial("G4_Mo"),
              * Vanium = man->FindOrBuildMaterial("G4_V"),
              * Cronium = man->FindOrBuildMaterial("G4_Cr"),
              * Titanium = man->FindOrBuildMaterial("G4_Ti");

    G4Material* Ti_alloy = new G4Material("Ti_alloy", 4.429*g/cm3, 6);
    G4double al_frac = 0.06, zr_frac = 0.02, mo_frac = 0.01, v_frac = 0.01, cr_frac = 0.02,
             ti_frac = 1-al_frac-zr_frac-mo_frac-v_frac-cr_frac;
    Ti_alloy->AddMaterial(Aluminium, al_frac);
    Ti_alloy->AddMaterial(Zrinium, zr_frac);
    Ti_alloy->AddMaterial(Monium, mo_frac);
    Ti_alloy->AddMaterial(Vanium, v_frac);
    Ti_alloy->AddMaterial(Cronium, cr_frac);
    Ti_alloy->AddMaterial(Titanium, ti_frac);

    G4Material* test_mat = new G4Material("test_mat", 5*g/cm3, 2);
    G4double tiAlloy_frac = 0.7, iron_frac = 0.3;
    test_mat->AddMaterial(Ti_alloy, tiAlloy_frac);
    test_mat->AddMaterial(iron, iron_frac);

}




void BMSDetectorConstruction::DefineCommands(){

    fMessenger = new G4GenericMessenger(this,
                                        "/BMSCamera/detector/",
                                        "Detector control");

    //  command setting the thickness of Layer1
    // auto& Layer1ThicknessCmd
    //   = fMessenger->DeclareMethodWithUnit("Layer1Thickness","mm",
    //                               &BMSDetectorConstruction::SetLayer1Thickness,
    //                               "Set thickness of the Layer1.");
    // Layer1ThicknessCmd.SetParameterName("Thickness", true);
    // Layer1ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    // Layer1ThicknessCmd.SetDefaultValue("8.");

    // // command setting the thickness of layer2
    // auto& Layer2ThicknessCmd
    //   = fMessenger->DeclareMethodWithUnit("Layer2Thickness","mm",
    //                               &BMSDetectorConstruction::SetLayer2Thickness,
    //                               "Set thickness of the Layer2.");
    // Layer2ThicknessCmd.SetParameterName("Thickness", true);
    // Layer2ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    // Layer2ThicknessCmd.SetDefaultValue("8.");

    // // command setting the thickness of layer3
    // auto& Layer3ThicknessCmd
    //   = fMessenger->DeclareMethodWithUnit("Layer3Thickness","mm",
    //                               &BMSDetectorConstruction::SetLayer3Thickness,
    //                               "Set thickness of the Layer3.");
    // Layer3ThicknessCmd.SetParameterName("Thickness", true);
    // Layer3ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    // Layer3ThicknessCmd.SetDefaultValue("8.");

    // auto& Cell1SizeCmd
    //   = fMessenger->DeclareMethodWithUnit("Cell1Size","mm",
    //                               &BMSDetectorConstruction::SetCell1Size,
    //                               "Set size of the Cell1.");
    // Cell1SizeCmd.SetParameterName("cellsize", true);
    // Cell1SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    // Cell1SizeCmd.SetDefaultValue("2.");

    // auto& Cell2SizeCmd
    //   = fMessenger->DeclareMethodWithUnit("Cell2Size","mm",
    //                               &BMSDetectorConstruction::SetCell2Size,
    //                               "Set size of the Cell2.");
    // Cell2SizeCmd.SetParameterName("cellsize", true);
    // Cell2SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    // Cell2SizeCmd.SetDefaultValue("2.");

    // auto& Cell3SizeCmd
    //   = fMessenger->DeclareMethodWithUnit("Cell3Size","mm",
    //                               &BMSDetectorConstruction::SetCell3Size,
    //                               "Set size of the Cell3.");
    // Cell3SizeCmd.SetParameterName("cellsize", true);
    // Cell3SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    // Cell3SizeCmd.SetDefaultValue("2.");

    // auto& Layer1ZPositionCmd
    //   = fMessenger->DeclareMethodWithUnit("Layer1Z","mm",
    //                               &BMSDetectorConstruction::SetLayer1ZPosition,
    //                               "Set z position of the Layer1.");
    // Layer1ZPositionCmd.SetParameterName("ZPosition", true);
    // Layer1ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
    // Layer1ZPositionCmd.SetDefaultValue("-18");

    // auto& Layer2ZPositionCmd
    //   = fMessenger->DeclareMethodWithUnit("Layer2Z","mm",
    //                               &BMSDetectorConstruction::SetLayer2ZPosition,
    //                               "Set z position of the Layer2.");
    // Layer2ZPositionCmd.SetParameterName("ZPosition", true);
    // Layer2ZPositionCmd.SetRange("ZPosition>=-50. && ZPosition<50.");
    // Layer2ZPositionCmd.SetDefaultValue("0");

    // auto& Layer3ZPositionCmd
    //   = fMessenger->DeclareMethodWithUnit("Layer3Z","mm",
    //                               &BMSDetectorConstruction::SetLayer3ZPosition,
    //                               "Set z position of the Layer3.");
    // Layer3ZPositionCmd.SetParameterName("ZPosition", true);
    // Layer3ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
    // Layer3ZPositionCmd.SetDefaultValue("18");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
