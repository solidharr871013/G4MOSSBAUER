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
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      File name:     G4MOSSNuclearLevelManager
//
//      Author:        Micheal zhang
//      Creation date: October 2021
//
//      Adapted from   G4NuclearLevelManager
//      CERN, Geneva, Switzerland
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//      Creation date: 24 October 1998
//
//      Modification history (G4NuclearLevelManager):
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions
//
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//        02 May 2003,   Vladimir Ivanchenko remove rublic copy contructor
//
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVELMANAGER_HH
#define G4NUCLEARLEVELMANAGER_HH

#include <fstream>

#include "globals.hh"
#include "G4MOSSPtrLevelVector.hh"
#include "G4MOSSNuclearLevel.hh"
#include "G4ios.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Integrator.hh"
#include "G4Exp.hh"

using std::ofstream;

class G4MOSSNuclearLevelManager {
 public:
  G4MOSSNuclearLevelManager(G4bool Verbose = false);
  G4MOSSNuclearLevelManager(const G4int Z, const G4int A, const G4String& filename,
         G4bool Verbose = false);
  G4MOSSNuclearLevelManager(const G4MOSSNuclearLevelManager & right);

  ~G4MOSSNuclearLevelManager();

  void SetNucleus(const G4int Z, const G4int A, const G4String& filename, G4bool standalone = false);

  G4bool IsValid() const;

  G4int NumberOfLevels() const;

  const G4MOSSPtrLevelVector* GetLevels() const;

  G4MOSSNuclearLevel* NearestLevel(const G4double energy, const G4double eDiffMax = 9999.*GeV);

  G4MOSSNuclearLevel* NearestLevelRecoilAbsorb(const G4double GammaEnergy, const G4double eDiffMax = 9999.*GeV);

  G4MOSSNuclearLevel* NearestLevelRecoilEmit(const G4double LevelEnergy,
               const G4double GammaEnergy, const G4double eDiffMax = 9999.*GeV);

  const G4MOSSNuclearLevel* LowestLevel()  const;
  const G4MOSSNuclearLevel* HighestLevel() const;

  G4double MinLevelEnergy() const;
  G4double MaxLevelEnergy() const;

  G4double GetGroundStateSpin() const; 
  G4double GetGroundStateParity() const;

  G4double GetTotalCC() const; //added in feb. 10th, 2022

  void PrintAll();

  void PrintAllTabular(ofstream& file);

  void PrintLevelEnergies();

  void ReadTDebyeData(G4bool standalone);
  void SetTDebye(G4double TDebye);
  G4double GetTDebye() const;

  G4double TeffIntegrand(G4double) const;
  G4double CalcTeff(G4double, G4double = 300*kelvin);
  G4double GetTeff();
  void SetTeff(G4double Teff);

 private:
  const G4MOSSNuclearLevelManager& operator=(const G4MOSSNuclearLevelManager &right);
  G4bool operator==(const G4MOSSNuclearLevelManager &right) const;
  G4bool operator!=(const G4MOSSNuclearLevelManager &right) const;

  G4bool Read(std::ifstream& aDataFile);

  void ReadGroundStateProperties(G4bool standalone = false);

  void MakeLevels();

  void delete_bad_levels();

  void print_to_standalone();

  G4int _nucleusA;
  G4int _nucleusZ;
  G4String _fileName;
  G4bool _validity;
  G4MOSSPtrLevelVector* _levels;

  G4double _gsAngularMomentum;
  G4double _gsParity;

  G4double _levelEnergy;
  G4double _gammaEnergy;
  G4double _probability;
  G4double _polarity;
  G4double _halfLife;
  G4double _angularMomentum;
  G4double _kCC;
  G4double _l1CC;
  G4double _l2CC;
  G4double _l3CC;
  G4double _m1CC;
  G4double _m2CC;
  G4double _m3CC;
  G4double _m4CC;
  G4double _m5CC;
  G4double _nPlusCC;
  G4double _totalCC;

  G4bool   _Verbose;

  G4double _TDebye;
  G4double _Teff;
};

#endif
