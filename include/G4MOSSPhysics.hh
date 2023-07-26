#ifndef G4MOSSPhysics_h
#define G4MOSSPhysics_h 1

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsListHelper.hh"
#include "G4LossTableManager.hh"
#include "G4BuilderType.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

#include "G4Gamma.hh"

#include "G4MOSS.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4VModularPhysicsList.hh"


class G4MOSSPhysics : public G4VPhysicsConstructor {
 public:
  G4MOSSPhysics(const G4String &name, G4bool, G4bool, G4bool);

  virtual ~G4MOSSPhysics();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

 private:
  G4bool use_xsec_tables;
  G4bool use_xsec_integration;
  G4bool force_isotropic;
};

#endif
