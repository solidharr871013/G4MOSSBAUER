// C++ source file to encapsulate the G4MOSS physics process


#include "G4MOSSPhysics.hh"

#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"

G4MOSSPhysics::G4MOSSPhysics(const G4String &name, G4bool use_xsec_tables_in,
  G4bool use_xsec_integration_in, G4bool force_isotropic_in)
  : use_xsec_tables(use_xsec_tables_in),
    use_xsec_integration(use_xsec_integration_in),
    force_isotropic(force_isotropic_in) {
  G4LossTableManager::Instance();
  SetPhysicsType(0);
}

G4MOSSPhysics::~G4MOSSPhysics()
{}

void G4MOSSPhysics::ConstructParticle()
{}

void G4MOSSPhysics::ConstructProcess() {
  G4MOSS *nrf = new G4MOSS("NRF", false, use_xsec_tables=false, use_xsec_integration, force_isotropic);


// Check if Geant4 version is >= 10.3.0
//#if G4MAJV * 10000 + G4MINV * 100 + G4SUBV > 100299
  auto aParticleIterator=GetParticleIterator();
//#endif
  
  aParticleIterator->reset();
  while( (*aParticleIterator)() ) {
    G4ParticleDefinition       *particle = aParticleIterator->value();
    G4ProcessManager *particleProcessMgr = particle->GetProcessManager();
    G4String                particleName = particle->GetParticleName();
    /*
    G4cout << " G4MOSSPhysics: particle = "           << particle << G4endl;
    G4cout << " G4MOSSPhysics: particleProcessMgr = " << particleProcessMgr << G4endl;
    G4cout << " G4MOSSPhysics: particleName = "       << particleName << G4endl;
    G4cout << G4endl;
    */
    if (particleName == "gamma") {
      particleProcessMgr->AddDiscreteProcess(nrf);
    }
  }
}
