#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "GPD3D_PhysicsList.hh"
#include "GPD3D_Physics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"

#include "G4StepLimiterPhysics.hh"
// #include "G4StepLimiterPhysics.hh"
// #include "G4UserSpecialCuts.hh"
/*!
 */
GPD3D_PhysicsList::GPD3D_PhysicsList() :
  G4VModularPhysicsList()
{
  G4int verb = 0;
  SetVerboseLevel(verb);
  // RegisterPhysics(new G4StepLimiterPhysics()); // step cut
  G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
  stepLimitPhys->SetApplyToAll(true); // activates step limit for ALL particles
  RegisterPhysics(stepLimitPhys);

  // RegisterPhysics(new G4UserSpecialCuts()); // track cut
  RegisterPhysics(new GPD3D_Physics(verb));
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4HadronElasticPhysics());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  RegisterPhysics(new G4StoppingPhysics());
  RegisterPhysics(new G4IonPhysics());
}


/*!
 */
GPD3D_PhysicsList::~GPD3D_PhysicsList()
{}


/*!
  Set cuts on secondary particles production
 */
void GPD3D_PhysicsList::SetCuts()
{
  // G4double lowEnergyRange = 15. * eV;
  G4double lowEnergyRange  = 10. * eV;
  G4double highEnergyRange = 200. * GeV; // comfortably above 100 GeV
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowEnergyRange, highEnergyRange);

  static G4double gammaCutRange = 50.  * um;
  // static G4double electronCutRange = 2. * um;
  // static G4double electronCutRange = 0.02 * um;
  static G4double electronCutRange = 0.2 * nm;
  SetCutValue(gammaCutRange, "gamma");
  SetCutValue(electronCutRange, "e-");
}
