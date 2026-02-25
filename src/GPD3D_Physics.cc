#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
// #include "G4VEmModel.hh"
// #include "G4ElementData.hh"
// #include "G4LPhysicsFreeVector.hh"

// Particles
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Alpha.hh"
#include "G4He3.hh"
#include "G4GenericIon.hh"
#include "G4Neutron.hh"

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePolarizedPhotoElectricGDModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermorePolarizedRayleighModel.hh"

// e+-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eplusAnnihilation.hh"

// mu+-
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

// hadrons
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4alphaIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// Multiple scattering models
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

// Interfaces
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

// Factory
#include "G4PhysicsConstructorFactory.hh"

#include "GPD3D_Physics.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(GPD3D_Physics);

/*!
 */
GPD3D_Physics::GPD3D_Physics(G4int ver, const G4String &) :
  G4VPhysicsConstructor("GPD3D_Physics"),
  verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  //REMEMBER
  param->SetMinEnergy(10*eV);
  // param->SetMinEnergy(0.1* keV);
  param->SetMaxEnergy(100*GeV);
  param->SetLowestElectronEnergy(50*eV);
  // Set a greater number of bin for model evaluation
  param->SetNumberOfBinsPerDecade(40);
  param->ActivateAngularGeneratorForIonisation(true);
  // Activate fluorescence and Auger emission, ignoring the production cuts
  param->SetFluo(true);
  param->SetAuger(true);
  param->SetDeexcitationIgnoreCut(true);
  SetPhysicsType(bElectromagnetic);
}


/*!
 */
GPD3D_Physics::~GPD3D_Physics()
{}


/*!
  @todo we need to define also muons.
 */
void GPD3D_Physics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
  G4Proton::Proton();
  G4AntiProton::AntiProton();
  G4Alpha::Alpha();
  G4He3::He3();
  G4GenericIon::GenericIonDefinition();
  G4Neutron::NeutronDefinition();
}


/*!
  We changed the ionization model in order to use the G4MollerBhabhaModel also
  for lower energies, because it seems to give a more accurate energy loss in
  gas, see:
  https://indico.fnal.gov/getFile.py/access?contribId=93&sessionId=9&resId=0&materialId=slides&confId=4535

  We tweaked also the step size for the ionization model using the corresponding
  SetStepFunction, see:
  http://hypernews.slac.stanford.edu/HyperNews/geant4/get/emprocess/1092/1/1/1.html

  @todo we need to define also muons and protons processes.
 */
void GPD3D_Physics::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // muon & hadron bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
  G4MuPairProduction* mup = new G4MuPairProduction();
  //G4hBremsstrahlung* pib = new G4hBremsstrahlung();
  //G4hPairProduction* pip = new G4hPairProduction();
  //G4hBremsstrahlung* kb = new G4hBremsstrahlung();
  //G4hPairProduction* kp = new G4hPairProduction();
  G4hBremsstrahlung* pb = new G4hBremsstrahlung();
  G4hPairProduction* pp = new G4hPairProduction();

  // muon & hadron multiple scattering
  G4MuMultipleScattering* mumsc = new G4MuMultipleScattering();
  mumsc->AddEmModel(0, new G4WentzelVIModel());
  //G4MuMultipleScattering* pimsc = new G4MuMultipleScattering();
  //pimsc->AddEmModel(0, new G4WentzelVIModel());
  //G4MuMultipleScattering* kmsc = new G4MuMultipleScattering();
  //kmsc->AddEmModel(0, new G4WentzelVIModel());
  G4MuMultipleScattering* pmsc = new G4MuMultipleScattering();
  pmsc->AddEmModel(0, new G4WentzelVIModel());
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

  // nuclear stopping
  G4NuclearStopping* pnuc = new G4NuclearStopping();

  // Add Livermore EM Processes
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();

  // High energy limit for e+- scattering models
  G4double highEnergyLimit = 100*MeV;

  // Applicability range for Livermore models (for higher energies, the
  // Standard models are used)
  G4double LivermoreHighEnergyLimit = GeV;

  while( (*myParticleIterator)() ){

    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      // ELIMINATE THIS PART DUE TO THE G4LivermorePolarizedPhotoElectricGDModel 
      G4LivermorePolarizedPhotoElectricGDModel* theLivermorePhotoElectricModel = new G4LivermorePolarizedPhotoElectricGDModel();
      theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      ph->RegisterProcess(thePhotoElectricEffect, particle);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermorePolarizedComptonModel* theLivermoreComptonModel = new G4LivermorePolarizedComptonModel();
      theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
      ph->RegisterProcess(theComptonScattering, particle);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermorePolarizedGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermorePolarizedGammaConversionModel();
      theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      ph->RegisterProcess(theGammaConversion, particle);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermorePolarizedRayleighModel* theRayleighModel = new G4LivermorePolarizedRayleighModel();
      theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theRayleigh->AddEmModel(0, theRayleighModel);
      ph->RegisterProcess(theRayleigh, particle);

    }

    else if (particleName == "e-") {
      // Multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4UrbanMscModel* msc1 = new G4UrbanMscModel();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);
      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(ss, particle);

      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      // G4LivermoreIonisationModel* theIoniLivermore = new G4LivermoreIonisationModel();
      // theIoniLivermore->SetHighEnergyLimit(0.1*MeV);
      // eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
      // eIoni->SetStepFunction(0.01, 0.1*um);
      eIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(eIoni, particle);

      // Bremsstrahlung from standard
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      ph->RegisterProcess(eBrem, particle);
    }

    else if (particleName == "e+") {
      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4UrbanMscModel* msc1 = new G4UrbanMscModel();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
      ph->RegisterProcess(ss, particle);

    } else if (particleName == "mu+" || particleName == "mu-") {
      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.2, 50*um);

      ph->RegisterProcess(mumsc, particle);
      ph->RegisterProcess(muIoni, particle);
      ph->RegisterProcess(mub, particle);
      ph->RegisterProcess(mup, particle);
      ph->RegisterProcess(new G4CoulombScattering(), particle);
    }

    else if (particleName == "proton" ||
                particleName == "anti_proton") {
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.2, 50*um);

      ph->RegisterProcess(pmsc, particle);
      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(pb, particle);
      ph->RegisterProcess(pp, particle);
      ph->RegisterProcess(pnuc, particle);
      ph->RegisterProcess(new G4CoulombScattering(), particle);
    }

    else if (particleName == "alpha" ||
                particleName == "He3" ) {

       // Identical to G4EmStandardPhysics_option3
       G4hMultipleScattering* msc = new G4hMultipleScattering();
       G4ionIonisation* ionIoni = new G4ionIonisation();
       ionIoni->SetStepFunction(0.1, 10*um);

       ph->RegisterProcess(msc, particle);
       ph->RegisterProcess(ionIoni, particle);
       ph->RegisterProcess(pnuc, particle);
    }

    else if (particleName == "GenericIon") {

      // Identical to G4EmStandardPhysics_option3
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);
      // ionIoni->SetStepFunction(rc::ionizationStepSlope,
                               // rc::ionizationStepPivot);
      ph->RegisterProcess(hmsc, particle);
      ph->RegisterProcess(ionIoni, particle);
      ph->RegisterProcess(pnuc, particle);
    }
  }
  // Deexcitation
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);

  G4EmModelActivator mact(GetPhysicsName());
}
