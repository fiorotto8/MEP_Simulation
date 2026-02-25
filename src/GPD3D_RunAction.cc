#include "GPD3D_RunAction.hh"

GPD3D_RunAction::GPD3D_RunAction()
: G4UserRunAction()
{
  fName = "result";
  SetAnalysis();
}

GPD3D_RunAction::GPD3D_RunAction(const char *name)
: G4UserRunAction()
{
  fName = name;
  SetAnalysis();
}

GPD3D_RunAction::~GPD3D_RunAction()
{
  
  // G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
  // analysisManager -> Write();
  // analysisManager -> CloseFile();

  delete G4CsvAnalysisManager::Instance();
}

void GPD3D_RunAction::BeginOfRunAction(const G4Run*)
{
}

void GPD3D_RunAction::EndOfRunAction(const G4Run*)
{
}

void GPD3D_RunAction::SetAnalysis()
{
  // G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
  // analysisManager -> OpenFile(fName);

  // analysisManager-> SetFirstNtupleId(1);
  // analysisManager -> CreateNtuple("TEST_CSV", "STEP_ALL");
  // analysisManager -> CreateNtupleIColumn("eventID");
  // analysisManager -> CreateNtupleIColumn("volumeID");
  // analysisManager -> CreateNtupleIColumn("trackID");
  // analysisManager -> CreateNtupleIColumn("parentID");
  // analysisManager -> CreateNtupleDColumn("position_x");
  // analysisManager -> CreateNtupleDColumn("position_y");
  // analysisManager -> CreateNtupleDColumn("position_z");
  // analysisManager -> CreateNtupleDColumn("track_length");
  // analysisManager -> CreateNtupleDColumn("PE_phi");
  // analysisManager -> CreateNtupleDColumn("PE_theta");
  // analysisManager -> CreateNtupleDColumn("TIME");
  // analysisManager -> CreateNtupleDColumn("KineticE");
  // analysisManager -> FinishNtuple();
}
