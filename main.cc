#include "globals.hh"

#include <vector>
#include <chrono>

// #include "G4RunManagerFactory.hh"
#include "G4String.hh"

#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#include "GPD3D_DetectorConstruction.hh"
#include "GPD3D_PrimaryGeneratorAction.hh"
#include "GPD3D_RunAction.hh"
#include "GPD3D_EventAction.hh"
#include "GPD3D_SteppingAction.hh"
#include "GPD3D_PhysicsList.hh"

#include "Randomize.hh"

using namespace std;

int main(int argc, char** argv)
{
  auto start = std::chrono::high_resolution_clock::now();
  G4cout << "Application starting..." << G4endl;

  vector<G4String> macros;
  G4bool interactive = false;

  // Parse command line arguments
  if (argc == 1) {
    interactive = true; // default: interactive
  } else {
    for (int i = 1; i < argc; i++) {
      G4String arg = argv[i];
      if (arg == "-i" || arg == "--interactive") {
        interactive = true;
        continue;
      } else {
        macros.push_back(arg); // treat as macro file
      }
    }
  }

  // If no macro is given, run the default visualization macro
  if (interactive && macros.empty()) {
    macros.push_back("macros/vis.mac");
  }

  // Choose the Random engine and set seed
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(time(NULL));

  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new GPD3D_PhysicsList());
  auto det = new GPD3D_DetectorConstruction();
  runManager->SetUserInitialization(det);
  runManager->SetUserAction(new GPD3D_PrimaryGeneratorAction());

  if (!macros.empty())
    runManager->SetUserAction(new GPD3D_RunAction(macros.front()));
  else
    runManager->SetUserAction(new GPD3D_RunAction());

  auto ea = new GPD3D_EventAction();
  runManager->SetUserAction(ea);


  runManager->SetUserAction(new GPD3D_SteppingAction(ea, det));
  runManager->Initialize();

  // Vis + UI manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  G4UIExecutive* ui = nullptr;
  if (interactive) {
    // Create Qt UI first
    ui = new G4UIExecutive(argc, argv);
  }

  // Execute macros (default vis.mac included if you decided so)
  for (const auto& macro : macros) {
    UImanager->ApplyCommand("/control/execute " + macro);
  }

  // Start UI session last
  if (ui) {
    ui->SessionStart();
    delete ui;
    ui = nullptr;
  }


  delete visManager;
  delete runManager;

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "++++++++++++++TIME++++++++++++++" << std::endl;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

  return 0;
}
