#include "GPD3D_RootWriter.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

GPD3D_RootWriter::GPD3D_RootWriter()
: m_file(nullptr),
  m_events(nullptr),
  m_hits(nullptr),
  m_runInfo(nullptr),
  m_runID(-1),
  m_filename(""),
  m_runInfoFilled(false)
{
  ResetBuffers_();
}

GPD3D_RootWriter::~GPD3D_RootWriter() {
  Close();
}

void GPD3D_RootWriter::ResetBuffers_() {
  // RunInfo
  b_run_runID = 0;
  b_run_nBeamOnRequested = 0;
  b_run_controlVerbose = 0;
  b_run_runVerbose = 0;
  b_run_eventVerbose = 0;
  b_run_trackingVerbose = 0;
  b_run_nThreads = 1;

  b_run_bgMode = false;
  b_run_bgWriteFITS = false;
  b_run_bgWriteROOT = false;
  b_run_bgRootFile.clear();

  b_run_genUseSphere = false;
  b_run_genSphereRadius_mm = 0.0;
  b_run_genSphereCenterX_mm = 0.0;
  b_run_genSphereCenterY_mm = 0.0;
  b_run_genSphereCenterZ_mm = 0.0;

  b_run_gpsParticle.clear();
  b_run_gpsEnergy_keV = 0.0;

  // Events
  b_evt_runID = 0;
  b_evt_eventID = 0;
  b_evt_nBeamOn = 0;
  b_evt_primaryPDG = 0;
  b_evt_primaryE0_keV = 0.0;
  b_evt_vtxX_mm = b_evt_vtxY_mm = b_evt_vtxZ_mm = 0.0;
  b_evt_dirX = b_evt_dirY = b_evt_dirZ = 0.0;
  b_evt_edepGasTotal_keV = 0.0;
  b_evt_edepGasIon_keV = 0.0;
  b_evt_edepGasNonIon_keV = 0.0;
  b_evt_nHitsEdep = 0;
  b_evt_enteredGas = b_evt_exitedGas = b_evt_stoppedInGas = b_evt_containedGas = false;
  b_GAGGHit = 0;

  // Hits
  b_hit_runID = b_hit_eventID = 0;
  b_hit_trackID = b_hit_parentID = b_hit_stepNum = 0;
  b_hit_pdg = 0;
  b_hit_x_mm = b_hit_y_mm = b_hit_z_mm = 0.0;
  b_hit_edep_keV = b_hit_edepIon_keV = b_hit_edepNonIon_keV = 0.0;
  b_hit_energy_keV = 0.0;
  b_hit_globalTime_ns = 0.0;
  b_hit_trkLength_mm = 0.0;
  b_hit_firstStep = b_hit_lastStep = false;
  b_hit_dirX = b_hit_dirY = b_hit_dirZ = 0.0;
  b_hit_vtxX_mm = b_hit_vtxY_mm = b_hit_vtxZ_mm = 0.0;
  b_hit_vtxDirX = b_hit_vtxDirY = b_hit_vtxDirZ = 0.0;
  b_hit_vtxEnergy_keV = 0.0;
  b_hit_creatorProcType = b_hit_creatorProcSubType = 0;
  b_hit_stepProcType = b_hit_stepProcSubType = 0;
}

bool GPD3D_RootWriter::Open(const std::string& filename, int runID) {
  Close();

  m_filename = filename;
  m_runID = runID;

  m_file = TFile::Open(filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie()) {
    std::cerr << "GPD3D_RootWriter::Open: cannot open file " << filename << "\n";
    m_file = nullptr;
    return false;
  }

  m_runInfoFilled = false;
  BookTrees_();
  return true;
}

void GPD3D_RootWriter::BookTrees_() {
  if (!m_file) return;

  // RunInfo
  m_runInfo = new TTree("RunInfo", "Run-level configuration and counters");
  m_runInfo->Branch("runID", &b_run_runID, "runID/I");
  m_runInfo->Branch("nBeamOnRequested", &b_run_nBeamOnRequested, "nBeamOnRequested/I");
  m_runInfo->Branch("controlVerbose", &b_run_controlVerbose, "controlVerbose/I");
  m_runInfo->Branch("runVerbose", &b_run_runVerbose, "runVerbose/I");
  m_runInfo->Branch("eventVerbose", &b_run_eventVerbose, "eventVerbose/I");
  m_runInfo->Branch("trackingVerbose", &b_run_trackingVerbose, "trackingVerbose/I");
  m_runInfo->Branch("nThreads", &b_run_nThreads, "nThreads/I");

  m_runInfo->Branch("bgMode", &b_run_bgMode, "bgMode/O");
  m_runInfo->Branch("bgWriteFITS", &b_run_bgWriteFITS, "bgWriteFITS/O");
  m_runInfo->Branch("bgWriteROOT", &b_run_bgWriteROOT, "bgWriteROOT/O");
  m_runInfo->Branch("bgRootFile", &b_run_bgRootFile);

  m_runInfo->Branch("genUseSphere", &b_run_genUseSphere, "genUseSphere/O");
  m_runInfo->Branch("genSphereRadius_mm", &b_run_genSphereRadius_mm, "genSphereRadius_mm/D");
  m_runInfo->Branch("genSphereCenterX_mm", &b_run_genSphereCenterX_mm, "genSphereCenterX_mm/D");
  m_runInfo->Branch("genSphereCenterY_mm", &b_run_genSphereCenterY_mm, "genSphereCenterY_mm/D");
  m_runInfo->Branch("genSphereCenterZ_mm", &b_run_genSphereCenterZ_mm, "genSphereCenterZ_mm/D");

  m_runInfo->Branch("gpsParticle", &b_run_gpsParticle);
  m_runInfo->Branch("gpsEnergy_keV", &b_run_gpsEnergy_keV, "gpsEnergy_keV/D");

  // Events
  m_events = new TTree("Events", "Event-level summary (background mode)");
  m_events->Branch("runID", &b_evt_runID, "runID/I");
  m_events->Branch("eventID", &b_evt_eventID, "eventID/I");
  m_events->Branch("nBeamOn", &b_evt_nBeamOn, "nBeamOn/I");

  m_events->Branch("primaryPDG", &b_evt_primaryPDG, "primaryPDG/I");
  m_events->Branch("primaryE0_keV", &b_evt_primaryE0_keV, "primaryE0_keV/D");

  m_events->Branch("vtxX_mm", &b_evt_vtxX_mm, "vtxX_mm/D");
  m_events->Branch("vtxY_mm", &b_evt_vtxY_mm, "vtxY_mm/D");
  m_events->Branch("vtxZ_mm", &b_evt_vtxZ_mm, "vtxZ_mm/D");

  m_events->Branch("dirX", &b_evt_dirX, "dirX/D");
  m_events->Branch("dirY", &b_evt_dirY, "dirY/D");
  m_events->Branch("dirZ", &b_evt_dirZ, "dirZ/D");

  m_events->Branch("edepGasTotal_keV", &b_evt_edepGasTotal_keV, "edepGasTotal_keV/D");
  m_events->Branch("edepGasIon_keV", &b_evt_edepGasIon_keV, "edepGasIon_keV/D");
  m_events->Branch("edepGasNonIon_keV", &b_evt_edepGasNonIon_keV, "edepGasNonIon_keV/D");
  m_events->Branch("nHitsEdep", &b_evt_nHitsEdep, "nHitsEdep/I");

  m_events->Branch("enteredGas", &b_evt_enteredGas, "enteredGas/O");
  m_events->Branch("exitedGas", &b_evt_exitedGas, "exitedGas/O");
  m_events->Branch("stoppedInGas", &b_evt_stoppedInGas, "stoppedInGas/O");
  m_events->Branch("containedGas", &b_evt_containedGas, "containedGas/O");
  m_events->Branch("GAGGHit", &b_GAGGHit, "b_GAGGHit/I");

  // Hits
  m_hits = new TTree("Hits", "Hit-level details (Edep>0)");
  m_hits->Branch("runID", &b_hit_runID, "runID/I");
  m_hits->Branch("eventID", &b_hit_eventID, "eventID/I");
  m_hits->Branch("trackID", &b_hit_trackID, "trackID/I");
  m_hits->Branch("parentID", &b_hit_parentID, "parentID/I");
  m_hits->Branch("stepNum", &b_hit_stepNum, "stepNum/I");
  m_hits->Branch("pdg", &b_hit_pdg, "pdg/I");

  m_hits->Branch("x_mm", &b_hit_x_mm, "x_mm/D");
  m_hits->Branch("y_mm", &b_hit_y_mm, "y_mm/D");
  m_hits->Branch("z_mm", &b_hit_z_mm, "z_mm/D");

  m_hits->Branch("edep_keV", &b_hit_edep_keV, "edep_keV/D");
  m_hits->Branch("edepIon_keV", &b_hit_edepIon_keV, "edepIon_keV/D");
  m_hits->Branch("edepNonIon_keV", &b_hit_edepNonIon_keV, "edepNonIon_keV/D");

  m_hits->Branch("energy_keV", &b_hit_energy_keV, "energy_keV/D");
  m_hits->Branch("globalTime_ns", &b_hit_globalTime_ns, "globalTime_ns/D");
  m_hits->Branch("trkLength_mm", &b_hit_trkLength_mm, "trkLength_mm/D");

  m_hits->Branch("firstStep", &b_hit_firstStep, "firstStep/O");
  m_hits->Branch("lastStep", &b_hit_lastStep, "lastStep/O");

  m_hits->Branch("dirX", &b_hit_dirX, "dirX/D");
  m_hits->Branch("dirY", &b_hit_dirY, "dirY/D");
  m_hits->Branch("dirZ", &b_hit_dirZ, "dirZ/D");

  m_hits->Branch("vtxX_mm", &b_hit_vtxX_mm, "vtxX_mm/D");
  m_hits->Branch("vtxY_mm", &b_hit_vtxY_mm, "vtxY_mm/D");
  m_hits->Branch("vtxZ_mm", &b_hit_vtxZ_mm, "vtxZ_mm/D");

  m_hits->Branch("vtxDirX", &b_hit_vtxDirX, "vtxDirX/D");
  m_hits->Branch("vtxDirY", &b_hit_vtxDirY, "vtxDirY/D");
  m_hits->Branch("vtxDirZ", &b_hit_vtxDirZ, "vtxDirZ/D");

  m_hits->Branch("vtxEnergy_keV", &b_hit_vtxEnergy_keV, "vtxEnergy_keV/D");

  m_hits->Branch("creatorProcType", &b_hit_creatorProcType, "creatorProcType/I");
  m_hits->Branch("creatorProcSubType", &b_hit_creatorProcSubType, "creatorProcSubType/I");
  m_hits->Branch("stepProcType", &b_hit_stepProcType, "stepProcType/I");
  m_hits->Branch("stepProcSubType", &b_hit_stepProcSubType, "stepProcSubType/I");
}

void GPD3D_RootWriter::Close() {
  if (!m_file) return;

  m_file->cd();
  if (m_runInfo) m_runInfo->Write("", TObject::kOverwrite);
  if (m_events)  m_events->Write("", TObject::kOverwrite);
  if (m_hits)    m_hits->Write("", TObject::kOverwrite);

  m_file->Close();

  delete m_file;
  m_file = nullptr;

  m_runInfo = nullptr;
  m_events  = nullptr;
  m_hits    = nullptr;

  m_runID = -1;
  m_filename.clear();
  m_runInfoFilled = false;
}

void GPD3D_RootWriter::FillRunInfo(
  int runID,
  int nBeamOnRequested,
  int controlVerbose,
  int runVerbose,
  int eventVerbose,
  int trackingVerbose,
  int nThreads,
  bool bgMode,
  bool bgWriteFITS,
  bool bgWriteROOT,
  const std::string& bgRootFile,
  bool genUseSphere,
  double genSphereRadius_mm,
  double genSphereCenterX_mm,
  double genSphereCenterY_mm,
  double genSphereCenterZ_mm,
  const std::string& gpsParticle,
  double gpsEnergy_keV
) {
  if (!m_file || !m_runInfo) return;

  // Fill only once per file
  if (m_runInfoFilled) return;

  b_run_runID = runID;
  b_run_nBeamOnRequested = nBeamOnRequested;
  b_run_controlVerbose = controlVerbose;
  b_run_runVerbose = runVerbose;
  b_run_eventVerbose = eventVerbose;
  b_run_trackingVerbose = trackingVerbose;
  b_run_nThreads = nThreads;

  b_run_bgMode = bgMode;
  b_run_bgWriteFITS = bgWriteFITS;
  b_run_bgWriteROOT = bgWriteROOT;
  b_run_bgRootFile = bgRootFile;

  b_run_genUseSphere = genUseSphere;
  b_run_genSphereRadius_mm = genSphereRadius_mm;
  b_run_genSphereCenterX_mm = genSphereCenterX_mm;
  b_run_genSphereCenterY_mm = genSphereCenterY_mm;
  b_run_genSphereCenterZ_mm = genSphereCenterZ_mm;

  b_run_gpsParticle = gpsParticle;
  b_run_gpsEnergy_keV = gpsEnergy_keV;

  m_runInfo->Fill();
  m_runInfoFilled = true;
}

void GPD3D_RootWriter::FillEventRow(
  int runID,
  int eventID,
  int nBeamOn,
  int primaryPDG,
  double primaryE0_keV,
  double vtxX_mm, double vtxY_mm, double vtxZ_mm,
  double dirX, double dirY, double dirZ,
  double edepGasTotal_keV,
  double edepGasIon_keV,
  double edepGasNonIon_keV,
  int nHitsEdep,
  bool enteredGas, bool exitedGas, bool stoppedInGas, bool containedGas, int GAGGHit
) {
  if (!m_file || !m_events) return;

  b_evt_runID = runID;
  b_evt_eventID = eventID;
  b_evt_nBeamOn = nBeamOn;

  b_evt_primaryPDG = primaryPDG;
  b_evt_primaryE0_keV = primaryE0_keV;

  b_evt_vtxX_mm = vtxX_mm;
  b_evt_vtxY_mm = vtxY_mm;
  b_evt_vtxZ_mm = vtxZ_mm;

  b_evt_dirX = dirX;
  b_evt_dirY = dirY;
  b_evt_dirZ = dirZ;

  b_evt_edepGasTotal_keV = edepGasTotal_keV;
  b_evt_edepGasIon_keV = edepGasIon_keV;
  b_evt_edepGasNonIon_keV = edepGasNonIon_keV;
  b_evt_nHitsEdep = nHitsEdep;

  b_evt_enteredGas = enteredGas;
  b_evt_exitedGas = exitedGas;
  b_evt_stoppedInGas = stoppedInGas;
  b_evt_containedGas = containedGas;
  b_GAGGHit = GAGGHit;

  m_events->Fill();
}

void GPD3D_RootWriter::FillHitRow(
  int runID, int eventID,
  int trackID, int parentID, int stepNum,
  int pdg,
  double x_mm, double y_mm, double z_mm,
  double edep_keV, double edepIon_keV, double edepNonIon_keV,
  double energy_keV,
  double globalTime_ns,
  double trkLength_mm,
  bool firstStep, bool lastStep,
  double dirX, double dirY, double dirZ,
  double vtxX_mm, double vtxY_mm, double vtxZ_mm,
  double vtxDirX, double vtxDirY, double vtxDirZ,
  double vtxEnergy_keV,
  int creatorProcType, int creatorProcSubType,
  int stepProcType, int stepProcSubType
) {
  if (!m_file || !m_hits) return;

  b_hit_runID = runID;
  b_hit_eventID = eventID;

  b_hit_trackID = trackID;
  b_hit_parentID = parentID;
  b_hit_stepNum = stepNum;
  b_hit_pdg = pdg;

  b_hit_x_mm = x_mm;
  b_hit_y_mm = y_mm;
  b_hit_z_mm = z_mm;

  b_hit_edep_keV = edep_keV;
  b_hit_edepIon_keV = edepIon_keV;
  b_hit_edepNonIon_keV = edepNonIon_keV;

  b_hit_energy_keV = energy_keV;
  b_hit_globalTime_ns = globalTime_ns;
  b_hit_trkLength_mm = trkLength_mm;

  b_hit_firstStep = firstStep;
  b_hit_lastStep = lastStep;

  b_hit_dirX = dirX;
  b_hit_dirY = dirY;
  b_hit_dirZ = dirZ;

  b_hit_vtxX_mm = vtxX_mm;
  b_hit_vtxY_mm = vtxY_mm;
  b_hit_vtxZ_mm = vtxZ_mm;

  b_hit_vtxDirX = vtxDirX;
  b_hit_vtxDirY = vtxDirY;
  b_hit_vtxDirZ = vtxDirZ;

  b_hit_vtxEnergy_keV = vtxEnergy_keV;

  b_hit_creatorProcType = creatorProcType;
  b_hit_creatorProcSubType = creatorProcSubType;
  b_hit_stepProcType = stepProcType;
  b_hit_stepProcSubType = stepProcSubType;

  m_hits->Fill();
}
