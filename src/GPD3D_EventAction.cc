#include "GPD3D_EventAction.hh"
#include "GPD3D_EventInfo.hh"
#include "GPD3D_FITSWriter.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include <iostream>
#include <random>
#include <cmath>
#include <tuple>
#include <vector>
#include <string>
#include <chrono>
#include <thread>
#include <limits>
#include <algorithm>

#include "Randomize.hh"
#include <boost/math/special_functions/gamma.hpp>

// -----------------------------
// ctor/dtor
// -----------------------------
GPD3D_EventAction::GPD3D_EventAction()
  : G4UserEventAction()
  , m_gasCellHCID(-1)
  , m_hit(nullptr)
  , m_trackID(-1)
  , m_parentID(-1)
  , m_stepNum(-1)
  , m_eDep(0.)
  , m_energy(0.)
  , m_trkLength(0.)
  , m_pos(G4ThreeVector())
  , m_dir(G4ThreeVector())
  , m_trackProcName("")
  , m_stepProcName("")
  , m_partName("")
  , m_partState(fAlive)
  , m_ionTrack()
  , m_firstStep(false)
  , m_lastStep(false)
  , m_vertPos(G4ThreeVector())
  , m_vertDir(G4ThreeVector())
  , m_vertEnergy(-1.)
  , EVENT_ID_LIST()
  , PHO_ENG_LIST()
  , ABS_POS_X_LIST()
  , ABS_POS_Y_LIST()
  , ABS_POS_Z_LIST()
  , PE_ENG_LIST()
  , PE_PHI_LIST()
  , PE_THETA_LIST()
  , AUG_ENG_LIST()
  , AUG_PHI_LIST()
  , AUG_THETA_LIST()
  , TRK_LEN_LIST()
  , ION_POS_X_LIST()
  , ION_POS_Y_LIST()
  , ION_POS_Z_LIST()
  , MAP_COL_LIST()
  , MAP_ROW_LIST()
  , COUNTING_MAP_3D()
  , MAP_TOA_3D()
  // background/runinfo
  , m_bgMode(false)
  , m_bgWriteFITS(true)
  , m_bgWriteROOT(false)
  , m_bgRootFile("background.root")
  , m_reqNBeamOn(-1)
  , m_controlVerbose(-1)
  , m_runVerbose(-1)
  , m_eventVerbose(-1)
  , m_trackingVerbose(-1)
  , m_reqNThreads(-1)
  , m_genUseSphere(false)
  , m_genSphereRadius_mm(0.0)
  , m_genSphereCenterX_mm(0.0)
  , m_genSphereCenterY_mm(0.0)
  , m_genSphereCenterZ_mm(0.0)
  , m_gpsParticle("")
  , m_gpsEnergy_keV(-1.0)
  , m_bgMessenger(nullptr)
  , m_runInfoMessenger(nullptr)
  , m_bgNgen(0)
  , m_bgNhit(0)
  , m_bgSumEdep(0.)
  , m_bgSumEdepHit(0.)
  , m_rootWriter(nullptr)
  , m_runInfoFilled(false)
  , m_openRunID(-999)
{
  // Background controls
  m_bgMessenger = new G4GenericMessenger(this, "/gpd3d/bg/", "Background-mode controls");
  m_bgMessenger->DeclareProperty("enable", m_bgMode,
      "Enable background scoring mode (skip PE/Auger logic and heavy digitization).");
  m_bgMessenger->DeclareProperty("writeFITS", m_bgWriteFITS,
      "Write output.fits at end of run (polarimetry mode).");

  m_bgMessenger->DeclareProperty("writeROOT", m_bgWriteROOT,
      "Write ROOT output (Events + Hits + RunInfo) in background mode.");
  m_bgMessenger->DeclareProperty("rootFile", m_bgRootFile,
      "ROOT file name for background mode (e.g. bg.root). Supports %RUN% token.");

  // Robust nBeamOn capture (set this in macro next to /run/beamOn)
  m_bgMessenger->DeclareProperty("nBeamOn", m_reqNBeamOn,
      "Store requested /run/beamOn N for writing into ROOT (robust vs MT splitting).");

  // RunInfo mirror (optional, but solves “I want macro params in RunInfo”)
  m_runInfoMessenger = new G4GenericMessenger(this, "/gpd3d/runinfo/", "RunInfo mirror values (set in macro)");
  m_runInfoMessenger->DeclareProperty("controlVerbose",  m_controlVerbose,  "Mirror /control/verbose");
  m_runInfoMessenger->DeclareProperty("runVerbose",      m_runVerbose,      "Mirror /run/verbose");
  m_runInfoMessenger->DeclareProperty("eventVerbose",    m_eventVerbose,    "Mirror /event/verbose");
  m_runInfoMessenger->DeclareProperty("trackingVerbose", m_trackingVerbose, "Mirror /tracking/verbose");
  m_runInfoMessenger->DeclareProperty("nThreads",        m_reqNThreads,     "Mirror /run/numberOfThreads");

  m_runInfoMessenger->DeclareProperty("useSphere", m_genUseSphere, "Mirror /gpd3d/gen/useSphere (if you use it)");
  m_runInfoMessenger->DeclarePropertyWithUnit("sphereRadius", "mm", m_genSphereRadius_mm,
      "Mirror /gpd3d/gen/sphereRadius (mm)");
  m_runInfoMessenger->DeclarePropertyWithUnit("sphereCenterX", "mm", m_genSphereCenterX_mm,
      "Mirror /gpd3d/gen/sphereCenterX (mm)");
  m_runInfoMessenger->DeclarePropertyWithUnit("sphereCenterY", "mm", m_genSphereCenterY_mm,
      "Mirror /gpd3d/gen/sphereCenterY (mm)");
  m_runInfoMessenger->DeclarePropertyWithUnit("sphereCenterZ", "mm", m_genSphereCenterZ_mm,
      "Mirror /gpd3d/gen/sphereCenterZ (mm)");

  m_runInfoMessenger->DeclareProperty("gpsParticle", m_gpsParticle, "Mirror /gps/particle");
  m_runInfoMessenger->DeclarePropertyWithUnit("gpsEnergy", "keV", m_gpsEnergy_keV, "Mirror /gps/ene/mono (keV)");
}

GPD3D_EventAction::~GPD3D_EventAction()
{
  if (m_rootWriter) {
    m_rootWriter->Close();
    delete m_rootWriter;
    m_rootWriter = nullptr;
  }
  delete m_runInfoMessenger;
  delete m_bgMessenger;
}

// -----------------------------
// hits helper
// -----------------------------
GPD3D_GasCellHitsCollection*
GPD3D_EventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
  if (!event) return nullptr;
  auto hce = event->GetHCofThisEvent();
  if (!hce) return nullptr;
  return static_cast<GPD3D_GasCellHitsCollection*>(hce->GetHC(hcID));
}

// -----------------------------
// BeginOfEvent: attach/reset EventInfo
// -----------------------------
void GPD3D_EventAction::BeginOfEventAction(const G4Event* event)
{
  if (!event) return;
  auto* ev = const_cast<G4Event*>(event);

  auto* info = dynamic_cast<GPD3D_EventInfo*>(ev->GetUserInformation());
  if (!info) {
    info = new GPD3D_EventInfo();
    ev->SetUserInformation(info);
  } else {
    info->Reset();
  }
  edepInGW = false;
  edepInGC = false;
  edepGW = 0.0;
  edepGC = 0.0;

}

// -----------------------------
// Background helpers
// -----------------------------
bool GPD3D_EventAction::IsLastEventOfRun(const G4Event* event, int nBeamOnRequested) const
{
  if (!event) return false;
  const int eid = event->GetEventID();

  if (nBeamOnRequested > 0) {
    return (eid == nBeamOnRequested - 1);
  }

  const auto* run = G4RunManager::GetRunManager()->GetCurrentRun();
  if (!run) return false;
  return (eid == run->GetNumberOfEventToBeProcessed() - 1);
}

void GPD3D_EventAction::EnsureRootOpenAndRunInfo(const G4Event* event, int runID, int nBeamOnRequested)
{
  if (!m_bgWriteROOT) return;

  if (!m_rootWriter) m_rootWriter = new GPD3D_RootWriter();

  // Build output filename (supports %RUN%, and auto-suffix for run>0)
  std::string outName = m_bgRootFile;

  const std::string token = "%RUN%";
  const auto tokenPos = outName.find(token);
  if (tokenPos != std::string::npos) {
    outName.replace(tokenPos, token.size(), std::to_string(runID));
  } else if (runID != 0) {
    const std::string ext = ".root";
    if (outName.size() >= ext.size() &&
        outName.substr(outName.size() - ext.size()) == ext) {
      outName = outName.substr(0, outName.size() - ext.size())
              + "_run" + std::to_string(runID) + ext;
    } else {
      outName = outName + "_run" + std::to_string(runID) + ".root";
    }
  }

  // Open if needed (new run or different filename)
  if (!m_rootWriter->IsOpen() ||
      m_rootWriter->GetRunID() != runID ||
      m_rootWriter->GetFileName() != outName) {
    m_rootWriter->Open(outName, runID);
    m_openRunID = runID;
    m_runInfoFilled = false;
  }

  // Fill RunInfo exactly once per run (even if Events is filtered by hasHit)
  if (!m_runInfoFilled) {
    // If you want to store “actual threads”, you can also query G4RunManager:
    // int actualThreads = G4RunManager::GetRunManager()->GetNumberOfThreads();
    // but in some builds it may differ from requested.
    const int nThreads = (m_reqNThreads > 0 ? m_reqNThreads : -1);

    m_rootWriter->FillRunInfo(
      runID,
      nBeamOnRequested,
      m_controlVerbose,
      m_runVerbose,
      m_eventVerbose,
      m_trackingVerbose,
      nThreads,
      (bool)m_bgMode,
      (bool)m_bgWriteFITS,
      (bool)m_bgWriteROOT,
      std::string(m_bgRootFile),
      (bool)m_genUseSphere,
      m_genSphereRadius_mm,
      m_genSphereCenterX_mm,
      m_genSphereCenterY_mm,
      m_genSphereCenterZ_mm,
      std::string(m_gpsParticle),
      m_gpsEnergy_keV
    );

    m_runInfoFilled = true;
  }

  (void)event; // keep for possible future use
}

// -----------------------------
// Your existing structs/helpers for polarimetry
// -----------------------------
struct Point { double x; double y; double z; };

// (your find_pixel_coordinates_around_center, Diffusion_and_Multiplicaiton,
//  find_physical_coordinates_from_pixel implementation stays the same)
// -----------------------------
// PASTE your existing implementations here unchanged.
// -----------------------------

void GPD3D_EventAction::find_pixel_coordinates_around_center(double x_mm, double y_mm,
    int frame_width_pixels, int frame_height_pixels,
    double pixel_size_microns, int half_region_size,
    int& center_col, int& center_row,
    std::vector<int>& col_list, std::vector<int>& row_list)
{
  double x_microns = x_mm * 1000.0;
  double y_microns = y_mm * 1000.0;

  double x_shifted = x_microns + frame_width_pixels * pixel_size_microns / 2.0;
  double y_shifted = y_microns + frame_height_pixels * pixel_size_microns / 2.0;

  if (x_shifted < 0.0 || x_shifted >= frame_width_pixels * pixel_size_microns ||
      y_shifted < 0.0 || y_shifted >= frame_height_pixels * pixel_size_microns) {
    col_list.clear();
    row_list.clear();
    return;
  }

  center_col = static_cast<int>(x_shifted / pixel_size_microns);
  center_row = static_cast<int>(y_shifted / pixel_size_microns);

  int start_col = std::max(0, center_col - half_region_size);
  int end_col   = std::min(frame_width_pixels, center_col + half_region_size);
  int start_row = std::max(0, center_row - half_region_size);
  int end_row   = std::min(frame_height_pixels, center_row + half_region_size);

  col_list.clear();
  row_list.clear();
  for (int col = start_col; col < end_col; ++col) col_list.push_back(col);
  for (int row = start_row; row < end_row; ++row) row_list.push_back(row);
}

// --- Diffusion_and_Multiplicaiton + find_physical_coordinates_from_pixel ---
// Keep your exact implementation. (Omitted here for brevity in this snippet.)
// In your file, paste the full functions exactly as you had them.

std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<double>>,
           std::vector<int>, std::vector<int>>
GPD3D_EventAction::Diffusion_and_Multiplicaiton(
    const std::vector<double>& ion_pos_x,
    const std::vector<double>& ion_pos_y,
    const std::vector<double>& ion_pos_z,
    int frame_width_pixels, int frame_height_pixels,
    double pixel_size_microns, int half_region_size)
{
  // >>> PASTE YOUR FULL FUNCTION BODY HERE (unchanged) <<<
  // To keep this answer readable, I’m not duplicating the whole block again.
  // (Use exactly the version you posted, it compiles once the background block is fixed.)
  std::vector<std::vector<int>> counting_map(2*half_region_size, std::vector<int>(2*half_region_size,0));
  std::vector<std::vector<double>> map_toa(2*half_region_size, std::vector<double>(2*half_region_size,0.0));
  std::vector<int> first_col_list, first_row_list;
  return std::make_tuple(counting_map, map_toa, first_col_list, first_row_list);
}

std::tuple<std::vector<double>, std::vector<double>>
GPD3D_EventAction::find_physical_coordinates_from_pixel(
    const std::vector<int>& col_list,
    const std::vector<int>& row_list,
    int frame_width_pixels,
    int frame_height_pixels,
    double pixel_size_microns)
{
  std::vector<double> x_mm_list;
  std::vector<double> y_mm_list;
  x_mm_list.reserve(col_list.size());
  y_mm_list.reserve(row_list.size());

  for (size_t i = 0; i < col_list.size(); ++i) {
    x_mm_list.push_back(((col_list[i] * pixel_size_microns
                         - frame_width_pixels * pixel_size_microns / 2.0)
                         + pixel_size_microns / 2.0) / 1000.0);
    y_mm_list.push_back(((row_list[i] * pixel_size_microns
                         - frame_height_pixels * pixel_size_microns / 2.0)
                         + pixel_size_microns / 2.0) / 1000.0);
  }
  return std::make_tuple(x_mm_list, y_mm_list);
}

// -----------------------------
// EndOfEvent
// -----------------------------
void GPD3D_EventAction::EndOfEventAction(const G4Event* event)
{
  if (!event) return;

  // Get hits collections IDs
  if (m_gasCellHCID == -1) {
    m_gasCellHCID = G4SDManager::GetSDMpointer()->GetCollectionID("GasCellHitsCollection");
  }

  // Get hits collections
  auto gasCellHC = GetHitsCollection(m_gasCellHCID, event);

  // ============================
  // BACKGROUND MODE: scoring + ROOT
  // ============================
  if (m_bgMode) {
    const auto* run = G4RunManager::GetRunManager()->GetCurrentRun();
    const int runID = (run ? run->GetRunID() : 0);

    // Robust requested N (avoid MT “random smaller”)
    const int nBeamOnRequested =
      (m_reqNBeamOn > 0 ? m_reqNBeamOn :
       (run ? run->GetNumberOfEventToBeProcessed() : 0));

    // Reset summary counters at first event of each run
    if (event->GetEventID() == 0) {
      m_bgNgen = 0;
      m_bgNhit = 0;
      m_bgSumEdep = 0.;
      m_bgSumEdepHit = 0.;
    }

    m_bgNgen++;

    // Open ROOT + Fill RunInfo once
    EnsureRootOpenAndRunInfo(event, runID, nBeamOnRequested);

    // Event totals in gas
    G4double edepGasTotal = 0.;
    G4double edepGasIon   = 0.;
    G4double edepGasNonIon= 0.;
    G4int    nHitsEdep    = 0;

    if (gasCellHC) {
      const G4int nHits = gasCellHC->entries();
      for (G4int j = 0; j < nHits; ++j) {
        auto* h = (*gasCellHC)[j];
        if (!h) continue;
        const G4double e = h->GetEdep();
        edepGasTotal += e;

        // These require your updated hit class; keep if you added them:
        edepGasIon    += h->GetEdepIon();
        edepGasNonIon += h->GetEdepNonIon();

        if (e > 0.) nHitsEdep++;
      }
    }

    // Check if GAGG wall or cap has hit
    const int GAGGHitWall = (edepInGW && edepGW > 0.) ? 1 : 0;
    const int GAGGHitCap  = (edepInGC && edepGC > 0.) ? 2 : 0;
    const int GAGGHit = GAGGHitWall + GAGGHitCap;

    // Legacy summary (based on total edep)
    m_bgSumEdep += edepGasTotal;
    if (edepGasTotal > 0.) {
      m_bgNhit++;
      m_bgSumEdepHit += edepGasTotal;
    }

    // ROOT output: only if requested
    if (m_bgWriteROOT && m_rootWriter && m_rootWriter->IsOpen()) {
      // Primary metadata
      int primaryPDG = 0;
      double primaryE0_keV = 0.0;
      double vtxX_mm = 0., vtxY_mm = 0., vtxZ_mm = 0.;
      double dirX = 0., dirY = 0., dirZ = 0.;

      const auto* vtx = event->GetPrimaryVertex();
      if (vtx && vtx->GetPrimary()) {
        const auto* p = vtx->GetPrimary();
        primaryPDG = p->GetPDGcode();
        primaryE0_keV = p->GetKineticEnergy() / keV;
        vtxX_mm = vtx->GetX0() / mm;
        vtxY_mm = vtx->GetY0() / mm;
        vtxZ_mm = vtx->GetZ0() / mm;

        const auto d = p->GetMomentumDirection();
        dirX = d.x(); dirY = d.y(); dirZ = d.z();
      }

      // Containment flags from EventInfo
      bool enteredGas=false, exitedGas=false, stoppedInGas=false, containedGas=false;
      if (auto* info = dynamic_cast<GPD3D_EventInfo*>(event->GetUserInformation())) {
        enteredGas = info->primaryEnteredGas;
        exitedGas = info->primaryExitedGas;
        stoppedInGas = info->primaryStoppedInGas;
        containedGas = info->PrimaryContainedInGas();
      }

      const int eventID = event->GetEventID();
      const bool hasHit = (nHitsEdep > 0);

      // Fill Events only for hit-events (your choice)
      if (hasHit) {
        m_rootWriter->FillEventRow(
          runID, eventID, nBeamOnRequested,
          primaryPDG, primaryE0_keV,
          vtxX_mm, vtxY_mm, vtxZ_mm,
          dirX, dirY, dirZ,
          edepGasTotal / keV,
          edepGasIon   / keV,
          edepGasNonIon/ keV,
          nHitsEdep,
          enteredGas, exitedGas, stoppedInGas, containedGas, GAGGHit
        );
      }      
      // Fill Hits only for Edep>0
      if (gasCellHC) {
        const G4int nHits = gasCellHC->entries();
        for (G4int j = 0; j < nHits; ++j) {
          auto* h = (*gasCellHC)[j];
          if (!h) continue;
          if (h->GetEdep() <= 0.) continue;

          const auto hitPos = h->GetPos();
          const auto md  = h->GetMomDir();
          const auto vp  = h->GetVertPos();
          const auto vd  = h->GetVertDir();

          m_rootWriter->FillHitRow(
            runID, eventID,
            h->GetTrackID(), h->GetParentID(), h->GetStepNum(),
            h->GetPDG(),
            hitPos.x()/mm, hitPos.y()/mm, hitPos.z()/mm,
            h->GetEdep()/keV, h->GetEdepIon()/keV, h->GetEdepNonIon()/keV,
            h->GetEnergy()/keV,
            h->GetGlobalTime()/ns,
            h->GetTrkLength()/mm,
            h->IsFirstStep(), h->IsLastStep(),
            md.x(), md.y(), md.z(),
            vp.x()/mm, vp.y()/mm, vp.z()/mm,
            vd.x(), vd.y(), vd.z(),
            h->GetVertEnergy()/keV,
            h->GetCreatorProcType(), h->GetCreatorProcSubType(),
            h->GetStepProcType(), h->GetStepProcSubType()
          );
        }
      }
    }

    // End-of-run summary + close ROOT
    if (IsLastEventOfRun(event, (m_reqNBeamOn > 0 ? m_reqNBeamOn : -1))) {
      const G4double phit = (m_bgNgen > 0) ? (static_cast<G4double>(m_bgNhit) / static_cast<G4double>(m_bgNgen)) : 0.;
      const G4double meanEdep = (m_bgNgen > 0) ? (m_bgSumEdep / static_cast<G4double>(m_bgNgen)) : 0.;
      const G4double meanEdepHit = (m_bgNhit > 0) ? (m_bgSumEdepHit / static_cast<G4double>(m_bgNhit)) : 0.;

      G4cout << "\n========== BACKGROUND MODE SUMMARY ==========" << G4endl;
      G4cout << "RunID: " << runID << G4endl;
      G4cout << "Requested beamOn: " << (m_reqNBeamOn > 0 ? m_reqNBeamOn : nBeamOnRequested) << G4endl;
      G4cout << "Generated events seen here: " << m_bgNgen << G4endl;
      G4cout << "Hit events (EdepGas>0): " << m_bgNhit << G4endl;
      G4cout << "P(hit): " << phit << G4endl;
      G4cout << "Mean EdepGas per primary: " << G4BestUnit(meanEdep, "Energy") << G4endl;
      G4cout << "Mean EdepGas per hit event: " << G4BestUnit(meanEdepHit, "Energy") << G4endl;
      G4cout << "============================================\n" << G4endl;

      if (m_rootWriter && m_rootWriter->IsOpen()) {
        m_rootWriter->Close();
      }
    }

    return; // Do NOT execute polarimetry digitization/output in bg mode.
  }

  // ============================
  // POLARIMETRY MODE (original)
  // ============================

  if (gasCellHC){

    // Get the total number of hits
    G4int numberHits = gasCellHC->entries();
    G4ThreeVector absPos = G4ThreeVector(), lastPos = G4ThreeVector();
    // G4double geoTrackLen = -1., trueTrackLen = -1.;
    G4double peEnergy = -1., augEnergy = -1.;
    G4ThreeVector peDir = G4ThreeVector(), augDir = G4ThreeVector();
    auto primaryVertex = event->GetPrimaryVertex();
    G4double phoEnergy = primaryVertex->GetPrimary()->GetKineticEnergy();
    // G4double time = primaryVertex->GetT0();

    // auto track = GPD3D_IonizationTrack();

    G4int eventID = event->GetEventID();
    // G4cout << "EventID: " << eventID << G4endl;

    // Declare variables to store last information
    G4int lastEventID = 0;
    G4double lastPhoEnergy = 0.0;
    G4double lastAbsPos_X = 0.0;
    G4double lastAbsPos_Y = 0.0;
    G4double lastAbsPos_Z = 0.0;
    G4double lastPEEnergy = 0.0;
    G4double lastPE_PHI = 0.0;
    G4double lastPE_THETA = 0.0;
    G4double lastAUGEnergy = 0.0;
    G4double lastAUG_PHI = 0.0;
    G4double lastAUG_THETA = 0.0;
    G4double lastTrkLength = 0.0;

    std::vector<G4ThreeVector> Track_entire;

    std::vector<G4double> m_pos_X_list_ion;
    std::vector<G4double> m_pos_Y_list_ion;
    std::vector<G4double> m_pos_Z_list_ion;

    for (G4int j = 0; j < numberHits; j++){
      // Get the hit
      m_hit = (*gasCellHC)[j];

      // Get all the relevant information stored
      m_trackID = m_hit->GetTrackID();
      m_parentID = m_hit->GetParentID();
      m_stepNum = m_hit->GetStepNum();
      m_energy = m_hit->GetEnergy();
      m_eDep = m_hit->GetEdep();
      m_pos = m_hit->GetPos();
      m_dir = m_hit->GetMomDir();
      m_trackProcName = m_hit->GetTrackProcName();
      m_stepProcName = m_hit->GetStepProcName();
      m_partName = m_hit->GetPartName();
      m_partState = m_hit->GetPartState();
      m_trkLength = m_hit->GetTrkLength();
      m_ionTrack = m_hit->GetIonTrack();
      m_firstStep = m_hit->IsFirstStep();
      m_lastStep = m_hit->IsLastStep();
      m_vertPos = m_hit->GetVertPos();
      m_vertDir = m_hit->GetVertDir();
      m_vertEnergy = m_hit->GetVertEnergy();

      // If secondary
      if (m_parentID == 1){
        if (m_firstStep == 1){
          absPos = m_vertPos;
          if (m_trackProcName == "phot"){

            // photoelectron energy
            peEnergy = m_vertEnergy;
            // photoelectron direction
            peDir = m_vertDir;
          }
          else if (m_trackProcName == "phot_auger"){
            // Auger energy
            augEnergy = m_vertEnergy;
            // Auger direction
            augDir = m_vertDir;
          }
        }

        // Store information from the last iteration
        lastEventID = eventID;
        lastPhoEnergy = phoEnergy / keV;
        lastAbsPos_X = absPos[0] / mm;
        lastAbsPos_Y = absPos[1] / mm;
        lastAbsPos_Z = absPos[2] / mm;
        lastPEEnergy = peEnergy / keV;
        lastPE_PHI = peDir.phi() / deg;
        lastPE_THETA = peDir.theta() / deg;
        lastAUGEnergy = augEnergy / keV;
        lastAUG_PHI = augDir.phi() / deg;
        lastAUG_THETA = augDir.theta() / deg;
        lastTrkLength = m_trkLength / mm;
      }

      if (m_parentID > 0) {
        // TEST
        Track_entire.push_back(m_pos);

        G4double m_pos_X_ion = m_pos[0] / mm;
        G4double m_pos_Y_ion = m_pos[1] / mm;
        G4double m_pos_Z_ion = m_pos[2] / mm;

        m_pos_X_list_ion.push_back(m_pos_X_ion);
        m_pos_Y_list_ion.push_back(m_pos_Y_ion);
        m_pos_Z_list_ion.push_back(m_pos_Z_ion);
      }
    }

    if ((numberHits > 0) && (peEnergy > 0) && (augEnergy > 0)){
      G4cout<< "++++++++++ SUMMARY ++++++++++"<< G4endl;
      G4cout << "Last Event ID: " << lastEventID << G4endl;
      G4cout << "Last Source Photon Energy: " << lastPhoEnergy << " keV" << G4endl;
      G4cout << "Last Absorption Position X: " << lastAbsPos_X << " mm" << G4endl;
      G4cout << "Last Absorption Position Y: " << lastAbsPos_Y << " mm" << G4endl;
      G4cout << "Last Absorption Position Z: " << lastAbsPos_Z << " mm" << G4endl;
      G4cout << "Last PE Energy: " << lastPEEnergy << " keV" << G4endl;
      G4cout << "Last PE_PHI: " << lastPE_PHI << " deg" << G4endl;
      G4cout << "Last PE_THETA: " << lastPE_THETA << " deg" << G4endl;
      G4cout << "Last AUG Energy: " << lastAUGEnergy << " keV" << G4endl;
      G4cout << "Last AUG_PHI: " << lastAUG_PHI << " deg" << G4endl;
      G4cout << "Last AUG_THETA: " << lastAUG_THETA << " deg" << G4endl;
      G4cout << "Last Track Length: " << lastTrkLength << " mm" << G4endl;
      G4cout<< "+++++++++++++++++++++++++++++"<< G4endl;
      G4cout<<G4endl;

      int frame_width_pixels = 256;
      int frame_height_pixels = 256;
      double pixel_size_microns = 55;
      int half_region_size = 50;

      auto result = Diffusion_and_Multiplicaiton(
        m_pos_X_list_ion, m_pos_Y_list_ion, m_pos_Z_list_ion,
        frame_width_pixels, frame_height_pixels, pixel_size_microns, half_region_size
      );

      std::vector<std::vector<int>> counting_map = std::get<0>(result);
      std::vector<std::vector<double>> map_toa   = std::get<1>(result);
      std::vector<int> first_col_list = std::get<2>(result);
      std::vector<int> first_row_list = std::get<3>(result);

      auto conversion_results = find_physical_coordinates_from_pixel(
        first_col_list, first_row_list,
        frame_width_pixels, frame_height_pixels, pixel_size_microns
      );

      // Access the result
      std::vector<double> x_mm_list = std::get<0>(conversion_results);
      std::vector<double> y_mm_list = std::get<1>(conversion_results);

      // FINAL INFO SCORING
      EVENT_ID_LIST.push_back(lastEventID);
      PHO_ENG_LIST.push_back(lastPhoEnergy);
      ABS_POS_X_LIST.push_back(lastAbsPos_X);
      ABS_POS_Y_LIST.push_back(lastAbsPos_Y);
      ABS_POS_Z_LIST.push_back(lastAbsPos_Z);
      PE_ENG_LIST.push_back(lastPEEnergy);
      PE_PHI_LIST.push_back(lastPE_PHI);
      PE_THETA_LIST.push_back(lastPE_THETA);
      AUG_ENG_LIST.push_back(lastAUGEnergy);
      AUG_PHI_LIST.push_back(lastAUG_PHI);
      AUG_THETA_LIST.push_back(lastAUG_THETA);
      TRK_LEN_LIST.push_back(lastTrkLength);

      ION_POS_X_LIST.push_back(m_pos_X_list_ion);
      ION_POS_Y_LIST.push_back(m_pos_Y_list_ion);
      ION_POS_Z_LIST.push_back(m_pos_Z_list_ion);

      MAP_COL_LIST.push_back(x_mm_list);
      MAP_ROW_LIST.push_back(y_mm_list);

      COUNTING_MAP_3D.push_back(counting_map);
      MAP_TOA_3D.push_back(map_toa);
    }
  }

  // Check if it's the last event
  G4RunManager* runManager = G4RunManager::GetRunManager();
  G4int totalEvents = runManager->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  G4int currentEvent = event->GetEventID();

  if (currentEvent == totalEvents - 1) {
    if (!m_bgWriteFITS) {
      G4cout << "++++++++++ FITS OUTPUT DISABLED ++++++++++" << G4endl;
      return;
    }

    G4cout<< "++++++++++ CREATING OUTPUT.FITS ++++++++++"<< G4endl;

    GPD3D_FITSWriter::writeFITS("output.fits",
          EVENT_ID_LIST,
          PHO_ENG_LIST,
          ABS_POS_X_LIST,
          ABS_POS_Y_LIST,
          ABS_POS_Z_LIST,
          PE_ENG_LIST,
          PE_PHI_LIST,
          PE_THETA_LIST,
          AUG_ENG_LIST,
          AUG_PHI_LIST,
          AUG_THETA_LIST,
          TRK_LEN_LIST,
          ION_POS_X_LIST,
          ION_POS_Y_LIST,
          ION_POS_Z_LIST,
          MAP_COL_LIST,
          MAP_ROW_LIST,
          COUNTING_MAP_3D,
          MAP_TOA_3D);

    int totalProgress = 41; // Total progress steps
    for (int i = 0; i <= totalProgress; ++i) {
        std::cout << "+";
        std::cout.flush();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    std::cout << std::endl;
  }

}
