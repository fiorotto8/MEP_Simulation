#ifndef GPD3D_EVENTACTION_HH
#define GPD3D_EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"

#include <tuple>
#include <vector>
#include <string>

#include "GPD3D_GasCellHit.hh"
#include "GPD3D_IonizationTrack.hh"

// ROOT writer (your updated one)
#include "GPD3D_RootWriter.hh"

// Forward declare
class G4Event;

class GPD3D_EventAction : public G4UserEventAction {
public:
  GPD3D_EventAction();
  ~GPD3D_EventAction() override;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  // --- Your existing helpers (polarimetry path) ---
  void find_pixel_coordinates_around_center(
      double x_mm, double y_mm,
      int frame_width_pixels, int frame_height_pixels,
      double pixel_size_microns, int half_region_size,
      int& center_col, int& center_row,
      std::vector<int>& col_list, std::vector<int>& row_list);

  std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<double>>,
             std::vector<int>, std::vector<int>>
  Diffusion_and_Multiplicaiton(
      const std::vector<double>& ion_pos_x,
      const std::vector<double>& ion_pos_y,
      const std::vector<double>& ion_pos_z,
      int frame_width_pixels, int frame_height_pixels,
      double pixel_size_microns, int half_region_size);

  std::tuple<std::vector<double>, std::vector<double>>
  find_physical_coordinates_from_pixel(
      const std::vector<int>& col_list,
      const std::vector<int>& row_list,
      int frame_width_pixels,
      int frame_height_pixels,
      double pixel_size_microns);

  bool edepInGW = false;
  bool edepInGC = false;
  double edepGW = 0.0;
  double edepGC = 0.0;

private:
  // --- Hits collection helpers ---
  GPD3D_GasCellHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;

  // --- Background mode helpers ---
  void EnsureRootOpenAndRunInfo(const G4Event* event, int runID, int nBeamOnRequested);
  bool IsLastEventOfRun(const G4Event* event, int nBeamOnRequested) const;

  // =========================
  // 1) NORMAL (polarimetry) state (your existing members)
  // =========================
  G4int m_gasCellHCID;
  GPD3D_GasCellHit* m_hit;

  G4int m_trackID;
  G4int m_parentID;
  G4int m_stepNum;
  G4double m_eDep;
  G4double m_energy;
  G4double m_trkLength;
  G4ThreeVector m_pos;
  G4ThreeVector m_dir;
  G4String m_trackProcName;
  G4String m_stepProcName;
  G4String m_partName;
  G4TrackStatus m_partState;
  GPD3D_IonizationTrack m_ionTrack;
  G4bool m_firstStep;
  G4bool m_lastStep;
  G4ThreeVector m_vertPos;
  G4ThreeVector m_vertDir;
  G4double m_vertEnergy;

  std::vector<G4int>    EVENT_ID_LIST;
  std::vector<G4double> PHO_ENG_LIST;
  std::vector<G4double> ABS_POS_X_LIST;
  std::vector<G4double> ABS_POS_Y_LIST;
  std::vector<G4double> ABS_POS_Z_LIST;
  std::vector<G4double> PE_ENG_LIST;
  std::vector<G4double> PE_PHI_LIST;
  std::vector<G4double> PE_THETA_LIST;
  std::vector<G4double> AUG_ENG_LIST;
  std::vector<G4double> AUG_PHI_LIST;
  std::vector<G4double> AUG_THETA_LIST;
  std::vector<G4double> TRK_LEN_LIST;

  std::vector<std::vector<G4double>> ION_POS_X_LIST;
  std::vector<std::vector<G4double>> ION_POS_Y_LIST;
  std::vector<std::vector<G4double>> ION_POS_Z_LIST;

  std::vector<std::vector<G4double>> MAP_COL_LIST;
  std::vector<std::vector<G4double>> MAP_ROW_LIST;

  std::vector<std::vector<std::vector<int>>>    COUNTING_MAP_3D;
  std::vector<std::vector<std::vector<G4double>>> MAP_TOA_3D;

  // =========================
  // 2) BACKGROUND MODE CONTROLS + RUNINFO CAPTURE
  // =========================
  G4bool   m_bgMode;
  G4bool   m_bgWriteFITS;
  G4bool   m_bgWriteROOT;
  G4String m_bgRootFile;

  // requested beamOn (robust vs MT splitting)
  G4int    m_reqNBeamOn;

  // Optional “runinfo mirror” values (set them in macro via /gpd3d/runinfo/*)
  G4int    m_controlVerbose;
  G4int    m_runVerbose;
  G4int    m_eventVerbose;
  G4int    m_trackingVerbose;

  G4int    m_reqNThreads;

  G4bool   m_genUseSphere;
  G4double m_genSphereRadius_mm;
  G4double m_genSphereCenterX_mm;
  G4double m_genSphereCenterY_mm;
  G4double m_genSphereCenterZ_mm;

  // You can mirror /gps/particle and /gps/ene/mono here if you want
  G4String m_gpsParticle;
  G4double m_gpsEnergy_keV;

  G4GenericMessenger* m_bgMessenger;
  G4GenericMessenger* m_runInfoMessenger;

  // background summary counters
  G4long   m_bgNgen;
  G4long   m_bgNhit;
  G4double m_bgSumEdep;
  G4double m_bgSumEdepHit;

  // ROOT writer
  GPD3D_RootWriter* m_rootWriter;
  bool m_runInfoFilled;
  int  m_openRunID;
};

#endif
