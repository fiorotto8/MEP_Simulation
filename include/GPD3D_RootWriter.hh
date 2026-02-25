#ifndef GPD3D_ROOTWRITER_HH
#define GPD3D_ROOTWRITER_HH

#include <string>

class TFile;
class TTree;

class GPD3D_RootWriter {
public:
  GPD3D_RootWriter();
  ~GPD3D_RootWriter();

  // Open/close
  bool Open(const std::string& filename, int runID);
  void Close();

  bool IsOpen() const { return (m_file != nullptr); }
  int  GetRunID() const { return m_runID; }
  const std::string& GetFileName() const { return m_filename; }

  // Run-level info (filled once)
  void FillRunInfo(
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
  );

  // Event-level row (fill only for “hit” events if you want)
  void FillEventRow(
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
  );

  // Hit-level row
  void FillHitRow(
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
  );

private:
  void BookTrees_();
  void ResetBuffers_();

private:
  // ROOT handles
  TFile* m_file;
  TTree* m_events;
  TTree* m_hits;
  TTree* m_runInfo;

  int m_runID;
  std::string m_filename;

  // -------- RunInfo buffers --------
  bool m_runInfoFilled;

  int  b_run_runID;
  int  b_run_nBeamOnRequested;
  int  b_run_controlVerbose;
  int  b_run_runVerbose;
  int  b_run_eventVerbose;
  int  b_run_trackingVerbose;
  int  b_run_nThreads;

  bool b_run_bgMode;
  bool b_run_bgWriteFITS;
  bool b_run_bgWriteROOT;

  std::string b_run_bgRootFile;

  bool b_run_genUseSphere;
  double b_run_genSphereRadius_mm;
  double b_run_genSphereCenterX_mm;
  double b_run_genSphereCenterY_mm;
  double b_run_genSphereCenterZ_mm;

  std::string b_run_gpsParticle;
  double b_run_gpsEnergy_keV;

  // -------- Events buffers --------
  int  b_evt_runID;
  int  b_evt_eventID;
  int  b_evt_nBeamOn;
  int  b_evt_primaryPDG;
  double b_evt_primaryE0_keV;
  double b_evt_vtxX_mm, b_evt_vtxY_mm, b_evt_vtxZ_mm;
  double b_evt_dirX, b_evt_dirY, b_evt_dirZ;
  double b_evt_edepGasTotal_keV;
  double b_evt_edepGasIon_keV;
  double b_evt_edepGasNonIon_keV;
  int  b_evt_nHitsEdep;
  bool b_evt_enteredGas, b_evt_exitedGas, b_evt_stoppedInGas, b_evt_containedGas;
  int b_GAGGHit;

  // -------- Hits buffers --------
  int b_hit_runID, b_hit_eventID;
  int b_hit_trackID, b_hit_parentID, b_hit_stepNum;
  int b_hit_pdg;
  double b_hit_x_mm, b_hit_y_mm, b_hit_z_mm;
  double b_hit_edep_keV, b_hit_edepIon_keV, b_hit_edepNonIon_keV;
  double b_hit_energy_keV;
  double b_hit_globalTime_ns;
  double b_hit_trkLength_mm;
  bool b_hit_firstStep, b_hit_lastStep;
  double b_hit_dirX, b_hit_dirY, b_hit_dirZ;
  double b_hit_vtxX_mm, b_hit_vtxY_mm, b_hit_vtxZ_mm;
  double b_hit_vtxDirX, b_hit_vtxDirY, b_hit_vtxDirZ;
  double b_hit_vtxEnergy_keV;
  int b_hit_creatorProcType, b_hit_creatorProcSubType;
  int b_hit_stepProcType, b_hit_stepProcSubType;
};

#endif
