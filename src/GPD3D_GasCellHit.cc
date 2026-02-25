#include <iomanip>

#include "G4UnitsTable.hh"
#include "GPD3D_GasCellHit.hh"

G4ThreadLocal G4Allocator<GPD3D_GasCellHit>* GPD3D_GasCellHitAllocator=0;

GPD3D_GasCellHit::GPD3D_GasCellHit() :
  G4VHit(),
  m_trackID(-1),
  m_parentID(-1),
  m_stepNum(-1),
  m_eDep(0.),
  m_energy(0.),
  m_trkLength(0.),
  m_pos(G4ThreeVector()),
  m_momDir(G4ThreeVector()),
  m_trackProcName(""),
  m_stepProcName(""),
  m_partName(""),
  m_partState(fAlive),
  m_ionTrack(),
  m_firstStep(false),
  m_lastStep(false),
  m_vertPos(G4ThreeVector()),
  m_vertDir(G4ThreeVector()),
  m_vertEnergy(-1.),
  // new
  m_eDepIon(0.),
  m_eDepNonIon(0.),
  m_pdg(0),
  m_globalTime(0.),
  m_creatorProcType(-1),
  m_creatorProcSubType(-1),
  m_stepProcType(-1),
  m_stepProcSubType(-1)
{}

GPD3D_GasCellHit::~GPD3D_GasCellHit() {}

GPD3D_GasCellHit::GPD3D_GasCellHit(const GPD3D_GasCellHit& right) :
  G4VHit()
{
  m_trackID = right.m_trackID;
  m_parentID = right.m_parentID;
  m_stepNum = right.m_stepNum;
  m_eDep = right.m_eDep;
  m_energy = right.m_energy;
  m_trkLength = right.m_trkLength;
  m_pos = right.m_pos;
  m_momDir = right.m_momDir;
  m_trackProcName = right.m_trackProcName;
  m_stepProcName = right.m_stepProcName;
  m_partName = right.m_partName;
  m_partState = right.m_partState;
  m_ionTrack = right.m_ionTrack;
  m_firstStep = right.m_firstStep;
  m_lastStep = right.m_lastStep;
  m_vertPos = right.m_vertPos;
  m_vertDir = right.m_vertDir;
  m_vertEnergy = right.m_vertEnergy;

  m_eDepIon = right.m_eDepIon;
  m_eDepNonIon = right.m_eDepNonIon;
  m_pdg = right.m_pdg;
  m_globalTime = right.m_globalTime;
  m_creatorProcType = right.m_creatorProcType;
  m_creatorProcSubType = right.m_creatorProcSubType;
  m_stepProcType = right.m_stepProcType;
  m_stepProcSubType = right.m_stepProcSubType;
}

const GPD3D_GasCellHit& GPD3D_GasCellHit::operator=(const GPD3D_GasCellHit& right)
{
  m_trackID = right.m_trackID;
  m_parentID = right.m_parentID;
  m_stepNum = right.m_stepNum;
  m_eDep = right.m_eDep;
  m_energy = right.m_energy;
  m_trkLength = right.m_trkLength;
  m_pos = right.m_pos;
  m_momDir = right.m_momDir;
  m_trackProcName = right.m_trackProcName;
  m_stepProcName = right.m_stepProcName;
  m_partName = right.m_partName;
  m_partState = right.m_partState;
  m_ionTrack = right.m_ionTrack;
  m_firstStep = right.m_firstStep;
  m_lastStep = right.m_lastStep;
  m_vertPos = right.m_vertPos;
  m_vertDir = right.m_vertDir;
  m_vertEnergy = right.m_vertEnergy;

  m_eDepIon = right.m_eDepIon;
  m_eDepNonIon = right.m_eDepNonIon;
  m_pdg = right.m_pdg;
  m_globalTime = right.m_globalTime;
  m_creatorProcType = right.m_creatorProcType;
  m_creatorProcSubType = right.m_creatorProcSubType;
  m_stepProcType = right.m_stepProcType;
  m_stepProcSubType = right.m_stepProcSubType;

  return *this;
}

G4int GPD3D_GasCellHit::operator==(const GPD3D_GasCellHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}
