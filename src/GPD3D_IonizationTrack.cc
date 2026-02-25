#include "GPD3D_IonizationTrack.hh"


/*!
 */
GPD3D_IonizationTrack::GPD3D_IonizationTrack()
{}


/*!
 */
void GPD3D_IonizationTrack::append(G4ThreeVector position)
{
  m_positionVector.push_back(position);
}


/*!
 */
void GPD3D_IonizationTrack::append(const std::vector<G4ThreeVector>& positions)
{
  m_positionVector.insert(m_positionVector.end(), positions.begin(),
                          positions.end());
}


/*!
 */
void GPD3D_IonizationTrack::append(const GPD3D_IonizationTrack& track)
{
  append(track.positionVector());
}


/*!
 */
G4ThreeVector GPD3D_IonizationTrack::position(int index) const
{
  return m_positionVector.at(index);
}
