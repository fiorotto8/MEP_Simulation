/*!
  @file 
  @brief Primary ionization track.
*/

#ifndef GPD3D_IONIZATIONTRACK_H
#define GPD3D_IONIZATIONTRACK_H


#include <vector>

#include "G4ThreeVector.hh"


//! Class describing the ionization created in the gas by a charged particle.
/*!
  This is one of the main data structures used in the simulation and 
  encapsulate the properties of a generic ionizaition track. The basic 
  class member is a standard vector of three-dimensional positions
  (precisely, G4ThreeVector objects) indicating the points where the 
  ionization electrons are located.

  The class provide facilities to access the position of a given electron in
  the track, as well as functions to append to an existing track either
  a single three-dimensional position or another track. In addition, the 
  += operator is overloaded so that one can concatenate track by just
  adding them.
 */
class GPD3D_IonizationTrack
{
 public:

  //! Constructor.
  GPD3D_IonizationTrack();

  //! Return the underline vector.
  const std::vector<G4ThreeVector>& positionVector() const {
    return  m_positionVector;
  }

  //! Clear the underlying vector of positions.
  void clear() {m_positionVector.clear();}

  //! Append a position to the track.
  void append(G4ThreeVector position);

  //! Append a position to the track passing the three coordinates explicitely.
  void append(double x, double y, double z) {append(G4ThreeVector(x, y, z));}
  
  //! Append a vector of positions to the track.
  void append(const std::vector<G4ThreeVector>& positions);

  //! Append another GPD3D_Ionizationtrack object.
  void append(const GPD3D_IonizationTrack& track);

  //! Overload of the += operator.
  GPD3D_IonizationTrack operator+=(const GPD3D_IonizationTrack& t) {
    this->append(t);
    return *this;
  }
  
  //! Retrieve a position from the track.
  G4ThreeVector position(int index) const;

  double x(int index) const {return position(index).x();}
  double y(int index) const {return position(index).y();}
  double z(int index) const {return position(index).z();}

  //! Return the number of primary electrons in the track.
  int size() const {return m_positionVector.size();}


 private:

  std::vector<G4ThreeVector> m_positionVector;

};


#endif //GPD3D_IONIZATIONTRACK_H
