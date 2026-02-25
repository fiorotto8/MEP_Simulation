/*!
  @file
  @brief User action initialization.
*/

#ifndef GPD3D_ACTIONINITIALIZATION_HH
#define GPD3D_ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"


//! Action initialization class.
class GPD3D_ActionInitialization : public G4VUserActionInitialization
{
  
 public:

  //! Constructor.
  GPD3D_ActionInitialization();

  //! Overloaded method to instantiate the user run action for master thread.
  virtual void BuildForMaster() const;

  //! Overloaded method to instantiate user action class objects.
  virtual void Build() const;

};

#endif //GPD3D_ACTIONINITIALIZATION_HH
