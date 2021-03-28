!!****if* source/Simulation/SimulationMain/Sod/Simulation_levelFunc
!!
!! NAME
!!
!!  Simulation_levelFunc
!!
!!***

real function FrontTracking_levelFunc(tmp, coords)

  use Simulation_data, ONLY: sim_posn, sim_xCos, sim_yCos, sim_zCos  

#include "Flash.h"
#include "constants.h"

  implicit none
  
  real, intent(in), dimension(MDIM) :: coords
  real, intent(in) :: tmp
  real :: xx, yy, zz
!  real, intent(out) :: level  
  real :: lPosn0, lPosn

  zz = coords(3) 
  lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
  yy = coords(2)
  lPosn = lPosn0 - yy*sim_yCos/sim_xCos
  xx = coords(1) 

  FrontTracking_levelFunc = yy - lPosn
  return

end function FrontTracking_levelFunc
