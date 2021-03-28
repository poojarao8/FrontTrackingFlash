! This function will propagate the front based on the interior velocities
! The velocities are accessed using the grid interface.
subroutine FrontTracking_getVel(%VAL(0),coords,vel,lcomp,rcomp)
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Driver_interface, ONLY : Driver_abortFlash

#include "Flash.h"
#include "constants.h"

    implicit none

    REAL,dimension(MDIM),intent(in) :: coords
    REAL,dimension(MDIM),intent(out) :: vel
    REAL,intent(in) :: lcomp
    REAL,intent(out) :: rcomp
    logical, save :: testUseFrontTracking
    
    call RuntimeParameters_get ("useFrontTracking", testUseFrontTracking)
    
    if (testUseFrontTracking) then
        call Driver_abortFlash("FrontTracking unit seems not to be compiled in, and the FrontTracking_init stub does not allow the value of useFrontTracking to be TRUE.")
    end if

end subroutine FrontTracking_getVel
