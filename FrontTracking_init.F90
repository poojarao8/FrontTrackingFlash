! This is the Front Tracking initialization function
! to load the data structures necessary for Front Tracking.
! This function uses the grid and driver interface to get information about
! the simulation domain and parallel decomposition, etc., and thus should only
! be called after those modules have initialized.

subroutine FrontTracking_init()
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

!!$  logical, INTENT(in) :: restart

  logical, save :: testUseFrontTracking

  !! It is a failure to invoke the stub when useParticles is set TRUE.

  call RuntimeParameters_get ("useFrontTracking", testUseFrontTracking)
  if (testUseFrontTracking) then
     call Driver_abortFlash("FrontTracking unit seems not to be compiled in,& 
     	  &and the FrontTracking_init stub does not &
          &allow the value of useFrontTracking to be TRUE.")
  end if
end subroutine FrontTracking_init
