! This is the Front Tracking output function
! to print front information for Front Tracking.

subroutine FrontTracking_output(step)
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none

    integer, intent(in) :: step
    logical, save :: testUseFrontTracking

    call RuntimeParameters_get ("useFrontTracking", testUseFrontTracking)
    
    if (testUseFrontTracking) then
        call Driver_abortFlash("FrontTracking unit seems not to be compiled in, and the FrontTracking_init stub does not allow the value of useFrontTracking to be TRUE.")
    end if

end subroutine FrontTracking_output
