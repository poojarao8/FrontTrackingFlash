! Print out the front restart in format.

subroutine FrontTracking_restartWrite(filename)
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none

    character :: filename
    logical, save :: testUseFrontTracking
    
    call RuntimeParameters_get ("useFrontTracking", testUseFrontTracking)
    
    if (testUseFrontTracking) then
        call Driver_abortFlash("FrontTracking unit seems not to be compiled in, and the FrontTracking_init stub does not allow the value of useFrontTracking to be TRUE.")
    end if

end subroutine FrontTracking_restartWrite

