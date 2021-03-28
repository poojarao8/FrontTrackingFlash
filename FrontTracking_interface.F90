Module FrontTracking_interface

  implicit none

  interface
     subroutine FrontTracking_init()
     end subroutine FrontTracking_init
  end interface

  interface
     subroutine FrontTracking_propagate(dt)
         REAL,intent(in) :: dt
     end subroutine FrontTracking_propagate
  end interface

  interface
     subroutine FrontTracking_output(step,chk)
         integer, intent(in) :: step
         integer, intent(in) :: chk
     end subroutine FrontTracking_output
  end interface

  interface
     subroutine FrontTracking_restartWrite(checkpointID)
         integer, intent(inout) :: checkpointid
     end subroutine FrontTracking_restartWrite
  end interface

  interface
     subroutine FrontTracking_levelFunc(tmp,coords)
         real, dimension(3), intent(in) :: coords
         real, POINTER :: tmp
     end subroutine FrontTracking_levelFunc
  end interface

  interface
     subroutine FrontTracking_getVel(tmp,coords,vel,lcomp,rcomp)
         real, dimension(3), intent(in) :: coords
         real, dimension(3), intent(out) :: vel
         real, POINTER :: tmp
         real, intent(in) :: lcomp
         real, intent(in) :: rcomp
     end subroutine FrontTracking_getVel
  end interface

end Module FrontTracking_interface

