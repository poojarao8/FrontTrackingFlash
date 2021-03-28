subroutine FrontTracking_getVel(tmp, coords, velocity, lcomp, rcomp)

#include "Flash.h"
#include "constants.h"

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                        Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none
    real, POINTER, INTENT(IN)                       ::  tmp
    real, INTENT(IN), dimension(MDIM)               ::  coords
    real, INTENT(OUT), dimension(MDIM)              ::  velocity
    real, INTENT(IN)                                ::  lcomp   !UNUSED HERE
    real, INTENT(IN)                                ::  rcomp   !UNUSED HERE

    real, pointer                           ::  solnData(:,:,:,:)

#ifdef FIXEDBLOCKSIZE
    real,DIMENSION(MAXCELLS),target         ::  iCoord, jCoord, kCoord
#else
    real, allocatable,DIMENSION(:),target   ::  iCoord, jCoord, kCoord
#endif

    integer                                 ::  blkID, isizeGC, jsizeGC, ksizeGC
    integer                                 ::  icell, jcell, kcell
    integer, dimension(2, MDIM)             ::  blkLimits, blkLimitsGC
    logical                                 ::  gcell = .true.
    real, dimension(MDIM)                   ::  del

    ! JAM -- Added for FTAPI   FIXME: I would rather have this passed in than
    !                                 hardcoded, it is not compatible with 
    !                                 a lot of the options in FLASH
    blkID = 1

    call Grid_getDeltas(blkID,del)
    call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC)
    iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
    jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
    kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
    numCells = max(isizeGC,jsizeGC)
    numCells = max(numCells,ksizeGC)
 
    allocate(iCoord(numCells))
    allocate(jCoord(numCells))
    allocate(kCoord(numCells))
#endif

    call Grid_getCellCoords(IAXIS,blkID,CENTER,gcell,iCoord,isizeGC)
    if (NDIM .GE. 2) then 
        call Grid_getCellCoords(JAXIS,blkID,CENTER,gcell,jCoord,jsizeGC)
    else
        jCoord(:) = 0.0
    end if        
    if (NDIM .GE. 2) then 
        call Grid_getCellCoords(KAXIS,blkID,CENTER,gcell,kCoord,ksizeGC)
    else
        kCoord(:) = 0.0
    end if

    call Grid_getBlkPtr(blkID,solnData)

    !! FIXME: This needs to be checked, it may not work?
    !! NOTE: IAXIS = 1, JAXIS = 2, KAXIS = 3, use these for FLASH syntax
    !! Need the + 1 to work, since FORTRAN starts indexing at 1
    icell = int(0.5 + (coords(IAXIS) - iCoord(1))/del(IAXIS)) + 1
    if (NDIM .GE. 2) then
        jcell = int(0.5 + (coords(JAXIS) - jCoord(1))/del(JAXIS)) + 1
    else 
        jcell = 1
    end if
    if (NDIM .GE. 3) then
        kcell = int(0.5 + (coords(KAXIS) - kCoord(1))/del(KAXIS)) + 1
    else
        kcell = 1
    end if

    ! To fix when the coord is on the up boundaries as the above code
    ! doesn't quite handle it properly
    icell = min(icell,MAXCELLS)
    jcell = min(jcell,MAXCELLS)
    kcell = min(kcell,MAXCELLS)

! Code below has rounding problems for the sanity check since FT can have coords
! which fall off the domain... FIXME: Why does FT have coords off the domain???
    !! Sanity check, make sure that the coords are in the cell we have selected
!    if (.NOT. (iCoord(icell) - 0.5000000000001*del(IAXIS) .LE. coords(IAXIS) .AND. &
!               iCoord(icell) + 0.5000000000001*del(IAXIS) .GE. coords(IAXIS))) then
!        call Driver_abortFlash('[FTI_getVelocity] ERROR: icell != icoord(IAXIS)') 
!    end if
!    if (NDIM .GE. 2) then
!        if (.NOT. (jCoord(jcell) - 0.5000000000001*del(JAXIS) .LE. coords(JAXIS) .AND. &
!                   jCoord(jcell) + 0.5000000000001*del(JAXIS) .GE. coords(JAXIS))) then
!            call Driver_abortFlash('[FTI_getVelocity] ERROR: jcell != icoord(JAXIS)') 
!        end if
!    end if
!    if (NDIM .GE. 3) then
!        if (.NOT. (kCoord(kcell) - 0.5000000000001*del(KAXIS) .LE. coords(KAXIS) .AND. &
!                   kCoord(kcell) + 0.5000000000001*del(KAXIS) .GE. coords(KAXIS))) then
!            call Driver_abortFlash('[FTI_getVelocity] ERROR: kcell != icoord(KAXIS)') 
!        end if
!    end if

    !! If we make it to here, then we have chosen the right cell
    velocity(IAXIS) = solnData(VELX_VAR, icell, jcell, kcell)
    velocity(JAXIS) = solnData(VELY_VAR, icell, jcell, kcell)
    velocity(KAXIS) = solnData(VELZ_VAR, icell, jcell, kcell)

    call Grid_releaseBlkPtr(blkID,solnData)

#ifndef FIXEDBLOCKSIZE
    deallocate(iCoord)
    deallocate(jCoord)
    deallocate(kCoord)
#endif

    return

end subroutine FrontTracking_getVel
