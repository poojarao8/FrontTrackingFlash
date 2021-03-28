subroutine FrontTracking_updateComp()

#include "Flash.h"
#include "constants.h"
#include "fatpi.h"

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                        Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none

    real, pointer                           ::  solnData(:,:,:,:)
    integer                                 ::  comp

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
    real, dimension(MDIM)                   ::  coords

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
    if (NDIM .GE. 3) then 
        call Grid_getCellCoords(KAXIS,blkID,CENTER,gcell,kCoord,ksizeGC)
    else
        kCoord(:) = 0.0
    end if

    call Grid_getBlkPtr(blkID,solnData)

    do icell = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        coords(IAXIS) = iCoord(icell)
        do jcell = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
            coords(JAXIS) = jCoord(jcell)
            do kcell = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                coords(KAXIS) = kCoord(kcell)
                ! TODO: FIXME    
                call FTAPI_getComponent(coords, comp)
                !comp = 2
                solnData(COMP_VAR, icell, jcell, kcell) = comp
            end do
        end do
    end do

end subroutine
