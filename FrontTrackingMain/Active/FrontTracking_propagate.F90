#define aindex(i,j,k,isize,jsize,ksize) (((((k)-1)*(jsize)+((j)-1))*(isize))+((i)-1) + 1 )
subroutine FrontTracking_propagate(grid_dt)
    use Grid_interface, only : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
        Grid_getLocalNumBlks, Grid_getListOfBlocks, Grid_fillGuardCells
    use FrontTracking_data, ONLY: fr_useFrontTracking
    use Grid_data, ONLY : gr_delta
    use Driver_data, ONLY : dr_simTime, dr_dt, dr_dtOld, dr_dtNew, dr_nstep

    implicit none
#include "Flash.h"
#include "constants.h"
#include "ftapi.h"

    real, INTENT(IN) :: grid_dt
    real, dimension(:), allocatable :: xvel,yvel,zvel, xCoord
    real*8, dimension(3) :: point
    real :: maxRho , minRho, front_dt
    real*8, pointer :: solnData(:,:,:,:)
    integer,dimension(MAXBLOCKS) :: blockList
    integer :: gridBlock, counter
    integer :: localNumBlocks, blockCount
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: isize,jsize,ksize
    integer :: iisize,ijsize,iksize
    integer :: i,j,k
    character(len=1024) :: filename
    integer :: ind,aind, surfType

    ! JAM - Added to have ghost states for mass scalars stored
    real, DIMENSION(NSPECIES+NMASS_SCALARS,   &
                    GRID_ILO_GC:GRID_IHI_GC,  &
                    GRID_JLO_GC:GRID_JHI_GC,  &
                    GRID_KLO_GC:GRID_KHI_GC,  &
                    2) :: xnGhostVal 

    ! Not sure what this is used for on Front Tracking end
    surfType = 1

    if(.NOT. fr_useFrontTracking) then
        return
    end if

    call Grid_fillGuardCells(CENTER_FACES,ALLDIR)

    call Grid_getLocalNumBlks(localNumBlocks)
    call Grid_getListOfBlocks(LEAF,blockList,blockCount)

    if(localNumBlocks .ne. 1) then
        write(6,*) "Error in FrontTracking_init(): Too many blocks. ",&
            "Only single block supported!!"
        call exit(-1)
    end if

    gridBlock = blockList(1)

    call Grid_getBlkIndexLimits(gridBlock,blkLimits,blkLimitsGC)

    ! TEMP
    iisize=blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
    ijsize=blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
    iksize=blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

    isize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
    jsize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
    ksize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

    allocate(xvel( (isize+1)*(jsize+1)*(ksize+1)))
    allocate(yvel( (isize+1)*(jsize+1)*(ksize+1)))
    allocate(zvel( (isize+1)*(jsize+1)*(ksize+1)))

    call Grid_getBlkPtr(gridBlock,solnData)

    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        ind=aindex(i,j,k,isize,jsize,ksize)
        xvel(ind) = solnData(VELX_VAR, i, j, k)
        yvel(ind) = solnData(VELY_VAR, i, j, k)
        zvel(ind) = solnData(VELZ_VAR, i, j, k)
    end do
    end do
    end do

    ! JAM - 5/9/16  We may need to subcycle the front if it's dt is smaller 
    ! than FLASH's dt, so we will do that here (it needs to be on the client
    ! side so that updateComp can be called before the next iteration through)
    front_dt = 0.0
    counter = 0
    do while (front_dt .lt. grid_dt)
        call FrontTracking_fillGhostStates(xnGhostVal)

        call Grid_fillGuardCells(CENTER_FACES, ALLDIR)

        call FTAPI_propagate_cartesian(grid_dt,xvel,yvel,zvel,front_dt)

        call FrontTracking_updateComp(xnGhostVal)

        counter = counter + 1
        if (counter .gt. 20) then
            print *, 'The front has tried to propagate more than  &
			& 20 times, this suggests an issue with the front that may  &
			& need to be resolved.  Sometimes a restart may be all that &
			& is needed'
           call Driver_abortFlash('See logfile for additional info, FrontTracking_Propagate.F90')
        end if
    end do

    !call ft_get_point(point)
    allocate(xCoord(blkLimitsGC(HIGH, IAXIS)))

    call Grid_getCellCoords(IAXIS, blockList(1), CENTER, .true., xCoord,&
                 blkLimitsGC(HIGH, IAXIS))
                 
    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        if( 0.001 > xCoord(i)) then
                maxRho = maxval(solnData(DENS_VAR,i:blkLimits(HIGH,IAXIS),1,1))
                minRho = minval(solnData(DENS_VAR,i:blkLimits(HIGH,IAXIS),1,1))
        end if
        exit
    end do
    
    deallocate(xvel) 
    deallocate(yvel) 
    deallocate(zvel) 

end subroutine FrontTracking_propagate
