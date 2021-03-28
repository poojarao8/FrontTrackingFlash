! NAME: FrontTracking init
! DESCRIPTION: This initializes the Front Tracking module, FronTier
! SYNOPSIS: This routine accesses the domain parameters via the grid interface
!               and initializes the FrontTracking domain and grid, and also
!               the initial front.


#define aindex(i,j,k,isize,jsize,ksize) (((((k)-1)*(jsize)+((j)-1))*(isize))+((i)-1) + 1)
subroutine FrontTracking_init()
    use RuntimeParameters_interface, only : RuntimeParameters_get
    use Driver_interface, only : Driver_abortFlash
    use Grid_interface, only :    Grid_getBlkPtr, Grid_getBlkIndexLimits,&
                                  Grid_getLocalNumBlks, Grid_getListOfBlocks,&
                                  Grid_getBlkBoundBox, Grid_releaseBlkPtr, Grid_getGeometry, &
                                  Grid_getGlobalIndexLimits, Grid_fillGuardCells
    use FrontTracking_data, ONLY: fr_useFrontTracking, fr_restart, FrontTracking_compGrid
    use FrontTracking_interface, only : FrontTracking_levelFunc, FrontTracking_getVel
    
    implicit none

#include "constants.h"
#include "FTAPI.h"
#include "Flash.h"

    ! The list of blocks local to this processor
    integer,dimension(MAXBLOCKS) :: blockList 
    ! The number of blocks on this processor
    ! Should be =1 after function call, since FT works only on uniform grid
    integer :: blockCount

    ! The grid block with a front on it
    integer :: thisGridBlock

    ! Integer grid limits
    integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC 

    ! Pointer to interior states
    real, pointer :: solnData(:,:,:,:)
    ! Linear storage arrays to pass into FT for velocity field.
    real, dimension(:), allocatable :: xvel,yvel,zvel

    ! Number of cells in direction
    integer :: isize ,jsize,ksize
    ! Nunber of interior (non-buffer) cells in direction
    integer :: iisize ,ijsize,iksize
    ! Buffer size in direction
    integer :: ibuf, jbuf, kbuf
    
    ! Temporary storage for computed index
    integer :: ind
    real, dimension(MDIM) :: coords
    REAL , dimension(2,MDIM) :: boundBox
    integer, dimension(MDIM) :: indexLims
    INTEGER :: i,j,k
    integer,dimension(2,MDIM) :: bdry, bdry1

    integer, dimension(MDIM) :: procGrid
    !integer, dimension(MDIM) :: comm, procGrid,me
    integer ::    numDim

    !Inidicator for type of Geometry
    integer :: geometry, aind
    !logical :: fr_restart
    integer :: restart_num
    character*80 :: restart_base
    character*80 :: restart_filename
    character*80 :: restart_num_str
    integer :: rest_len 
    
    call RuntimeParameters_get('useFrontTracking', fr_useFrontTracking)
    call RuntimeParameters_get('restart', fr_restart)

    call Grid_getListOfBlocks(LEAF,blockList,blockCount)
    
    if(blockCount .eq. 0) then
        write(6,*) "WARNING in FrontTracking_init(): No blocks found." 
        return
    end if

    if(blockCount .ne. 1) then
        write(6,*) "Error in FrontTracking_init(): Too many blocks. (",&
                    blockCount, ") Only single block supported!!"
        call exit(-1)
    end if

    thisGridBlock = blockList(1) !FIXME: hardwired assumption of UG!

#ifdef COMP_VAR
    ! Get interior data pointer (for velocities)
    call Grid_getBlkPtr(thisGridBlock,solnData)

    ! Initialize the components to be -1
    if (.not. fr_restart) then
        solnData(COMP_VAR, :, :, :) = -1
    end if
#endif

    ! Get the domain information
    call Grid_getBlkBoundBox(thisGridBlock,boundBox)

    ! Grid spacing information, interior zone and buffer.
    call Grid_getBlkIndexLimits(thisGridBlock,blkLimits,blkLimitsGC)
  
    !Interior
    iisize=blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
    ijsize=blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
    iksize=blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

    !With buffer
    isize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
    jsize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
    ksize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

    ! Set the dimension
    if(ksize .eq. 1) then
        numDim=2
    else
        numDim=3
    end if                      !FIXME: convoluted way to get numDim... just use NDIM !

    ! FIXME: assumed symmetry, simply change each case: blkLimitsGC - blkLimits 
    ibuf = (isize - iisize)/2
    jbuf = (jsize - ijsize)/2
    kbuf = (ksize - iksize)/2
 
    ! Allocate storage arrays for the propagation velocity 
    ! (needed during init to compute front dt)
    allocate(xvel( (isize+1)*(jsize+1)*(ksize+1)))
    allocate(yvel( (isize+1)*(jsize+1)*(ksize+1)))
    allocate(zvel( (isize+1)*(jsize+1)*(ksize+1)))
    !allocate(FrontTracking_compGrid((isize+1)*(jsize+1)*(ksize+1)))
 
    ! CONVERT THE MULTI-DIMENSIONAL ARRAY TO A LINEAR ARRAY
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

    if(.NOT. fr_useFrontTracking) then
        return
    end if
    if (fr_restart) then
        call RuntimeParameters_get('checkPointFileNumber', restart_num)
        write(restart_num_str, "(I0.4)") restart_num
        restart_filename = "ft_chk"//restart_num_str
        rest_len = LEN(TRIM(restart_filename))
    endif

    ! Read the boundary conditions
    call Grid_getBlkBC(thisGridBlock, bdry1, bdry)
    ! JAM: If you have a periodic bdry, it is stored in bdry, not bdry1 so
    ! we need to use bdry in FTAPI call below  
 
    ! FIXME: Hardcoded for UG and standard FLASH procedure of gridsize being
    ! iProcs*nxb , jProcs*nyb , kProcs*nzb
    ! JAM : When running in parallel, isize, jsize and ksize are for the local
    ! processor, thus we should

    call Driver_getNumProcs(AXIS_COMM, procGrid(1), IAXIS) !FIXME: hardwired assumption of UG?
    call Driver_getNumProcs(AXIS_COMM, procGrid(2), JAXIS)
    call Driver_getNumProcs(AXIS_COMM, procGrid(3), KAXIS)

    ! FIXME: Hardcoded for UG and standard FLASH procedure of gridsize being
    ! iProcs*nxb , jProcs*nyb , kProcs*nzb
    ! JAM : When running in parallel, isize, jsize and ksize are for the local
    ! processor, thus we should
    
    call Grid_getGlobalIndexLimits(indexLims)
    call Grid_getGeometry(geometry)
    call FTAPI_init(numDim, procGrid,         &
                boundbox(1,1), boundbox(2,1), &
                boundbox(1,2), boundbox(2,2), &
                boundbox(1,3), boundbox(2,3), &
                indexLims(1), indexLims(2),   &
                indexLims(3),                 &
                ibuf, jbuf, kbuf,             &
                bdry(1,1), bdry(2,1),         &
                bdry(1,2), bdry(2,2),         &
                bdry(1,3), bdry(2,3),         &
                FrontTracking_levelFunc,%VAL(0),   &  !level_func -> unused?
                FrontTracking_getVel,%VAL(0),              &  !vel_func -> unused?
                geometry, fr_restart,         &
                TRIM(restart_filename), rest_len)
     
    ! JAM - Adding a call to Propagate front with a dt of 0.0

    call FrontTracking_propagate(0.0)
   
    ! JAM - We need to initialize the components after adding subcycling to the propagate step 
    call FrontTracking_updateComp()

    deallocate(xvel)
    deallocate(yvel)
    deallocate(zvel)

    call Grid_fillGuardCells(CENTER_FACES, ALLDIR)

    call Grid_releaseBlkPtr(thisGridBlock, solnData)

end subroutine FrontTracking_init

