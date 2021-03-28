#define aindex(i,j,k,isize,jsize,ksize) (((((k)-1)*(jsize)+((j)-1))*(isize))+((i)-1) + 1 )
subroutine FrontTracking_updateComp()
    use Grid_interface, only: Grid_getBlkPtr, Grid_getBlkIndexLimits,&
        Grid_getLocalNumBlks, Grid_getListOfBlocks, Grid_getCellCoords, &
        Grid_getBlkBoundBox, Grid_releaseBlkPtr, Grid_getGeometry 
    use FrontTracking_data, ONLY: fr_useFrontTracking, FrontTracking_compGrid 
    use Grid_data, only : gr_delta
        
implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

        integer,dimension(MAXBLOCKS) :: blockList
        integer :: blockCount
        integer :: thisGridBlock
        integer :: ind
        integer, dimension(3) :: icoord
        integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
        integer :: oldComp, newComp
        REAL , dimension(2,MDIM) :: boundBox
        real*8, dimension(3) :: coords
        INTEGER :: i,j,k
        integer :: isize ,jsize,ksize
        real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
        real :: dens, energy, velx, vely
        real, pointer :: solnData(:,:,:,:)

        call Grid_getListOfBlocks(LEAF,blockList,blockCount)
        
        thisGridBlock = blockList(1)

        call Grid_getBlkBoundBox(thisGridBlock,boundBox)
        call Grid_getBlkIndexLimits(thisGridBlock,blkLimits,blkLimitsGC)

        call Grid_getBlkPtr(thisGridBlock,solnData)
          
        isize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jsize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        ksize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        allocate(xCoord(isize))
        allocate(yCoord(jsize))
        allocate(zCoord(ksize))
        
        call Grid_getCellCoords(IAXIS, blockList(1), CENTER, .false., &
        xCoord, isize)
        call Grid_getCellCoords(JAXIS, blockList(1), CENTER, .false., &
        yCoord, jsize)
        call Grid_getCellCoords(KAXIS, blockList(1), CENTER, .false., &
        zCoord, ksize)
        
        !write(*,*) zCoord, gr_delta(3), ksize
        !call fr_print(FrontTracking_compGrid)
        
       ! do k = 1,1 !blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
               
               ind=aindex(i,j,1,isize,jsize,1)
               oldComp = FrontTracking_compGrid(ind)
               coords(1) = xCoord(i)
               coords(2) = yCoord(j)
               coords(3) = -0.25 

               icoord(1) = i
               icoord(2) = j
               icoord(3) = 1
               
               !call FTAPI_getComponent(coords, newComp)

               if(newComp .eq. oldComp) then
                    go to 80 !write(*,*) "in if", newComp, oldComp, coords(1), coords(2)
               else if(oldComp .eq. -1) then
                    FrontTracking_compGrid(ind) = newComp
               else
                    FrontTracking_compGrid(ind) = newComp
                    call FrontTracking_getGhostState(icoord, newComp, dens, energy, velx, vely)
                    solnData(VELX_VAR,i,j,1) = velx
                    solnData(VELY_VAR,i,j,1) = vely
                    solnData(EINT_VAR,i,j,1) = energy
                    solnData(DENS_VAR,i,j,1) = dens
                    write(*,*) "ghost"
               80 end if

        end do
        end do
        !end do

        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)

end subroutine
