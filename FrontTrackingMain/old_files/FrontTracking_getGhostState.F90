#define aindex(i,j,k,isize,jsize,ksize) (((((k)-1)*(jsize)+((j)-1))*(isize))+((i)-1) + 1 )
subroutine FrontTracking_getGhostState(coord, comp, dens, energy, velx, vely)
    use Grid_interface, only: Grid_getBlkPtr, Grid_getBlkIndexLimits,&
        Grid_getLocalNumBlks, Grid_getListOfBlocks, Grid_getCellCoords, &
        Grid_getBlkBoundBox, Grid_releaseBlkPtr, Grid_getGeometry
    use Grid_data, only : gr_delta
    use FrontTracking_data, ONLY: FrontTracking_compGrid
    use Eos_interface, ONLY : Eos_wrapped
    !use Hydro_data, ONLY : hy_unsplitEosMode
    use Driver_interface, ONLY : Driver_abortFlash

implicit none

#include "constants.h"
#include "Flash.h"

        integer, dimension(3), INTENT(IN) :: coord
        real*8, dimension(3) :: newCoord
        integer :: comp, compCell
        integer :: isize,jsize,ksize
        integer :: displace_i, displace_j, i,j,k
        integer,dimension(MAXBLOCKS) :: blockList
        integer :: blockCount
        integer :: thisGridBlock
        integer :: ind, counter
        integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
        real :: xInd, yInd, zInd
        REAL , dimension(2,MDIM) :: boundBox
        real, pointer :: solnData(:,:,:,:)
        real, INTENT(OUT) :: dens, velx, vely, energy

        dens = 0.
        energy = 0.
        velx = 0.
        vely = 0.
        counter = 0
        
        call Grid_getListOfBlocks(LEAF,blockList,blockCount)
        
        thisGridBlock = blockList(1)

        call Grid_getBlkBoundBox(thisGridBlock,boundBox)
        call Grid_getBlkIndexLimits(thisGridBlock,blkLimits,blkLimitsGC)

        call Grid_getBlkPtr(thisGridBlock,solnData)
          
        isize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jsize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        ksize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
        
        !do displace_k = 0, 0
        do displace_j = -1, 1
        do displace_i = -1, 1
                
           !xInd = ((coord(1) ) / gr_delta(1))
           !yInd = ((coord(2) + 0.50)  / gr_delta(2))
           !zInd = 0 !(coord(3) + 0.25) / gr_delta(3) 
               
           i = coord(1) + displace_i
           j = coord(2) + displace_j
           !k = 1 !zInd + displace_k
           !ind = aindex(i,j,1,isize,jsize,ksize)

           newCoord(1) = (gr_delta(1) * i) - 0.25
           newCoord(2) = (gr_delta(2) * j) - 0.50
           newCoord(3) = -0.25
               
           !call FTAPI_getComponent(newCoord, compCell)

           if((comp .ne. compCell)) then
              !if(displace_i .eq. 0 .AND. displace_j .eq. 0) then
              !  write(*,*) "in if: ", newCoord, i, j
              !  counter = counter + 0
              !else
                dens = solnData(DENS_VAR,i,j,1) + dens
                velx = solnData(VELX_VAR,i,j,1) + velx
                vely = solnData(VELY_VAR,i,j,1) + vely
                energy = solnData(EINT_VAR,i,j,1) + energy
                counter = counter + 1
              !end if
           end if

         end do
         end do
         !end do

         if (counter .eq. 0) then
                call Driver_abortFlash("cell has no neighbor of app. comp.")
         end if

         dens = dens / counter
         velx = velx / counter
         vely = vely / counter
         energy = energy / counter
        
end subroutine
