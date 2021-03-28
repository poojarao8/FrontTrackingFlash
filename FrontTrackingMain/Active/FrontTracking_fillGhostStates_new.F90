subroutine FrontTracking_fillGhostStates()

#include "Flash.h"
#include "constants.h"

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                        Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none

    real, pointer                       ::  solnData(:,:,:,:)

    integer                             ::  blkID
    integer                             ::  icell, jcell, kcell
    integer, dimension(2, MDIM)         ::  blkLimits, blkLimitsGC
    real                                ::  dens_sum, pres_sum
    real                                ::  game_sum, gamc_sum
    real                                ::  velx_sum, vely_sum, velz_sum
    integer                             ::  c, num, i, j, k, counter

    ! JAM -- Added for FTAPI   FIXME: I would rather have this passed in than
    !                                 hardcoded, it is not compatible with 
    !                                 a lot of the options in FLASH
    blkID = 1

    call Grid_getBlkPtr(blkID,solnData)
    call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC)

    ! Initialize Ghost States to 0.0 before setting
    do icell = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        do jcell = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
            do kcell = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                solnData(G0CGAM_VAR:G1ZVEL_VAR, icell, jcell, kcell) = 0.0
            end do
        end do
    end do

    ! FIXME: Hardcoded to two physical components -- will have to do for now
    ! FIXME: SUPER UGLY! Need to be using vectors and what not...

    ! Pass 1 through the domain, here we are going to use the real state of the
    ! opp component to interp a ghost state in the cell. 
    ! This gives us coverage one cells from the front, with the second pass
    ! below giving us the coverage two cells from the front which should be all
    ! that is needed for all simulations.
    do icell = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        do jcell = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
            do kcell = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                c = solnData(COMP_VAR, icell, jcell, kcell)
                
                ! During the first pass, we fill the ghost states for this
                ! component to match its real state.
                ! If component of cell is 1 then G0 should match the cell
                ! states, otherwise G1 should match the cell states
                if (c == 1) then 
                    solnData(G0DENS_VAR, icell, jcell, kcell) = &
                        solnData(DENS_VAR, icell, jcell, kcell)
                    solnData(G0PRES_VAR, icell, jcell, kcell) = &
                        solnData(PRES_VAR, icell, jcell, kcell)
                    solnData(G0EGAM_VAR, icell, jcell, kcell) = &
                        solnData(GAME_VAR, icell, jcell, kcell)
                    solnData(G0CGAM_VAR, icell, jcell, kcell) = &
                        solnData(GAMC_VAR, icell, jcell, kcell)
                    solnData(G0XVEL_VAR, icell, jcell, kcell) = &
                        solnData(VELX_VAR, icell, jcell, kcell)
                    solnData(G0YVEL_VAR, icell, jcell, kcell) = &
                        solnData(VELY_VAR, icell, jcell, kcell)
                    solnData(G0ZVEL_VAR, icell, jcell, kcell) = &
                        solnData(VELZ_VAR, icell, jcell, kcell)

                    solnData(G1DENS_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1PRES_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1EGAM_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1CGAM_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1XVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1YVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1ZVEL_VAR, icell, jcell, kcell) = 0.0
                else
                    solnData(G1DENS_VAR, icell, jcell, kcell) = &
                        solnData(DENS_VAR, icell, jcell, kcell)
                    solnData(G1PRES_VAR, icell, jcell, kcell) = &
                        solnData(PRES_VAR, icell, jcell, kcell)
                    solnData(G1EGAM_VAR, icell, jcell, kcell) = &
                        solnData(GAME_VAR, icell, jcell, kcell)
                    solnData(G1CGAM_VAR, icell, jcell, kcell) = &
                        solnData(GAMC_VAR, icell, jcell, kcell)
                    solnData(G1XVEL_VAR, icell, jcell, kcell) = &
                        solnData(VELX_VAR, icell, jcell, kcell)
                    solnData(G1YVEL_VAR, icell, jcell, kcell) = &
                        solnData(VELY_VAR, icell, jcell, kcell)
                    solnData(G1ZVEL_VAR, icell, jcell, kcell) = &
                        solnData(VELZ_VAR, icell, jcell, kcell)

                    solnData(G0DENS_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0PRES_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0XVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0YVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0ZVEL_VAR, icell, jcell, kcell) = 0.0
                end if

                ! If component of cell is 1 then G1 should be an average of
                ! neighbors with component 2 and vice versa
                dens_sum = 0.0
                pres_sum = 0.0
                game_sum = 0.0
                gamc_sum = 0.0
                velx_sum = 0.0
                vely_sum = 0.0
                velz_sum = 0.0
                num = 0
                counter = 0
                if (c == 1) then 
                    ! Loop over all neighbors
                    ! TODO: Verify FORTRAN loop acts as a strict <= upper bound
                    do i = icell - 1, icell + 1     
                        do j = jcell - 1, jcell + 1
                            do k = kcell - 1, kcell + 1
                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
                                    cycle
                                end if
                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
                                    cycle
                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
                                    cycle
                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
                                    cycle 
                                end if
                                ! We are only going to use cells of the opposite
                                ! component to fill the ghost state, i.e.
                                ! interpolated from real states, we will do a
                                ! second pass that interpolates it from ghost
                                ! states
                                if (solnData(COMP_VAR,i,j,k) .eq. 2) then
                                    dens_sum = dens_sum + solnData(DENS_VAR,i,j,k)
                                    pres_sum = pres_sum + solnData(PRES_VAR,i,j,k)
                                    game_sum = game_sum + solnData(GAME_VAR,i,j,k)
                                    gamc_sum = gamc_sum + solnData(GAMC_VAR,i,j,k)
                                    velx_sum = velx_sum + solnData(VELX_VAR,i,j,k)
                                    vely_sum = vely_sum + solnData(VELY_VAR,i,j,k)
                                    velz_sum = velz_sum + solnData(VELZ_VAR,i,j,k)
                                    num = num + 1
                                end if
                            end do
                        end do
                    end do
                    if (num .gt. 0) then
                        solnData(G1DENS_VAR, icell, jcell, kcell) = dens_sum/num
                        solnData(G1PRES_VAR, icell, jcell, kcell) = pres_sum/num
                        solnData(G1EGAM_VAR, icell, jcell, kcell) = game_sum/num
                        solnData(G1CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
                        solnData(G1XVEL_VAR, icell, jcell, kcell) = velx_sum/num
                        solnData(G1YVEL_VAR, icell, jcell, kcell) = vely_sum/num
                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
                    else
                        solnData(G1DENS_VAR, icell, jcell, kcell) = 0.0  
                        solnData(G1PRES_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1EGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1CGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1XVEL_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1YVEL_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = 0.0
                    end if
                else
                    ! Loop over all neighbors
                    ! TODO: Verify FORTRAN loop acts as a  <= upper bound
                    do i = icell - 1, icell + 1     
                        do j = jcell - 1, jcell + 1
                            do k = kcell - 1, kcell + 1
                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
                                    cycle 
                                end if
                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
                                    cycle
                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
                                    cycle
                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
                                    cycle
                                end if
`                               ! See above note
                                if (solnData(COMP_VAR,i,j,k) .eq. 1) then
                                    dens_sum = dens_sum + solnData(DENS_VAR,i,j,k)
                                    pres_sum = pres_sum + solnData(PRES_VAR,i,j,k)
                                    game_sum = game_sum + solnData(GAME_VAR,i,j,k)
                                    gamc_sum = gamc_sum + solnData(GAMC_VAR,i,j,k)
                                    velx_sum = velx_sum + solnData(VELX_VAR,i,j,k)
                                    vely_sum = vely_sum + solnData(VELY_VAR,i,j,k)
                                    velz_sum = velz_sum + solnData(VELZ_VAR,i,j,k)
                                    num = num + 1
                                end if
                            end do
                        end do
                    end do
                    if (num .gt. 0) then
                        solnData(G0DENS_VAR, icell, jcell, kcell) = dens_sum/num
                        solnData(G0PRES_VAR, icell, jcell, kcell) = pres_sum/num
                        solnData(G0EGAM_VAR, icell, jcell, kcell) = game_sum/num
                        solnData(G0CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
                        solnData(G0XVEL_VAR, icell, jcell, kcell) = velx_sum/num
                        solnData(G0YVEL_VAR, icell, jcell, kcell) = vely_sum/num
                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
                    else
                        solnData(G0DENS_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G0PRES_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0EGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0CGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0XVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0YVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = 0.0
                    end if 
                end if       
            end do
        end do
    end do     

    ! Pass 2 through the domain, here we are going to use the filled in ghost
    ! states to interpolate a ghost state in a new cell, this will only be used
    ! if none of the neighbor cells have a real state of the other component.
    ! This gives us coverage two cells from the front, which should be enough
    ! for all simulations.
    do icell = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        do jcell = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
            do kcell = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                c = solnData(COMP_VAR, icell, jcell, kcell)
                dens_sum = 0.0
                pres_sum = 0.0
                game_sum = 0.0
                gamc_sum = 0.0
                velx_sum = 0.0
                vely_sum = 0.0
                velz_sum = 0.0
                num = 0
                counter = 0
                if (c == 1) then
                    ! Loop over all neighbors and see if any were of the
                    ! opposite component, if so then they were filled in during
                    ! pass 1 and can be ignored here
                    do i = icell - 1, icell + 1
                        do j = jcell - 1, jcell + 1
                            do k = kcell - 1, kcell + 1
                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
                                    cycle
                                end if
                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
                                    cycle
                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
                                    cycle
                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
                                    cycle
                                end if
                                if (solnData(COMP_VAR,i,j,k) .eq. 2)
                                    counter = counter + 1
                                end if
                            end do
                        end do
                    end do

                    if (counter > 0) then
                        cycle
                    end if
                    ! Now we loop over the neighbors and see if they have a
                    ! non-zero ghost state, if so we use that to interpolate
                    ! this cells ghost state
                    do i = icell - 1, icell + 1
                        do j = jcell - 1, jcell + 1
                            do k = kcell - 1, kcell + 1
                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
                                    cycle
                                end if
                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
                                    cycle
                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
                                    cycle
                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
                                    cycle
                                end if
                                if (solnData(G1DENS_VAR,i,j,k) > 0.0) then
                                    dens_sum = dens_sum + solnData(G1DENS_VAR,i,j,k)
                                    pres_sum = pres_sum + solnData(G1PRES_VAR,i,j,k)
                                    game_sum = game_sum + solnData(G1EGAM_VAR,i,j,k)
                                    gamc_sum = gamc_sum + solnData(G1CGAM_VAR,i,j,k)
                                    velx_sum = velx_sum + solnData(G1XVEL_VAR,i,j,k)
                                    vely_sum = vely_sum + solnData(G1YVEL_VAR,i,j,k)
                                    velz_sum = velz_sum + solnData(G1ZVEL_VAR,i,j,k)
                                    num = num + 1
                                end if
                            end do
                        end do
                    end do
                    if (num .gt. 0) then
                        solnData(G1DENS_VAR, icell, jcell, kcell) = dens_sum/num
                        solnData(G1PRES_VAR, icell, jcell, kcell) = pres_sum/num
                        solnData(G1EGAM_VAR, icell, jcell, kcell) = game_sum/num
                        solnData(G1CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
                        solnData(G1XVEL_VAR, icell, jcell, kcell) = velx_sum/num
                        solnData(G1YVEL_VAR, icell, jcell, kcell) = vely_sum/num
                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
                    else
                        solnData(G1DENS_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1PRES_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1EGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1CGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1XVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1YVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = 0.0
                    end if
                else
                    ! Loop over all neighbors and see if any were of the
                    ! opposite component, if so then they were filled in during
                    ! pass 1 and can be ignored here
                    do i = icell - 1, icell + 1
                        do j = jcell - 1, jcell + 1
                            do k = kcell - 1, kcell + 1
                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
                                    cycle
                                end if
                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
                                    cycle
                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
                                    cycle
                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
                                    cycle
                                end if
                                if (solnData(COMP_VAR,i,j,k) .eq. 1)
                                    counter = counter + 1
                                end if
                            end do
                        end do
                    end do

                    if (counter > 0) then
                        cycle
                    end if
                    ! Now we loop over the neighbors and see if they have a
                    ! non-zero ghost state, if so we use that to interpolate
                    ! this cells ghost state
                    do i = icell - 1, icell + 1
                        do j = jcell - 1, jcell + 1
                            do k = kcell - 1, kcell + 1
                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
                                    cycle
                                end if
                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
                                    cycle
                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
                                    cycle
                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
                                    cycle
                                end if
`                               ! See above note
                                if (solnData(G0DENS_VAR,i,j,k) > 0.0) then
                                    dens_sum = dens_sum + solnData(G0DENS_VAR,i,j,k)
                                    pres_sum = pres_sum + solnData(G0PRES_VAR,i,j,k)
                                    game_sum = game_sum + solnData(G0EGAM_VAR,i,j,k)
                                    gamc_sum = gamc_sum + solnData(G0CGAM_VAR,i,j,k)
                                    velx_sum = velx_sum + solnData(G0XVEL_VAR,i,j,k)
                                    vely_sum = vely_sum + solnData(G0YVEL_VAR,i,j,k)
                                    velz_sum = velz_sum + solnData(G0ZVEL_VAR,i,j,k)
                                    num = num + 1
                                end if
                            end do
                        end do
                    end do
                    if (num .gt. 0) then
                        solnData(G0DENS_VAR, icell, jcell, kcell) = dens_sum/num
                        solnData(G0PRES_VAR, icell, jcell, kcell) = pres_sum/num
                        solnData(G0EGAM_VAR, icell, jcell, kcell) = game_sum/num
                        solnData(G0CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
                        solnData(G0XVEL_VAR, icell, jcell, kcell) = velx_sum/num
                        solnData(G0YVEL_VAR, icell, jcell, kcell) = vely_sum/num
                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
                    else
                        solnData(G0DENS_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0PRES_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0EGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0CGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0XVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0YVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = 0.0
                    end if
                end if
            end do
        end do
    end do

    call Grid_releaseBlkPtr(blkID,solnData)

    return

end subroutine FrontTracking_fillGhostStates
