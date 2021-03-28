subroutine FrontTracking_fillGhostStates(xnGhostVal)

#include "Flash.h"
#include "constants.h"

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                        Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
    use Driver_interface, ONLY : Driver_abortFlash

    implicit none

    ! JAM - Added to have ghost states for mass scalars stored
    ! 2 is for the number of components
    real, INTENT(INOUT), DIMENSION(NSPECIES+NMASS_SCALARS,   & 
                                   GRID_ILO_GC:GRID_IHI_GC,  &
                                   GRID_JLO_GC:GRID_JHI_GC,  &
                                   GRID_KLO_GC:GRID_KHI_GC,  &
                                   2)       :: xnGhostVal

    real, pointer                           ::  solnData(:,:,:,:)

    integer                                 ::  blkID
    integer                                 ::  icell, jcell, kcell, xx, step
    integer, dimension(2, MDIM)             ::  blkLimits, blkLimitsGC
    real                                    ::  dens_sum, pres_sum
    real                                    ::  game_sum, gamc_sum
    real                                    ::  velx_sum, vely_sum, velz_sum
    real, dimension(NSPECIES+NMASS_SCALARS) ::  xn_sum
    integer                                 ::  c, num, i, j, k, counter, numXN

#ifdef FLASH_3T
    real                                    ::  pion_sum, pele_sum, prad_sum
#endif



    ! JAM -- Added for FTAPI   FIXME: I would rather have this passed in than
    !                                 hardcoded, it is not compatible with 
    !                                 a lot of the options in FLASH
    blkID = 1

    numXN = NSPECIES+NMASS_SCALARS

    call Grid_getBlkPtr(blkID,solnData)
    call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC)

    ! Initialize Ghost States to 0.0 before setting
    do icell = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        do jcell = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
            do kcell = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                solnData(G0CGAM_VAR:G1ZVEL_VAR, icell, jcell, kcell) = 0.0
                xnGhostVal(:, icell, jcell, kcell, :) = 0.0
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
#ifdef FLASH_3T
                    solnData(G0IONP_VAR, icell, jcell, kcell) = &
                        solnData(PION_VAR, icell, jcell, kcell)
                    solnData(G0ELEP_VAR, icell, jcell, kcell) = &
                        solnData(PELE_VAR, icell, jcell, kcell)
                    solnData(G0RADP_VAR, icell, jcell, kcell) = &
                        solnData(PRAD_VAR, icell, jcell, kcell)
#endif
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
                    do xx = 1, numXN
                        xnGhostVal(xx, icell, jcell, kcell, 1) = &
                            solnData(SPECIES_BEGIN-1+xx, icell, jcell, kcell)
                    end do 

                    solnData(G1DENS_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1PRES_VAR, icell, jcell, kcell) = 0.0
#ifdef FLASH_3T     
                    solnData(G1IONP_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1ELEP_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1RADP_VAR, icell, jcell, kcell) = 0.0
#endif
                    solnData(G1EGAM_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1CGAM_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1XVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1YVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G1ZVEL_VAR, icell, jcell, kcell) = 0.0
                    do xx = 1, numXN
                        xnGhostVal(xx, icell, jcell, kcell, 2) = 0.0
                    end do
                else
                    solnData(G1DENS_VAR, icell, jcell, kcell) = &
                        solnData(DENS_VAR, icell, jcell, kcell)
                    solnData(G1PRES_VAR, icell, jcell, kcell) = &
                        solnData(PRES_VAR, icell, jcell, kcell)
#ifdef FLASH_3T     
                    solnData(G1IONP_VAR, icell, jcell, kcell) = &
                        solnData(PION_VAR, icell, jcell, kcell)
                    solnData(G1ELEP_VAR, icell, jcell, kcell) = &
                        solnData(PELE_VAR, icell, jcell, kcell)
                    solnData(G1RADP_VAR, icell, jcell, kcell) = &
                        solnData(PRAD_VAR, icell, jcell, kcell)
#endif
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
                    do xx = 1, numXN
                        xnGhostVal(xx, icell, jcell, kcell, 2) = &
                            solnData(SPECIES_BEGIN-1+xx, icell, jcell, kcell)
                    end do

                    solnData(G0DENS_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0PRES_VAR, icell, jcell, kcell) = 0.0
#ifdef FLASH_3T     
                    solnData(G0IONP_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0ELEP_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0RADP_VAR, icell, jcell, kcell) = 0.0
#endif
                    solnData(G0XVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0YVEL_VAR, icell, jcell, kcell) = 0.0
                    solnData(G0ZVEL_VAR, icell, jcell, kcell) = 0.0
                    do xx = 1, numXN
                        xnGhostVal(xx, icell, jcell, kcell, 1) = 0.0
                    end do
                end if

                ! If component of cell is 1 then G1 should be an average of
                ! neighbors with component 2 and vice versa
                dens_sum = 0.0
                pres_sum = 0.0
#ifdef FLASH_3T
                pion_sum = 0.0
                pele_sum = 0.0
                prad_sum = 0.0
#endif
                game_sum = 0.0
                gamc_sum = 0.0
                velx_sum = 0.0
                vely_sum = 0.0
                velz_sum = 0.0
                xn_sum(:) = 0.0
                num = 0
                counter = 0
                if (c == 1) then 
                    ! Loop over all neighbors
                    ! TODO: Verify FORTRAN loop acts as a strict <= upper bound
                    do step = 1, 5
                        do i = icell - step, icell + step     
                            do j = jcell - step, jcell + step
                                do k = kcell - step, kcell + step
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
                                    ! This checks to see if we are within the distance
                                    ! of one of the previous checks, so that we don't 
                                    ! check something we have already eliminated
                                    if (abs(i - icell) .le. (step - 1) .and. &
                                        abs(j - jcell) .le. (step - 1) .and. &
                                        abs(k - kcell) .le. (step - 1)) then
                                        cycle
                                    end if
                                    ! We are only going to use cells of the opposite
                                    ! component to fill the ghost state, i.e.
                                    ! interpolated from real states
                                    if (solnData(COMP_VAR,i,j,k) .eq. 2) then
                                        dens_sum = dens_sum + solnData(DENS_VAR,i,j,k)
                                        pres_sum = pres_sum + solnData(PRES_VAR,i,j,k)
#ifdef FLASH_3T     
                                        pion_sum = pion_sum + solnData(PION_VAR,i,j,k)
                                        pele_sum = pele_sum + solnData(PELE_VAR,i,j,k)
                                        prad_sum = prad_sum + solnData(PRAD_VAR,i,j,k)
#endif
                                        game_sum = game_sum + solnData(GAME_VAR,i,j,k)
                                        gamc_sum = gamc_sum + solnData(GAMC_VAR,i,j,k)
                                        velx_sum = velx_sum + solnData(VELX_VAR,i,j,k)
                                        vely_sum = vely_sum + solnData(VELY_VAR,i,j,k)
                                        velz_sum = velz_sum + solnData(VELZ_VAR,i,j,k)
                                        do xx = 1, numXN
                                            xn_sum(xx) = xn_sum(xx) + &
                                                    solnData(SPECIES_BEGIN-1+xx,i,j,k)
                                        end do
                                        num = num + 1
                                    end if
                                end do
                            end do
                        end do
                        if (num .gt. 0) then
                            solnData(G1DENS_VAR, icell, jcell, kcell) = dens_sum/num
                            solnData(G1PRES_VAR, icell, jcell, kcell) = pres_sum/num
#ifdef FLASH_3T     
                            solnData(G1IONP_VAR, icell, jcell, kcell) = pion_sum/num
                            solnData(G1ELEP_VAR, icell, jcell, kcell) = pele_sum/num
                            solnData(G1RADP_VAR, icell, jcell, kcell) = prad_sum/num
#endif
                            solnData(G1EGAM_VAR, icell, jcell, kcell) = game_sum/num
                            solnData(G1CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
                            solnData(G1XVEL_VAR, icell, jcell, kcell) = velx_sum/num
                            solnData(G1YVEL_VAR, icell, jcell, kcell) = vely_sum/num
                            solnData(G1ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
                            do xx = 1, numXN
                                xnGhostVal(xx, icell, jcell, kcell, 2) = xn_sum(xx)/num
                            end do
                            exit
                        end if
                    end do
                    ! If we didnt hit any opposite components going up to a box
                    ! 5 cells in each direction, assume we are too far from the
                    ! front to need ghost states, set to 0 and move to the next
                    ! cell
                    if (num .eq. 0) then
                        solnData(G1DENS_VAR, icell, jcell, kcell) = 0.0  
                        solnData(G1PRES_VAR, icell, jcell, kcell) = 0.0 
#ifdef FLASH_3T     
                        solnData(G1IONP_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1ELEP_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1RADP_VAR, icell, jcell, kcell) = 0.0 
#endif
                        solnData(G1EGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1CGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G1XVEL_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1YVEL_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = 0.0
                        do xx = 1, numXN
                            xnGhostVal(xx, icell, jcell, kcell, 2) = 0.0
                        end do
                    end if
                else
                    ! Loop over all neighbors
                    ! FIXME: This 5 is arbitrary, simply how far you want to
                    ! propagate the ghost cell data away from the front, it
                    ! should be set as an input option
                    do step = 1, 5
                        do i = icell - step, icell + step     
                            do j = jcell - step, jcell + step
                                do k = kcell - step, kcell + step
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
                                    ! This checks to see if we are within the
                                    ! distance
                                    ! of one of the previous checks, so that we
                                    ! don't 
                                    ! check something we have already eliminated
                                    if (abs(i - icell) .le. (step - 1) .and. &
                                        abs(j - jcell) .le. (step - 1) .and. &
                                        abs(k - kcell) .le. (step - 1)) then
                                        cycle
                                    end if

                                    ! See above note
                                    if (solnData(COMP_VAR,i,j,k) .eq. 1) then
                                        dens_sum = dens_sum + solnData(DENS_VAR,i,j,k)
                                        pres_sum = pres_sum + solnData(PRES_VAR,i,j,k)
#ifdef FLASH_3T     
                                        pion_sum = pion_sum + solnData(PION_VAR,i,j,k)
                                        pele_sum = pele_sum + solnData(PELE_VAR,i,j,k)
                                        prad_sum = prad_sum + solnData(PRAD_VAR,i,j,k)
#endif
                                        game_sum = game_sum + solnData(GAME_VAR,i,j,k)
                                        gamc_sum = gamc_sum + solnData(GAMC_VAR,i,j,k)
                                        velx_sum = velx_sum + solnData(VELX_VAR,i,j,k)
                                        vely_sum = vely_sum + solnData(VELY_VAR,i,j,k)
                                        velz_sum = velz_sum + solnData(VELZ_VAR,i,j,k)
                                        do xx = 1, numXN
                                            xn_sum(xx) = xn_sum(xx) + &
                                                    solnData(SPECIES_BEGIN-1+xx,i,j,k)
                                        end do
                                        num = num + 1
                                    end if
                                end do
                            end do
                        end do
                        if (num .gt. 0) then
                            solnData(G0DENS_VAR, icell, jcell, kcell) = dens_sum/num
                            solnData(G0PRES_VAR, icell, jcell, kcell) = pres_sum/num
#ifdef FLASH_3T     
                            solnData(G0IONP_VAR, icell, jcell, kcell) = pion_sum/num
                            solnData(G0ELEP_VAR, icell, jcell, kcell) = pele_sum/num
                            solnData(G0RADP_VAR, icell, jcell, kcell) = prad_sum/num
#endif
                            solnData(G0EGAM_VAR, icell, jcell, kcell) = game_sum/num
                            solnData(G0CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
                            solnData(G0XVEL_VAR, icell, jcell, kcell) = velx_sum/num
                            solnData(G0YVEL_VAR, icell, jcell, kcell) = vely_sum/num
                            solnData(G0ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
                            do xx = 1, numXN
                                xnGhostVal(xx, icell, jcell, kcell, 1) = xn_sum(xx)/num
                            end do
                            exit
                        end if
                    end do
                    ! If we didnt hit any opposite components going up to a box
                    ! 5 cells in each direction, assume we are too far from the
                    ! front to need ghost states, set to 0 and move to the next
                    ! cell
                    if (num .eq. 0) then
                        solnData(G0DENS_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G0PRES_VAR, icell, jcell, kcell) = 0.0
#ifdef FLASH_3T     
                        solnData(G0IONP_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G0ELEP_VAR, icell, jcell, kcell) = 0.0 
                        solnData(G0RADP_VAR, icell, jcell, kcell) = 0.0 
#endif
                        solnData(G0EGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0CGAM_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0XVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0YVEL_VAR, icell, jcell, kcell) = 0.0
                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = 0.0
                        do xx = 1, numXN
                            xnGhostVal(xx, icell, jcell, kcell, 1) = 0.0
                        end do
                    end if 
                end if       
            end do
        end do
    end do     

    ! JAM : We are no longer doing a second pass and filling from ghost states

    ! Pass 2 through the domain, here we are going to use the filled in ghost
    ! states to interpolate a ghost state in a new cell, this will only be used
    ! if none of the neighbor cells have a real state of the other component.
    ! This gives us coverage two cells from the front, which should be enough
    ! for all simulations.
!    do icell = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
!        do jcell = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
!            do kcell = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
!                c = solnData(COMP_VAR, icell, jcell, kcell)
!                dens_sum = 0.0
!                pres_sum = 0.0
!#ifdef FLASH_3T
!                pion_sum = 0.0
!                pele_sum = 0.0
!                prad_sum = 0.0
!#endif
!                game_sum = 0.0
!                gamc_sum = 0.0
!                velx_sum = 0.0
!                vely_sum = 0.0
!                velz_sum = 0.0
!                xn_sum(:) = 0.0
!                num = 0
!                counter = 0
!                if (c == 1) then
!                    ! Loop over all neighbors and see if any were of the
!                    ! opposite component, if so then they were filled in during
!                    ! pass 1 and can be ignored here
!                    do i = icell - 1, icell + 1
!                        do j = jcell - 1, jcell + 1
!                            do k = kcell - 1, kcell + 1
!                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
!                                    cycle
!                                end if
!                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
!                                    cycle
!                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
!                                    cycle
!                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
!                                    cycle
!                                end if
!                                if (solnData(COMP_VAR,i,j,k) .eq. 2) then
!                                    counter = counter + 1
!                                end if
!                            end do
!                        end do
!                    end do
!
!                    if (counter > 0) then
!                        cycle
!                    end if
!
!                    ! Now we loop over the neighbors and see if they have a
!                    ! non-zero ghost state, if so we use that to interpolate
!                    ! this cells ghost state
!                    do i = icell - 1, icell + 1
!                        do j = jcell - 1, jcell + 1
!                            do k = kcell - 1, kcell + 1
!                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
!                                    cycle
!                                end if
!                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
!                                    cycle
!                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
!                                    cycle
!                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
!                                    cycle
!                                end if
!                                if (solnData(G1DENS_VAR,i,j,k) > 0.0) then
!                                    dens_sum = dens_sum + solnData(G1DENS_VAR,i,j,k)
!                                    pres_sum = pres_sum + solnData(G1PRES_VAR,i,j,k)
!#ifdef FLASH_3T     
!                                    pion_sum = pion_sum + solnData(G1IONP_VAR,i,j,k)
!                                    pele_sum = pele_sum + solnData(G1ELEP_VAR,i,j,k)
!                                    prad_sum = prad_sum + solnData(G1RADP_VAR,i,j,k)
!#endif
!                                    game_sum = game_sum + solnData(G1EGAM_VAR,i,j,k)
!                                    gamc_sum = gamc_sum + solnData(G1CGAM_VAR,i,j,k)
!                                    velx_sum = velx_sum + solnData(G1XVEL_VAR,i,j,k)
!                                    vely_sum = vely_sum + solnData(G1YVEL_VAR,i,j,k)
!                                    velz_sum = velz_sum + solnData(G1ZVEL_VAR,i,j,k)
!                                    do xx = 1, numXN
!                                        xn_sum(xx) = xn_sum(xx) + xnGhostVal(xx,i,j,k,2)
!                                    end do
!                                    num = num + 1
!                                end if
!                            end do
!                        end do
!                    end do
!                    if (num .gt. 0) then
!                        solnData(G1DENS_VAR, icell, jcell, kcell) = dens_sum/num
!                        solnData(G1PRES_VAR, icell, jcell, kcell) = pres_sum/num
!#ifdef FLASH_3T     
!                        solnData(G1IONP_VAR, icell, jcell, kcell) = pion_sum/num
!                        solnData(G1ELEP_VAR, icell, jcell, kcell) = pele_sum/num
!                        solnData(G1RADP_VAR, icell, jcell, kcell) = prad_sum/num
!#endif
!                        solnData(G1EGAM_VAR, icell, jcell, kcell) = game_sum/num
!                        solnData(G1CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
!                        solnData(G1XVEL_VAR, icell, jcell, kcell) = velx_sum/num
!                        solnData(G1YVEL_VAR, icell, jcell, kcell) = vely_sum/num
!                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
!                        do xx = 1, numXN
!                            xnGhostVal(xx, icell, jcell, kcell, 2) = xn_sum(xx)/num
!                        end do
!                    else
!                        solnData(G1DENS_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1PRES_VAR, icell, jcell, kcell) = 0.0
!#ifdef FLASH_3T     
!                        solnData(G1IONP_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1ELEP_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1RADP_VAR, icell, jcell, kcell) = 0.0
!#endif
!                        solnData(G1EGAM_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1CGAM_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1XVEL_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1YVEL_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G1ZVEL_VAR, icell, jcell, kcell) = 0.0
!                        do xx = 1, numXN
!                            xnGhostVal(xx, icell, jcell, kcell, 2) = 0.0
!                        end do
!                    end if
!                else
!                    ! Loop over all neighbors and see if any were of the
!                    ! opposite component, if so then they were filled in during
!                    ! pass 1 and can be ignored here
!                    do i = icell - 1, icell + 1
!                        do j = jcell - 1, jcell + 1
!                            do k = kcell - 1, kcell + 1
!                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
!                                    cycle
!                                end if
!                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
!                                    cycle
!                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
!                                    cycle
!                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
!                                    cycle
!                                end if
!                                if (solnData(COMP_VAR,i,j,k) .eq. 1) then
!                                    counter = counter + 1
!                                end if
!                            end do
!                        end do
!                    end do
!
!                    if (counter > 0) then
!                        cycle
!                    end if
!                    ! Now we loop over the neighbors and see if they have a
!                    ! non-zero ghost state, if so we use that to interpolate
!                    ! this cells ghost state
!                    do i = icell - 1, icell + 1
!                        do j = jcell - 1, jcell + 1
!                            do k = kcell - 1, kcell + 1
!                                if (i .lt. 1 .or. j .lt. 1 .or. k .lt. 1) then
!                                    cycle
!                                end if
!                                if (i .gt. blkLimitsGC(HIGH,IAXIS)) then
!                                    cycle
!                                else if (j .gt. blkLimitsGC(HIGH,JAXIS)) then
!                                    cycle
!                                else if (k .gt. blkLimitsGC(HIGH,KAXIS)) then
!                                    cycle
!                                end if
!                                ! See above note
!                                if (solnData(G0DENS_VAR,i,j,k) > 0.0) then
!                                    dens_sum = dens_sum + solnData(G0DENS_VAR,i,j,k)
!                                    pres_sum = pres_sum + solnData(G0PRES_VAR,i,j,k)
!#ifdef FLASH_3T     
!                                    pion_sum = pion_sum + solnData(G0IONP_VAR,i,j,k)
!                                    pele_sum = pele_sum + solnData(G0ELEP_VAR,i,j,k)
!                                    prad_sum = prad_sum + solnData(G0RADP_VAR,i,j,k)
!#endif
!                                    game_sum = game_sum + solnData(G0EGAM_VAR,i,j,k)
!                                    gamc_sum = gamc_sum + solnData(G0CGAM_VAR,i,j,k)
!                                    velx_sum = velx_sum + solnData(G0XVEL_VAR,i,j,k)
!                                    vely_sum = vely_sum + solnData(G0YVEL_VAR,i,j,k)
!                                    velz_sum = velz_sum + solnData(G0ZVEL_VAR,i,j,k)
!                                    do xx = 1, numXN
!                                        xn_sum(xx) = xn_sum(xx) + xnGhostVal(xx,i,j,k,1)
!                                    end do
!                                    num = num + 1
!                                end if
!                            end do
!                        end do
!                    end do
!                    if (num .gt. 0) then
!                        solnData(G0DENS_VAR, icell, jcell, kcell) = dens_sum/num
!                        solnData(G0PRES_VAR, icell, jcell, kcell) = pres_sum/num
!#ifdef FLASH_3T     
!                        solnData(G0IONP_VAR, icell, jcell, kcell) = pion_sum/num
!                        solnData(G0ELEP_VAR, icell, jcell, kcell) = pele_sum/num
!                        solnData(G0RADP_VAR, icell, jcell, kcell) = prad_sum/num
!#endif
!                        solnData(G0EGAM_VAR, icell, jcell, kcell) = game_sum/num
!                        solnData(G0CGAM_VAR, icell, jcell, kcell) = gamc_sum/num
!                        solnData(G0XVEL_VAR, icell, jcell, kcell) = velx_sum/num
!                        solnData(G0YVEL_VAR, icell, jcell, kcell) = vely_sum/num
!                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = velz_sum/num
!                        do xx = 1, numXN
!                            xnGhostVal(xx, icell, jcell, kcell, 1) = xn_sum(xx)/num
!                        end do
!                    else
!                        solnData(G0DENS_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0PRES_VAR, icell, jcell, kcell) = 0.0
!#ifdef FLASH_3T     
!                        solnData(G0IONP_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0ELEP_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0RADP_VAR, icell, jcell, kcell) = 0.0
!#endif
!                        solnData(G0EGAM_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0CGAM_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0XVEL_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0YVEL_VAR, icell, jcell, kcell) = 0.0
!                        solnData(G0ZVEL_VAR, icell, jcell, kcell) = 0.0
!                        do xx = 1, numXN
!                            xnGhostVal(xx, icell, jcell, kcell, 1) = 0.0
!                        end do
!                    end if
!                end if
!            end do
!        end do
!    end do

    call Grid_releaseBlkPtr(blkID,solnData)

    return

end subroutine FrontTracking_fillGhostStates
