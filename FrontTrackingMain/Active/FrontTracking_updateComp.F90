subroutine FrontTracking_updateComp(xnGhostVal)

#include "Flash.h"
#include "constants.h"
#include "ftapi.h"

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                        Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
    use Driver_interface, ONLY : Driver_abortFlash
    use Eos_interface, ONLY : Eos_wrapped

    implicit none

    ! JAM - Added to have ghost states for mass scalars stored
    real, INTENT(IN), DIMENSION(NSPECIES+NMASS_SCALARS,   & 
                                GRID_ILO_GC:GRID_IHI_GC,  &
                                GRID_JLO_GC:GRID_JHI_GC,  &
                                GRID_KLO_GC:GRID_KHI_GC,  &
                                2)          :: xnGhostVal

    real, pointer                           ::  solnData(:,:,:,:)
    integer                                 ::  oldComp, newComp

#ifdef FIXEDBLOCKSIZE
    real,DIMENSION(MAXCELLS),target         ::  iCoord, jCoord, kCoord
    ! JAM : We need a temp array to store the updates so that when we get the
    ! ghost cells of components from their neighbors, their neighbors haven't
    ! already been updated to new states and components.  At the end we will
    ! copy over the changes to the solnData array
    real,DIMENSION(UNK_VARS_END, GRID_ILO_GC:GRID_IHI_GC, &
                                 GRID_JLO_GC:GRID_JHI_GC, &
                                 GRID_KLO_GC:GRID_KHI_GC) :: tmpSoln
#else
    real, allocatable,DIMENSION(:),target   ::  iCoord, jCoord, kCoord
    ! See note above
    real, allocatable, DIMENSION(:,:,:,:)   :: tmpSoln
#endif

    integer                                 ::  blkID, isizeGC, jsizeGC, ksizeGC
    integer                                 ::  icell, jcell, kcell, cnt, i, ivar
    integer, dimension(2, MDIM)             ::  blkLimits, blkLimitsGC
    logical                                 ::  gcell = .true.
    real, dimension(MDIM)                   ::  del
    real, dimension(MDIM)                   ::  coords
    real                                    ::  lnormvel, rnormvel, ltanvelx, rtanvelx
    real                                    ::  ltanvely, ltanvelz, rtanvely, rtanvelz
    real, dimension(MDIM)                   ::  normal
    real, dimension(NPROP_VARS)             ::  lstate, rstate, midstate
    real                                    ::  gV, gamc
    ! This will store a list of updated indices so we know which ones to update
    ! at the end.
    integer, dimension (GRID_IHI_GC*GRID_JHI_GC*GRID_KHI_GC, MDIM)   &
                                            ::  updatedIndices
#ifdef FLASH_3T
    real                                    ::  intEIon, intEEle, intERad, eneRat
    real                                    ::  pfracIon, pfracEle, pfracRad, sumEne
#endif

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
    ! See note above
    allocate(tmpSoln(UNK_VARS_END,isizeGC,jsizeGC,ksizeGC))
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
    cnt = 1
    updatedIndices(:,:) = 0
    tmpSoln(:,:,:,:) = 0.0

    ! We only need to loop over the interior cells and update components as we
    ! will update the guard cells after this function and that will overwrite
    ! what is in the boundaries
    do icell = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        coords(IAXIS) = iCoord(icell)
        do jcell = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
            coords(JAXIS) = jCoord(jcell)
            do kcell = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                coords(KAXIS) = kCoord(kcell)
                oldComp = solnData(COMP_VAR, icell, jcell, kcell)
                call FTAPI_getComponent(coords, newComp)
                ! We only want to use the comp array (below if condition), if the
                ! component hasn't changed or is changing from -1, so for others
                ! we will store in tmpSoln until the end so the ghost cell calc
                ! works properly
                tmpSoln(COMP_VAR, icell, jcell, kcell) = newComp
                if ((oldComp .ne. newComp) .and. (newComp .ne. -1) .and. (oldComp .ne. -1) &
                                           .and. (newComp .ne. 0)) then
                    lstate(DENS_VAR) = solnData(DENS_VAR, icell, jcell, kcell) 
                    lstate(PRES_VAR) = solnData(PRES_VAR, icell, jcell, kcell) 
                    lstate(GAME_VAR) = solnData(GAME_VAR, icell, jcell, kcell) 
                    lstate(GAMC_VAR) = solnData(GAMC_VAR, icell, jcell, kcell) 
                    lstate(VELX_VAR) = solnData(VELX_VAR, icell, jcell, kcell) 
                    lstate(VELY_VAR) = solnData(VELY_VAR, icell, jcell, kcell) 
                    lstate(VELZ_VAR) = solnData(VELZ_VAR, icell, jcell, kcell) 
                    ! If my old component was 1 and thus is becoming 2, I need
                    ! to solve a Riemann problem between the current state (comp
                    ! 1) and the new "ghost" state comp 2, thus I use the G1
                    ! array which stores the ghost states for comp 2
                    if (oldComp .eq. 1) then
                        rstate(DENS_VAR) = solnData(G1DENS_VAR, icell, jcell, kcell) 
                        rstate(PRES_VAR) = solnData(G1PRES_VAR, icell, jcell, kcell) 
                        rstate(GAME_VAR) = solnData(G1EGAM_VAR, icell, jcell, kcell) 
                        rstate(GAMC_VAR) = solnData(G1CGAM_VAR, icell, jcell, kcell) 
                        rstate(VELX_VAR) = solnData(G1XVEL_VAR, icell, jcell, kcell) 
                        rstate(VELY_VAR) = solnData(G1YVEL_VAR, icell, jcell, kcell) 
                        rstate(VELZ_VAR) = solnData(G1ZVEL_VAR, icell, jcell, kcell)
#ifdef FLASH_3T
                        pfracIon = solnData(G1IONP_VAR, icell, jcell, kcell) /   &
                                     solnData(G1PRES_VAR, icell, jcell, kcell)
                        pfracEle = solnData(G1ELEP_VAR, icell, jcell, kcell) /   &
                                     solnData(G1PRES_VAR, icell, jcell, kcell)
                        pfracRad = solnData(G1RADP_VAR, icell, jcell, kcell) /   &
                                     solnData(G1PRES_VAR, icell, jcell, kcell)
#endif 
                    else
                        rstate(DENS_VAR) = solnData(G0DENS_VAR, icell, jcell, kcell)
                        rstate(PRES_VAR) = solnData(G0PRES_VAR, icell, jcell, kcell)
                        rstate(GAME_VAR) = solnData(G0EGAM_VAR, icell, jcell, kcell)
                        rstate(GAMC_VAR) = solnData(G0CGAM_VAR, icell, jcell, kcell)
                        rstate(VELX_VAR) = solnData(G0XVEL_VAR, icell, jcell, kcell)
                        rstate(VELY_VAR) = solnData(G0YVEL_VAR, icell, jcell, kcell)
                        rstate(VELZ_VAR) = solnData(G0ZVEL_VAR, icell, jcell, kcell)
#ifdef FLASH_3T
                        pfracIon = solnData(G0IONP_VAR, icell, jcell, kcell) /   &
                                     solnData(G0PRES_VAR, icell, jcell, kcell)
                        pfracEle = solnData(G0ELEP_VAR, icell, jcell, kcell) /   &
                                     solnData(G0PRES_VAR, icell, jcell, kcell)
                        pfracRad = solnData(G0RADP_VAR, icell, jcell, kcell) /   &
                                     solnData(G0PRES_VAR, icell, jcell, kcell)
#endif
                    end if

                    if (rstate(DENS_VAR) .eq. 0.0) then
                        print*,   'ERROR, front moved across more                    &
                                  &than 1 cell in 1 timestep and the code in         &  
                                  &FrontTracking_updateComp.F90 isnt set up for that'
                        call Driver_abortFlash('Error')
                    end if
                                
                    call FTAPI_gridNormal(coords, del, normal)

                    lnormvel = normal(IAXIS)*lstate(VELX_VAR)
                    rnormvel = normal(IAXIS)*rstate(VELX_VAR)
                    if (NDIM .GE. 2) then
                        lnormvel = lnormvel + normal(JAXIS)*lstate(VELY_VAR)
                        rnormvel = rnormvel + normal(JAXIS)*rstate(VELY_VAR)
                    end if
                    if (NDIM .GE. 3) then
                        lnormvel = lnormvel + normal(KAXIS)*lstate(VELZ_VAR)
                        rnormvel = rnormvel + normal(KAXIS)*rstate(VELZ_VAR)    
                    end if        
               
                    ! JAM : XXX : We need to reverse the normals here if the
                    ! front is moving in the "negative" direction so that the
                    ! riemann solver takes the correct left and right states, I
                    ! don't fully understand this but it is what is in cFluid
                    ! and produces correct results so I am okay adding it at
                    ! this point but it should be better understood to make sure
                    ! it doesn't produce extraneous results for a more
                    ! complicated problem
                    if (rnormvel < 0.0) then
                        rnormvel = -rnormvel
                        lnormvel = -lnormvel
                        normal(IAXIS) = -normal(IAXIS)
                        normal(JAXIS) = -normal(JAXIS)
                        normal(KAXIS) = -normal(KAXIS)
                    end if
 
                    ltanvelx = lstate(VELX_VAR) - lnormvel * normal(IAXIS)
                    rtanvelx = rstate(VELX_VAR) - rnormvel * normal(IAXIS)
                    if (NDIM .GE. 2) then
                        ltanvely = lstate(VELY_VAR) - lnormvel * normal(JAXIS)
                        rtanvely = rstate(VELY_VAR) - rnormvel * normal(JAXIS)
                    end if
                    if (NDIM .GE. 3) then
                        ltanvelz = lstate(VELZ_VAR) - lnormvel * normal(KAXIS)
                        rtanvelz = rstate(VELZ_VAR) - rnormvel * normal(KAXIS)
                    end if
                    lstate(VELX_VAR) = lnormvel
                    rstate(VELX_VAR) = rnormvel

                    lstate(VELY_VAR) = 0.0
                    lstate(VELZ_VAR) = 0.0
                    rstate(VELY_VAR) = 0.0
                    rstate(VELZ_VAR) = 0.0
                
                    ! JAM : XXX: I don't quite get this, but for some reason the
                    ! way I have configured the rstate and the lstate, the rstate
                    ! needs to be put in first as the "left" state and the
                    ! lstate needs to go second as the "right" state.  This has
                    ! something to do with the velocities, for example if you
                    ! put them in the other order you get a different wave type
                    ! as the right side may be moving faster than the left or
                    ! vice versa.  This order seems to produce the correct
                    ! result, so I will stick with it, but as I said... I don't
                    ! completely understand why at this point
                    call RiemSolver(rstate, lstate, midstate)
                 
                    ! rtanvel is correct here since it the "ghost" state, i.e
                    ! what the component is becoming
                    solnData(VELX_VAR, icell, jcell, kcell) =             &
                            midstate(VELX_VAR)*normal(IAXIS) + rtanvelx
                    if (NDIM .GE. 2) then
                        solnData(VELY_VAR, icell, jcell, kcell) =         &
                            midstate(VELX_VAR)*normal(JAXIS) + rtanvely 
                    end if
                    if (NDIM .GE. 3) then
                        solnData(VELZ_VAR, icell, jcell, kcell) =         &
                            midstate(VELX_VAR)*normal(KAXIS) + rtanvelz 
                    end if

                    ! Since we are putting our "real" state in lstate and our
                    ! "ghost" state in rstate, and the component has changed to
                    ! match the "ghost" state, the correct dens is the right
                    ! dens from the riemann solution which is in G0DENS since in
                    ! the above we put rstate in first, i.e as "lstate"
                    solnData(DENS_VAR, icell, jcell, kcell) = midstate(G0DENS_VAR) 
                    solnData(PRES_VAR, icell, jcell, kcell) = midstate(PRES_VAR) 
                    solnData(GAME_VAR, icell, jcell, kcell) = midstate(GAME_VAR)
#ifdef FLASH_3T
                    ! Since the component has changed, we need to use the
                    ! pfractions from the "ghost state" which was set above    
                    ! See the notes in hy_ppm_block for more details
                    solnData(PION_VAR,icell,jcell,kcell) = midstate(PRES_VAR)*pfracIon
                    solnData(PELE_VAR,icell,jcell,kcell) = midstate(PRES_VAR)*pfracEle
                    solnData(PRAD_VAR,icell,jcell,kcell) = midstate(PRES_VAR)*pfracRad
                    call FT_eosFindE(icell, jcell, kcell,                            &
                         midstate(PRES_VAR)*pfracIon,midstate(G0DENS_VAR),0,intEIon)
                    call FT_eosFindE(icell, jcell, kcell,                            &
                         midstate(PRES_VAR)*pfracEle,midstate(G0DENS_VAR),1,intEEle)
                    call FT_eosFindE(icell, jcell, kcell,                            &
                         midstate(PRES_VAR)*pfracRad,midstate(G0DENS_VAR),2,intERad)
                    solnData(EION_VAR, icell, jcell, kcell) = intEIon
                    solnData(EELE_VAR, icell, jcell, kcell) = intEEle
                    solnData(ERAD_VAR, icell, jcell, kcell) = intERad
#endif
                    ! Now I need to fill in the rest of the variables that don't
                    ! come out of the Riemann solve (species/rad energy grps).  
                    ! I want to use the new state, i.e the ghost state of neighbors
                    ! of the rest of the variables.
                    do ivar = SPECIES_BEGIN, UNK_VARS_END
                        !call GhostState_1var(blkID, ivar, icell, jcell, kcell, gV) 
                        ! Changed to a stored array filled during FillGhostStates
                        gV = xnGhostVal(ivar-SPECIES_BEGIN+1, icell, jcell, kcell, newComp)
                        tmpSoln(ivar, icell, jcell, kcell) = gV     
                    end do
#ifdef FLASH_3T
                    ! There is no reason that the rad energy groups we have just
                    ! filled in will be consistent with the rad energy we have
                    ! added, so we need to scale the rad energy groups to make
                    ! them consistent
                    ! FIXME: The 2 is a hack to remoe SUMY and YE Mass scalars,
                    ! should be fixed to reference the end of MGDR variables
                    sumEne = 0.0
                    do ivar = SPECIES_BEGIN+NSPECIES, MASS_SCALARS_END-2
                        sumEne = sumEne + tmpSoln(ivar, icell, jcell, kcell)
                    end do
                    do ivar = SPECIES_BEGIN+NSPECIES, MASS_SCALARS_END-2
                        eneRat = intERad / sumEne
                        tmpSoln(ivar, icell, jcell, kcell) =                    &
                                   tmpSoln(ivar, icell, jcell, kcell) * eneRat
                    end do
#endif
                    ! Store the index of the updated cell, so we know which to
                    ! copy over after all cells have been updated
                    updatedIndices(cnt, 1) = icell
                    updatedIndices(cnt, 2) = jcell
                    updatedIndices(cnt, 3) = kcell
                    ! Keep a count of how many cells have been updated
                    cnt = cnt + 1
                else if (newComp .eq. 0) then
                    ! The front thinks its a boundary, verify this and if it is
                    ! a boundary then leave component alone as the boundary
                    ! update which is called next will handle it accordingly
                    if (icell .lt. blkLimits(LOW,IAXIS) .or.    &
                        icell .gt. blkLimits(HIGH,IAXIS)) then
                        solnData(COMP_VAR, icell, jcell, kcell) = oldComp
                    else if (jcell .lt. blkLimits(LOW,JAXIS) .or.    &
                             jcell .gt. blkLimits(HIGH,JAXIS)) then
                        solnData(COMP_VAR, icell, jcell, kcell) = oldComp
                    else if (kcell .lt. blkLimits(LOW,KAXIS) .or.    &
                             kcell .gt. blkLimits(HIGH,KAXIS)) then
                        solnData(COMP_VAR, icell, jcell, kcell) = oldComp
                    else
                        call Driver_abortFlash('Front thinks a boundary is in the domain')
                    end if
                else
                    solnData(COMP_VAR, icell, jcell, kcell) = newComp
                end if
            end do  
        end do
    end do

    if (cnt .gt. 1) then
        ! Since we started cnt at 1, we need to go up to and including cnt-1,
        ! not cnt
        do i = 1, cnt - 1
            icell = updatedIndices(i, 1)
            jcell = updatedIndices(i, 2)
            kcell = updatedIndices(i, 3)
            solnData(COMP_VAR, icell, jcell, kcell) = tmpSoln(COMP_VAR, icell, jcell, kcell)
            do ivar = SPECIES_BEGIN, UNK_VARS_END
                solnData(ivar, icell, jcell, kcell) = tmpSoln(ivar, icell, jcell, kcell)
            end do
        end do
    end if
    
    !JAM - Added EOS call (needed to update the state var before the hydro step
    ! FIXME: Ideally this would go over only those cell which were updated but
    ! for now we will just call the EOS update on the whole grid 
#ifdef FLASH_3T
    call Eos_wrapped(MODE_DENS_EI_GATHER, blkLimits, blkID)
#else
    call Eos_wrapped(MODE_DENS_PRES, blkLimits, blkID) 
#endif

#ifndef FIXEDBLOCKSIZE
    deallocate(iCoord)
    deallocate(jCoord)
    deallocate(kCoord)
    deallocate(tmpSoln)
#endif

contains
    subroutine GhostState_1var (blkID, ivar, ii, jj, kk, ghostVal)

      use Driver_interface, ONLY : Driver_abortFlash
      use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits
      implicit none

      integer, INTENT(IN)                 :: ivar
      integer, INTENT(IN)                 :: ii
      integer, INTENT(IN)                 :: jj
      integer, INTENT(IN)                 :: kk
      integer, INTENT(IN)                 :: blkID
      real, INTENT(INOUT)                 :: ghostVal

      integer                             :: i, j, k, comp, counter
      integer                             :: iii, jjj, kkk, ghostVal1, counter1
      integer, dimension(2, MDIM)         :: blkLimits, blkLimitsGC
      real, pointer                       :: solnData(:,:,:,:)

      call Grid_getBlkPtr(blkID,solnData)
      call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC)

      ! We set the comp to the opposite comp of the cell as we are trying to
      ! fill in the ghost state for that cell. We will use this to compare to
      ! neighbors to get an average of the other component neighbors for this
      ! variable.  **NOTE: the cell's component hasn't been updated yet, so this
      ! is pulling the correct ghost state!
      if (solnData(COMP_VAR, ii, jj, kk) .eq. 1) then
          comp = 2
      else if (solnData(COMP_VAR, ii, jj, kk) .eq. 2) then
          comp = 1
      else
          call Driver_abortFlash('Non-physical component in FrontTracking_updateComp')
      end if
      ghostVal = 0.0
      counter = 0

      ! First sweep through neighbors to take average of any with the same
      ! component as what we want the ghost cell to be
      do i = ii - 1, ii + 1
          do j = jj - 1, jj + 1
              do k = kk - 1, kk + 1
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
                  if (solnData(COMP_VAR, i, j, k) .eq. comp) then 
                      ghostVal = ghostVal + solnData(ivar, i, j, k) 
                      counter = counter + 1
                  end if
              end do
          end do
      end do
   
      ! If we didn't find any opposite components in the first sweep we need to
      ! expand our search area, by having each neighbor loop through their
      ! neighbors to determine their avg state and then use that 
      if (counter .eq. 0) then
          do i = ii - 1, ii + 1
            do j = jj - 1, jj + 1
              do k = kk - 1, kk + 1
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
                  ghostVal1 = 0.0
                  counter1 = 0
                  do iii = i - 1, i + 1
                    do jjj = j - 1, j + 1 
                      do kkk = k - 1, k + 1 
                          if (iii .lt. 1 .or. jjj .lt. 1 .or. kkk .lt. 1) then
                              cycle
                          end if
                          if (iii .gt. blkLimitsGC(HIGH,IAXIS)) then
                              cycle
                          else if (jjj .gt. blkLimitsGC(HIGH,JAXIS)) then
                              cycle
                          else if (kkk .gt. blkLimitsGC(HIGH,KAXIS)) then
                              cycle
                          end if
                          if (solnData(COMP_VAR, i, j, k) .eq. comp) then       
                              ghostVal1 = ghostVal1 + solnData(ivar, i, j, k)
                              counter1 = counter1 + 1
                          end if
                      end do
                    end do
                  end do
                  if (counter1 .ne. 0) then
                      ghostVal1 = ghostVal1 / counter1
                      ghostVal = ghostVal + ghostVal1
                      counter = counter + 1
                  end if
              end do
            end do
          end do
      end if

      ghostVal = ghostVal / counter
      
      if (counter .eq. 0 ) then
          call Driver_abortFlash('Error, no nearby same component neighbors found')
      end if

      return

    end subroutine GhostState_1var

end subroutine
