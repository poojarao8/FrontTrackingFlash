subroutine FrontTracking_getVel(tmp, coords, velocity, lcomp, rcomp)

#include "Flash.h"
#include "constants.h"
#include "ftapi.h"

    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                        Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
    use Driver_interface, ONLY : Driver_abortFlash
    use FrontTracking_interface, ONLY : FrontTracking_output
    
    implicit none
    real, POINTER, INTENT(IN)               ::  tmp
    real, INTENT(IN), dimension(MDIM)       ::  coords
    real, INTENT(OUT), dimension(MDIM)      ::  velocity
    integer, INTENT(IN)                     ::  lcomp
    integer, INTENT(IN)                     ::  rcomp

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
    real, dimension(MDIM)                   ::  normal
    real, dimension(NPROP_VARS)             ::  lstate, rstate, midstate
    integer                                 ::  imin, jmin, kmin
    real                                    :: lnormvel, rnormvel, ltanvelx, rtanvelx
    real                                    :: ltanvely, ltanvelz, rtanvely, rtanvelz

    ! JAM -- Added for FTAPI   FIXME: I would rather have this passed in than
    !                                 hardcoded, it is not compatible with 
    !                                 a lot of the options in FLASH
    blkID = 1

    velocity(:) = 0.0
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

    ! Determine the bounding box (imin, imin+1), (jmin, jmin+1), (kmin, kmin+1) of cells
    ! for the front pt, so that we know which cell states to use for the interpolation
    if (coords(IAXIS) .ge. iCoord(icell)) then
        imin = icell
    else
        imin = icell - 1 
    end if
    if (coords(JAXIS) .ge. jCoord(jcell)) then
        jmin = jcell
    else
        jmin = jcell - 1 
    end if
    if (coords(KAXIS) .ge. kCoord(kcell)) then
        kmin = kcell
    else
        kmin = kcell - 1 
    end if

    call interpolateStates(imin, jmin, kmin, coords, lstate, lcomp)
    call interpolateStates(imin, jmin, kmin, coords, rstate, rcomp)

    call FTAPI_normal(coords, normal)

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
    
    !FIXME: Probably need to make an EOS call here before I pass into RiemSolver??
    ! JAM: So far seems like it works without an EOS call
    call RiemSolver(lstate, rstate, midstate)

    ! JAM: We only propagate the front point in the normal direction, no tangential component
    ! Theory behind this is discussed in I-L Chern et al (1986) JCP 62 1 83-110
    velocity(IAXIS) = midstate(VELX_VAR)*normal(IAXIS) 
    if (NDIM .GE. 2) then
        velocity(JAXIS) = midstate(VELX_VAR)*normal(JAXIS)
    end if
    if (NDIM .GE. 3) then
        velocity(KAXIS) = midstate(VELX_VAR)*normal(KAXIS)
    end if

    ! Error Catching
    if (velocity(IAXIS) .ne. velocity(IAXIS)) then
        call Driver_abortFlash('NAN xvels in FrontTracking_getVel.F90')
    end if

    if (velocity(JAXIS) .ne. velocity(JAXIS)) then
        call Driver_abortFlash('NAN yvels in FrontTracking_getVel.F90')
    end if

    if (velocity(KAXIS) .ne. velocity(KAXIS)) then
        call Driver_abortFlash('NAN zvels in FrontTracking_getVel.F90')
    end if

    call Grid_releaseBlkPtr(blkID,solnData)

#ifndef FIXEDBLOCKSIZE
    deallocate(iCoord)
    deallocate(jCoord)
    deallocate(kCoord)
#endif

contains
    subroutine interpolateStates(imin, jmin, kmin, coords, state, comp)

        use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                            Grid_putBlkData, Grid_releaseBlkPtr, Grid_getCellCoords
        use Driver_interface, ONLY : Driver_abortFlash

        implicit none
        integer, INTENT(IN)                     ::  imin, jmin, kmin
        real,INTENT(IN),dimension(MDIM)         ::  coords
        real,INTENT(OUT),dimension(NPROP_VARS)  ::  state
        integer,INTENT(IN)                         ::  comp

        real, pointer                           ::  solnData(:,:,:,:)
        
#ifdef FIXEDBLOCKSIZE
        real,DIMENSION(MAXCELLS),target         ::  iCoord, jCoord, kCoord
#else
        real, allocatable,DIMENSION(:),target   ::  iCoord, jCoord, kCoord
#endif

        integer                                 ::  blkID, isizeGC, jsizeGC, ksizeGC
        integer                                 ::  icell, jcell, kcell
        integer                                 ::  i, j, k, ii, jj, kk
        integer, dimension(2, MDIM)             ::  blkLimits, blkLimitsGC
        logical                                 ::  gcell = .true.
        real, dimension(2)                      ::  fraci, fracj, frack
        real                                    ::  weight, dens, pres, intDens, intPres
        real                                    ::  xvel, yvel, zvel, game, gamc
        real                                    ::  intXvel, intYvel, intZvel
        real                                    ::  intGame, intGamc

        ! FIXME: This should not be hardcoded!
        blkID = 1

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
        fraci(1) = (coords(IAXIS) - iCoord(imin)) / (iCoord(imin+1) - iCoord(imin))
        if (NDIM .GE. 2) then 
            call Grid_getCellCoords(JAXIS,blkID,CENTER,gcell,jCoord,jsizeGC)
            fracj(1) = (coords(JAXIS) - jCoord(jmin)) / (jCoord(jmin+1) - jCoord(jmin))
        else
            jCoord(:) = 0.0
            fracj(1) = 1.0
        end if        
        if (NDIM .GE. 3) then 
            call Grid_getCellCoords(KAXIS,blkID,CENTER,gcell,kCoord,ksizeGC)
            frack(1) = (coords(KAXIS) - kCoord(kmin)) / (kCoord(kmin+1) - kCoord(kmin))
        else
            kCoord(:) = 0.0
            frack(1) = 1.0
        end if

        call Grid_getBlkPtr(blkID,solnData)

        fraci(2) = 1.0 - fraci(1)
        fracj(2) = 1.0 - fracj(1)
        frack(2) = 1.0 - frack(1)

        do i = 1, NPROP_VARS
            state(i) = 0.0
        end do

        intDens = 0.0
        intPres = 0.0
        intGame = 0.0
        intGamc = 0.0
        intXvel = 0.0
        intYvel = 0.0
        intZvel = 0.0
        ! JAM - Loop over bounding box (imin - imax, etc...)
        do i = imin, imin + 1
          do j = jmin, jmin + 1
            do k = kmin, kmin + 1
              if((i .lt. blkLimitsGC(LOW,IAXIS)) .or. (i .gt. blkLimitsGC(HIGH,IAXIS))) then
                cycle
              else if((j .lt. blkLimitsGC(LOW,JAXIS)) .or. (j .gt. blkLimitsGC(HIGH,JAXIS))) then
                cycle
              else if((k .lt. blkLimitsGC(LOW,KAXIS)) .or. (k .gt. blkLimitsGC(HIGH,KAXIS))) then
                cycle
              end if
              ! left bound will multiply by frac_(2) and right bound by frac_(1)
              ii = (imin - i) + 2
              if (NDIM .ge. 2) then
                  jj = (jmin - j) + 2
              else 
                  jj = 1.0
              end if
              if (NDIM .ge. 3) then
                  kk = (kmin - k) + 2
              else
                  kk = 1.0
              end if
              weight = fraci(ii)*fracj(jj)*frack(kk)
              if (comp .eq. 1) then
                dens = solnData(G0DENS_VAR, i, j, k) 
                pres = solnData(G0PRES_VAR, i, j, k)
                game = solnData(G0EGAM_VAR, i, j, k)
                gamc = solnData(G0CGAM_VAR, i, j, k)
                xvel = solnData(G0XVEL_VAR, i, j, k)
                yvel = solnData(G0YVEL_VAR, i, j, k)
                zvel = solnData(G0ZVEL_VAR, i, j, k)
              else
                dens = solnData(G1DENS_VAR, i, j, k) 
                pres = solnData(G1PRES_VAR, i, j, k)
                game = solnData(G1EGAM_VAR, i, j, k)
                gamc = solnData(G1CGAM_VAR, i, j, k)
                xvel = solnData(G1XVEL_VAR, i, j, k)
                yvel = solnData(G1YVEL_VAR, i, j, k)
                zvel = solnData(G1ZVEL_VAR, i, j, k)
              end if
              intDens = intDens + dens*weight
              intPres = intPres + pres*weight
              intGame = intGame + game*weight
              intGamc = intGamc + gamc*weight
              intXvel = intXvel + xvel*weight
              intYvel = intYvel + yvel*weight
              intZvel = intZvel + zvel*weight
            end do
          end do
        end do

        state(DENS_VAR) = intDens
        state(PRES_VAR) = intPres
        state(GAME_VAR) = intGame
        state(GAMC_VAR) = intGamc
        state(VELX_VAR) = intXvel
        state(VELY_VAR) = intYvel
        state(VELZ_VAR) = intZvel        

    end subroutine interpolateStates


end subroutine FrontTracking_getVel
