!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/rieman
!!
!! NAME
!! 
!!  RiemSolver: Riemann Solver for a left and right state
!!
!! SYNOPSIS
!!
!!       call RiemSolver  ( STATE(IN)    :: lstate, 
!!                          STATE(IN)    :: rstate, 
!!                          STATE(OUT)   :: midstate)
!!                 
!!
!! DESCRIPTION
!!  
!!  Solve riemann shock tube problem for a general equation of state using 
!!  the method of Colella and Glaz.  Use a two shock approximation, and
!!  linearly interpolation between the head and tail of a rarefaction to
!!  treat rarefactions.
!!
!!  Take as input the effective left and right states, obtained by 
!!  integrating the parabolic reconstructions of the data over the domain
!!  of dependence of each of the characteristics on the respective side
!!  of the interface.  This is accomplished by states().  Return the
!!  solution to Riemann's problem on axis -- this is used in computing
!!  the fluxes.
!!
!!  The Riemann problem for the Euler's equation produces 4 states, 
!!  separated by the three characteristics (u - cs, u, u + cs):
!!
!!
!!        l_1      t    l_2       l_3
!!         \       ^   .       /
!!          \  *L  |   . *R   /
!!           \     |  .     /
!!            \    |  .    /
!!        L    \   | .   /    R
!!              \  | .  /
!!               \ |. / 
!!                \|./
!!       ----------+----------------> x
!!
!!       l_1 = u - cs   eigenvalue
!!       l_2 = u        eigenvalue (contact)
!!       l_3 = u + cs   eigenvalue
!!
!!       only density jumps across l_2
!!
!!  References:
!!
!!   CG:   Colella & Glaz 1985, JCP, 59, 264.
!!
!!   CW:   Colella & Woodward 1984, JCP, 54, 174.
!!
!!   Fry:  Fryxell et al. 2000, ApJS, 131, 273.
!!
!!   Toro: Toro 1999, ``Riemann Solvers and Numerical Methods for Fluid
!!         Dynamcs: A Practical Introduction, 2nd Ed.'', Springer-Verlag
!!
!! ARGUMENTS
!!    lstate: Left State
!!    rstate: Right State
!!    midstate: Output midstate
!!
!!***
subroutine RiemSolver (lstate, rstate, midstate)                   

#include "Flash.h"

  use Hydro_data, ONLY: hy_numXn, &
                        hy_utlft,  hy_utrght, &
                        hy_uttlft, hy_uttrgt, &
                        hy_eiLft,  hy_eiRght, &
                        hy_xnlft,  hy_xnrght, &
                        hy_smallp, hy_smallu, &
                        hy_smlrho, hy_nriem,  &
                        hy_riemanTol

  use Driver_interface, ONLY : Driver_abortFlash


  implicit none
!! Arguments ---------------------------------- 

  real, intent(IN),  DIMENSION(NPROP_VARS) :: lstate, rstate
  real, intent(OUT), DIMENSION(NPROP_VARS) :: midstate

!! Local variable ---------------------------------
  real :: wlft, wrght, pstar, ustar, vstar, cestar, &
       rhostr, westar, ps, us, uts, utts, vs, rhos, ces, ws, wes, &
       gmstar, games, gamcs, lRhoC, rRhoc, lInvRho, rInvRho

  ! JAM - Added to do the wave curve dens expansions for left and right
  real :: ldens_ms, rdens_ms, coef1, coef2, coef3, coef4
  real :: mf_left, mf_rght, p_rat
  
  real :: pstar1, pstar2, gmstrl, gmstrr, &
       &   wlft1, wrght1, gmin, gmax, &
       &   gamfac
  real :: laux, raux
  real :: scrch1, scrch2, scrch3, scrch4

  real ::  ge, gc, ustrl1, ustrr1, ustrl2, ustrr2, &
       & delu1, delu2, pres_err
 
  real :: rhoav, pav, uav, gameav, ugrdl, urell

  integer :: i, j, k, n, ierr
  
  character(len=1), save :: dirs(3) = (/ 'x', 'y', 'z' /)
  
  real, parameter :: small_dp = 1.e2 * epsilon(1.e0)
  real, DIMENSION(hy_nriem+2) :: riem_pstor
  
    !--------------------------------------------------------------------------
    do i = 1, NPROP_VARS
        midstate(i) = 0.0
    end do


    ! FIXME: Need this to come from EOS
    lRhoc = sqrt(lstate(GAMC_VAR)*lstate(PRES_VAR)*lstate(DENS_VAR))
    rRhoc = sqrt(rstate(GAMC_VAR)*rstate(PRES_VAR)*rstate(DENS_VAR))
    lInvRho = 1./lstate(DENS_VAR)
    rInvRho = 1./rstate(DENS_VAR)

    ! FIXME: What is this used for?
    ugrdl = 0

    ! We carry around two adiabatic indices (gammas).  gamc is usually
    ! referred to as gamma_1, and is d(log P)/d(log rho) (CG Eq. 7).  game
    ! is CG Eq. 8.      

    ! calculate limits on gamma based on the values in the neighboring zones
    ! gamfac is the max gamma factor in CG Eq. 31, and is used over and over,
    ! so store it for efficiency.

    ! TODO: Determine if the 5 cell buffer is needed for the calculation of two
    ! input states NO!

     laux      = sqrt(0.5 * (lstate(GAME_VAR) - 1.0) / lstate(GAME_VAR))
     raux      = sqrt(0.5 * (rstate(GAME_VAR) - 1.0) / rstate(GAME_VAR))
     ge        = 0.5 * (lstate(GAME_VAR) + rstate(GAME_VAR))
     gc        = 0.5 * (lstate(GAMC_VAR) + rstate(GAMC_VAR))
     gamfac    = (1. - ge / gc) * (ge - 1.)
     gmin      = min (lstate(GAME_VAR), rstate(GAME_VAR))
     gmax      = max (lstate(GAME_VAR), rstate(GAME_VAR))
  
    ! construct first guess for secant iteration by assuming that the nonlinear 
    ! wave speed is equal to the sound speed -- the resulting expression is the
    ! same as Toro, Eq. 9.28 in the Primitive Variable Riemann Solver (PVRS).
    ! See also Fry Eq. 72.
  
    ! XXX: hy_crght and hy_clft are rho*c! 
     pstar1 = rstate(PRES_VAR)-lstate(PRES_VAR)-rRhoc*(rstate(VELX_VAR)-lstate(VELX_VAR))
     pstar1 = lstate(PRES_VAR) + pstar1 * (lRhoc / (lrhoC + rRhoc))
     pstar1 = max(hy_smallp, pstar1)

    ! calculate approximation jump in gamma acrosss the interface based on the 
    ! first guess for the pressure jump.  There is a left and right 'star' region,
    ! so we need gamma add both places.  Use CG Eq. 31 and 32, with definitions
    ! as in CG Eq. 33.
    
     gmstrl = gamfac * (pstar1 - lstate(PRES_VAR))
     gmstrl = lstate(GAME_VAR) + 2.e0 * gmstrl / (pstar1 + lstate(PRES_VAR))
     
     gmstrr = gamfac * (pstar1 - rstate(PRES_VAR))
     gmstrr = rstate(GAME_VAR) + 2.e0 * gmstrr / (pstar1 + rstate(PRES_VAR))
     
     gmstrl = max(gmin, min(gmstrl, gmax))
     gmstrr = max(gmin, min(gmstrr, gmax))

    ! calculate nonlinear wave speeds for the left and right moving waves based
    ! on the first guess for the pressure jump.  Again, there is a left and a 
    ! right wave speed.  Compute this using CG Eq. 34.
    
     scrch1 = pstar1 - (gmstrl - 1.e0) * lstate(PRES_VAR) &
          & / (lstate(GAME_VAR) - 1.e0)
     if (scrch1 .EQ. 0.e0) scrch1 = hy_smallp
     
     wlft1  = pstar1 + 0.5e0 * (gmstrl - 1.e0) * (pstar1 + lstate(PRES_VAR))
     wlft1  = (pstar1 - lstate(PRES_VAR)) * wlft1 / (lInvRho * scrch1)
     wlft1  = sqrt(abs(wlft1))
     
     scrch2 = pstar1-(gmstrr - 1.e0)*rstate(PRES_VAR)/(rstate(GAME_VAR)-1.e0)
     
     if (scrch2 .EQ. 0.e0) scrch2 = hy_smallp
     
     wrght1 = pstar1 + 0.5e0 * (gmstrr - 1.e0) * (pstar1 + rstate(PRES_VAR))
     wrght1 = (pstar1 - rstate(PRES_VAR)) * wrght1 / (rInvRho * scrch2)
     wrght1 = sqrt(abs(wrght1))
     
    ! if the pressure jump is small, the wave speed is just the sound speed

     if (abs(pstar1-lstate(PRES_VAR)) < small_dp*(pstar1+lstate(PRES_VAR))) wlft1  = lRhoc
     wlft1  = max (wlft1,  laux * lRhoc)
    
     if (abs(pstar1-rstate(PRES_VAR)) < small_dp*(pstar1+rstate(PRES_VAR))) wrght1 = rRhoc
     wrght1 = max (wrght1, raux * rRhoc)

    ! construct second guess for the pressure using the nonlinear wave speeds
    ! from the first guess.  This is basically the same thing we did to get
    ! pstar1, except now we are using the better wave speeds instead of the 
    ! sound speed.

     pstar2 = rstate(PRES_VAR)-lstate(PRES_VAR)-wrght1*(rstate(VELX_VAR)-lstate(VELX_VAR))
     pstar2 = lstate(PRES_VAR) + pstar2 * wlft1 / (wlft1 + wrght1)
     pstar2 = max (hy_smallp, pstar2)

    ! begin the secant iteration -- see CG Eq. 17 for details.  We will continue to
    ! interate for convergence until the error falls below tol (in which case, 
    ! things are good), or we hit hy_nriem iterations (in which case we have a 
    ! problem, and we spit out an error).

     riem_pstor(1) = pstar1
     riem_pstor(2) = pstar2
     
     do n = 1, hy_nriem
        
        ! new values for the gamma at the "star" state -- again, using CG Eq. 31
          
        gmstrl = gamfac * (pstar2 - lstate(PRES_VAR))
        gmstrl = lstate(GAME_VAR) + 2.e0 * gmstrl / (pstar2 + lstate(PRES_VAR))
        
        gmstrr = gamfac * (pstar2 - rstate(PRES_VAR))
        gmstrr = rstate(GAME_VAR) + 2.e0 * gmstrr / (pstar2 + rstate(PRES_VAR))
        
        gmstrl = max (gmin, min (gmax, gmstrl))
        gmstrr = max (gmin, min (gmax, gmstrr))
        
        ! new nonlinear wave speeds, using CG Eq. 34 and the updated gammas
          
        scrch1 = pstar2-(gmstrl-1.e0)*lstate(PRES_VAR)/(lstate(GAME_VAR)-1.e0)
        if (scrch1 .EQ. 0.e0) scrch1 = hy_smallp
        
        wlft  = pstar2 + 0.5e0 * (gmstrl - 1.e0) * (pstar2 + lstate(PRES_VAR))
        wlft  = (pstar2 - lstate(PRES_VAR)) * wlft  / (lInvRho * scrch1)
        wlft  = sqrt(abs(wlft))

        scrch2 = pstar2-(gmstrr-1.e0) *rstate(PRES_VAR)/(rstate(GAME_VAR)-1.e0)
        if (scrch2 .EQ. 0.e0) scrch2 = hy_smallp
        
        wrght = pstar2 + 0.5e0 * (gmstrr - 1.e0) * (pstar2 + rstate(PRES_VAR))
        wrght = (pstar2 - rstate(PRES_VAR)) * wrght / (rInvRho * scrch2)
        wrght = sqrt(abs(wrght))
        
        ! if the pressure jump is small, the wave speed is just the sound speed

        if (abs(pstar2-lstate(PRES_VAR)) < small_dp*(pstar2+lstate(PRES_VAR))) wlft  = lRhoc
        wlft  = max (wlft,  laux * lRhoc)
        
        if (abs(pstar2-rstate(PRES_VAR)) < small_dp*(pstar2+rstate(PRES_VAR))) wrght = rRhoc
        wrght = max (wrght, raux * rRhoc)

        ! compute the velocities in the "star" state -- using CG Eq. 18 -- ustrl2 and
        ! ustrr2 are the velocities they define there.  ustrl1 and ustrl2 seem to be
        ! the velocities at the last time, since pstar1 is the old 'star' pressure, and
        ! wlft1 is the old wave speed.
        
        ustrl1 = lstate(VELX_VAR) - (pstar1 - lstate(PRES_VAR)) /  wlft1
        ustrr1 = rstate(VELX_VAR) + (pstar1 - rstate(PRES_VAR)) / wrght1
        ustrl2 = lstate(VELX_VAR) - (pstar2 - lstate(PRES_VAR)) /   wlft
        ustrr2 = rstate(VELX_VAR) + (pstar2 - rstate(PRES_VAR)) /  wrght
        
        delu1  = ustrl1 - ustrr1
        delu2  = ustrl2 - ustrr2
        scrch1 = delu2  - delu1
        
        if (abs(pstar2-pstar1) .le. hy_smallp) scrch1 = 0.e0
        
        if (abs(scrch1) .lt. hy_smallu) then
           delu2  = 0.e0
           scrch1 = 1.e0
        endif

        ! pressure at the "star" state -- using CG Eq. 18

        pstar  = pstar2 - delu2 * (pstar2 - pstar1) / scrch1
        pstar  = max (hy_smallp, pstar)
        
        ! check for convergence of iteration, hy_riemanTol is a run-time parameter
        
        pres_err = abs(pstar-pstar2) / pstar
        if (pres_err .lt. hy_riemanTol) goto 10
        
        ! reset variables for next iteration
          
        pstar1 = pstar2
        pstar2 = pstar
        riem_pstor(n+2) = pstar
        
        wlft1  =  wlft
        wrght1 = wrght
        
     enddo
     
     n = n - 1
     
     ! print error message and stop code if iteration fails to converge
     
     print *, ' '
     print *, 'Nonconvergence in subroutine RiemSolver'
     print *, ' '
     print *, 'Iterations tried = ', n+2
     print *, 'Pressure error   = ', pres_err
     print *, 'rieman_tol       = ', hy_riemanTol
     print *, ' '
     print *, 'pL       = ', lstate(PRES_VAR),  ' pR       =', rstate(PRES_VAR)
     print *, 'uL       = ', lstate(VELX_VAR),  ' uR       =', rstate(VELX_VAR)
     print *, 'cL       = ', lRhoc,             ' cR       =', rRhoc
     print *, 'gamma_eL = ', lstate(GAME_VAR),  ' gamma_eR =', rstate(GAME_VAR)
     print *, 'gamma_cL = ', lstate(GAMC_VAR),  ' gamma_cR =', rstate(GAMC_VAR)
     print *, ' '
     print *, 'Iteration history:'
     print *, ' '
     print '(A4, 2X, A20)', 'n', 'p*'
     do j = 1, n+2
        print '(I4, 2X, E20.12)', j, riem_pstor(j)
     enddo
     print *, ' '
     print *, 'Terminating execution.'
     call Driver_abortFlash('Nonconvergence in subroutine rieman')
       
     ! land here if the iterations have converged
       
10   continue

! end of secant iteration

! calculate fluid velocity for the "star" state -- this comes from the shock
! jump equations, Fry Eq. 68 and 69.  The ustar velocity can be computed
! using either the jump eq. for a left moving or right moving shock -- we use
! the average of the two.
! NOTE: Also look at Fry Eqn. 75 and 76.

     scrch3 = lstate(VELX_VAR) - (pstar - lstate(PRES_VAR)) /  wlft
     scrch4 = rstate(VELX_VAR) + (pstar - rstate(PRES_VAR)) / wrght
     ustar  = 0.5e0 * (scrch3 + scrch4)

! account for grid velocity

     ! XXX: I don't understand what this is trying to do
     ! FIXME: Understand what is in ugrdl, seems like input?
     ! TODO: Do we need this?, what is a grid velocity?
     urell  = ustar - ugrdl
     scrch1 = sign (0.5, urell) - sign(0.5,-urell)     

! decide which state is located at the zone iterface based on the values 
! of the wave speeds.  This is just saying that if ustar > 0, then the state
! is U_L.  if ustar < 0, then the state on the axis is U_R.

     scrch2 = 0.5e0 * ( 1.e0 + scrch1)
     scrch3 = 0.5e0 * ( 1.e0 - scrch1)
     
     ps    = lstate(PRES_VAR) * scrch2 + rstate(PRES_VAR) * scrch3
     us    = lstate(VELX_VAR) * scrch2 + rstate(VELX_VAR) * scrch3
     vs    = lInvRho          * scrch2 + rInvRho          * scrch3 !v for v=1/rho
     games = lstate(GAME_VAR) * scrch2 + rstate(GAME_VAR) * scrch3
     gamcs = lstate(GAMC_VAR) * scrch2 + rstate(GAMC_VAR) * scrch3

     ! FIXME: Do I need these?? Don't think they matter for what I am doing?
     !uts(i)   = hy_utlft  * scrch2 + hy_utrght(i) * scrch3(i)
     !utts(i)  = hy_uttlft * scrch2 + hy_uttrgt(i) * scrch3(i)

     rhos  = 1.e0 / vs
     rhos  = max (hy_smlrho, rhos)
     
     vs    = 1.e0 / rhos
     ws    = wlft * scrch2 + wrght * scrch3
     ces   = sqrt (gamcs * ps * vs)
     
     ! compute rhostar, using the shock jump condition (Fry Eq. 80)
     
     vstar  = vs - (pstar - ps) / ws / ws
     rhostr = 1.e0 / vstar
     cestar = sqrt (gamcs * pstar * vstar)
     
! compute some factors, Fry Eq. 81 and 82       

     wes    = ces    - scrch1 * us
     westar = cestar - scrch1 * ustar
     
     scrch4 = ws * vs - scrch1 * us
      
     if (pstar - ps .ge. 0.e0) then
        wes    = scrch4
        westar = scrch4
     endif
     
     ! FIXME: Understand where ugrdl comes from and how it is used...?
     wes    = wes    + scrch1 * ugrdl
     westar = westar + scrch1 * ugrdl


  ! compute Fry Eq. 86
     gamfac = (1.e0 - games / gamcs) * (games - 1.e0)
     gmstar = gamfac * (pstar - ps)
     gmstar = games + 2.e0 * gmstar / (pstar + ps)
     gmstar = max (gmin, min (gmax, gmstar))

     ! JAM: For 3T implementations, we need these weights to determine the
     ! correct midstate fractions of each component for pressures and energies
     ! of the midstate
     midstate(G0PRES_VAR) = scrch2
     midstate(G1PRES_VAR) = scrch3

     ! FIXME: Do I need to calculate this???
     !eintAv = hy_eiLft(i) * scrch2(i) + hy_eiRght(i) * scrch3(i)
  
  !do n = 1, hy_numXn 
  !    xnav(n) = hy_xnlft(n) * scrch2 + hy_xnrght(n) * scrch3
  !enddo
  
! compute correct state for rarefaction fan by linear interpolation

     scrch1 = max (wes - westar, wes + westar, hy_smallu)
     scrch1 =     (wes + westar) / scrch1
     
     scrch1 = 0.5e0 * (1.e0 + scrch1)
     scrch2 =          1.e0 - scrch1
     
     rhoav  = scrch1 * rhostr + scrch2 * rhos
     uav    = scrch1 * ustar  + scrch2 * us
     ! FIXME: Do I need these?
     !utav   = uts
     !uttav  = utts
     pav    = scrch1 * pstar  + scrch2 * ps
     gameav = scrch1 * gmstar + scrch2 * games
  
     if (westar .ge. 0.e0) then
        rhoav  = rhostr
        uav    = ustar
        pav    = pstar
        gameav = gmstar
     endif
     
     if (wes .lt. 0.e0) then
        rhoav  = rhos
        uav    = us
        pav    = ps
        gameav = games
     endif
     
     ! TODO: Understand where ugrdl comes from and how it is used...?
     urell = uav - ugrdl

    !JAM - Calculate left state and right state using mass flux and wave types
    ! See cFgriemann.cpp:279 and tags within for details
    ! First we will calculate the left and right mass flux (cFeos.cpp:547)
    if (pstar < lstate(PRES_VAR)) then
        if (pstar .le. 0.0) then
           print*, "ERROR, pstar is negative"
           call Driver_abortFlash('ERROR, pstar is negative')
        end if
        coef3 = 0.5 * (lstate(GAMC_VAR) - 1.0) / lstate(GAMC_VAR)
        p_rat = pstar/lstate(PRES_VAR)
        mf_left = coef3 * (1.0 - p_rat) / (1.0 - p_rat**coef3)    
        mf_left = lRhoc*mf_left
    else
        coef1 = 0.5 * (lstate(GAMC_VAR) + 1.0)
        coef2 = 0.5 * (lstate(GAMC_VAR) - 1.0)
        mf_left = SQRT(lstate(DENS_VAR)*(coef1*pstar + coef2*lstate(PRES_VAR)))
    end if
    
    if (pstar < rstate(PRES_VAR)) then
        if (pstar .le. 0.0) then
           print*, "ERROR, pstar is negative"
           call Driver_abortFlash('ERROR, pstar is negative')
        end if
        coef3 = 0.5 * (rstate(GAMC_VAR) - 1.0) / rstate(GAMC_VAR)
        p_rat = pstar/rstate(PRES_VAR)
        mf_rght = coef3 * (1.0 - p_rat) / (1.0 - p_rat**coef3)
        mf_rght = rRhoc*mf_rght
    else
        coef1 = 0.5 * (rstate(GAMC_VAR) + 1.0)
        coef2 = 0.5 * (rstate(GAMC_VAR) - 1.0)
        mf_rght = SQRT(rstate(DENS_VAR)*(coef1*pstar + coef2*rstate(PRES_VAR)))
    end if

    ! Next we determine if the wave is a shock or rarefaction in that direction
    ! and then calculate the midstate appropriately, using dens_Hugoniot or 
    ! state on adiabat with pr functions from cFeos.cpp
    if (pstar < lstate(PRES_VAR)) then  ! RAREFACTION
        p_rat = pstar/lstate(PRES_VAR)
        ldens_ms = lstate(DENS_VAR)*p_rat**(1.0/lstate(GAMC_VAR))
    else                                ! SHOCK
        ! We use "1.0" as coeff on mf for lstate and "-1.0" for rstate
        ldens_ms = mf_left/(ustar - lstate(VELX_VAR) + (mf_left/lstate(DENS_VAR)))
        if ((pstar - lstate(PRES_VAR))*(ldens_ms - lstate(DENS_VAR)) < 0.0) then
            coef4 = (lstate(GAMC_VAR) - 1.0) / (lstate(GAMC_VAR) + 1.0)
            ldens_ms = lstate(DENS_VAR)*(pstar + lstate(PRES_VAR)*coef4) /         &
                                        (lstate(PRES_VAR) + pstar*coef4)
        end if
    end if
    
    if (pstar < rstate(PRES_VAR)) then  ! RAREFACTION
        p_rat = pstar/rstate(PRES_VAR)
        rdens_ms = rstate(DENS_VAR)*p_rat**(1.0/rstate(GAMC_VAR))
    else                                ! SHOCK
        ! We use "1.0" as coeff on mf for lstate and "-1.0" for rstate
        rdens_ms = -mf_rght/(ustar - rstate(VELX_VAR) - (mf_rght/rstate(DENS_VAR)))
        if ((pstar - rstate(PRES_VAR))*(rdens_ms - rstate(DENS_VAR)) < 0.0) then
            coef4 = (rstate(GAMC_VAR) - 1.0) / (rstate(GAMC_VAR) + 1.0)
            rdens_ms = rstate(DENS_VAR)*(pstar + rstate(PRES_VAR)*coef4) /         &
                                        (rstate(PRES_VAR) + pstar*coef4)
        end if
    end if

    !JAM - Put solutions into midstate vector (these are the star states)
    midstate(DENS_VAR) = rhostr
    midstate(PRES_VAR) = pstar
    midstate(VELX_VAR) = ustar
    midstate(GAME_VAR) = gmstar
    ! JAM - Store the left and right dens in G0 (left) and G1 (right) dens
    ! arrays which have been addd for other purposes, so in midstate they go
    ! unused... this will be easier than creating new variables to pass in
    midstate(G0DENS_VAR) = ldens_ms
    midstate(G1DENS_VAR) = rdens_ms
  return
end subroutine RiemSolver 
  
