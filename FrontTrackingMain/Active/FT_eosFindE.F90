subroutine FT_eosFindE(icell, jcell, kcell, pres, dens, comp, intEne)

  use Multispecies_interface, ONLY : Multispecies_getSumInv
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_data, ONLY : eos_gasConstant
  use Grid_interface, ONLY : Grid_getBlkPtr
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  real, INTENT(IN)            :: pres
  real, INTENT(IN)            :: dens
  integer, INTENT(IN)         :: comp
  integer, INTENT(IN)         :: icell
  integer, INTENT(IN)         :: jcell
  integer, INTENT(IN)         :: kcell
  real, INTENT(OUT)           :: intEne

  real, pointer               ::  solnData(:,:,:,:) 
  real, DIMENSION(NSPECIES)   :: massfrac
  real                        :: abarInv, ionTemp, gam
  integer                     :: i
  integer                     :: blkID

  ! FIXME: Should be passed in, not hardcoded... or something...
  blkID = 1
 
  select case (comp)
  case(0)                              ! Ion EOS (gamma is f(T))
      ! Loop over the cells and fill the species data (this
      ! assumes that
      ! the species data has already been filled in)
      call Grid_getBlkPtr(blkID,solnData)
      do i = 1, NSPECIES
          massfrac(i) = solnData(SPECIES_BEGIN + i - 1, icell, jcell, kcell)
      end do
      call Multispecies_getSumInv(A, abarinv, massfrac)
      ionTemp = (pres*abarinv/(dens*eos_gasConstant))/11604505. ! Convert to keV
      gam = (7./4. + 5./3.) / 2. - (7./4. - 5./3.) / 2. * erf(15. * (ionTemp - 0.5))
  !    gam = 5./3.
      intEne = pres / ((gam - 1.) * dens)
  case(1)                                    ! Ele EOS (gamma is 5/3)
      intEne = pres / ((5./3. - 1.) * dens)
  case(2)                                    ! Rad EOS (gamma is 4/3)
      intEne = pres * 3.0 / dens
  case default 
      call Driver_abortFlash('Energy component not known FT_eosFindE')
  end select
 
end subroutine FT_eosFindE

