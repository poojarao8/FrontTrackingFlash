subroutine FrontTracking_restartWrite(filename)

    use FrontTracking_data, ONLY: fr_useFrontTracking
    implicit none
#include "constants.h"

    character*MAX_STRING_LENGTH,intent(in) :: filename

    if(.NOT. fr_useFrontTracking) then
        return
    end if

    call FTAPI_WriteRestart(filename)
    write (6,*) "Wrote front restart output file."

end subroutine FrontTracking_restartWrite

