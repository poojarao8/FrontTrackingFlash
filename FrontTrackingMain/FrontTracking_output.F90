subroutine FrontTracking_output(step,chk)

#include "FTAPI.h"
#include "constants.h"

    use FrontTracking_data, ONLY: fr_useFrontTracking
    use IO_data, ONLY : io_checkpointFileNumber, io_plotFileNumber

    implicit none
    
    integer, INTENT(IN)         :: step
    character*MAX_STRING_LENGTH :: filename
    integer                     :: tmp
    integer, INTENT(IN)         :: chk

    if (.not. fr_useFrontTracking) then
        return
    end if

    if (chk .eq. 0) then
        write(filename, "(A6,I0.4)") "ft_plt", io_plotFileNumber
    else 
        write(filename, "(A6,I0.4)") "ft_chk", io_checkpointFileNumber
    end if

    tmp = LEN(TRIM(filename))

    if (chk .eq. 0) then
        call ftapi_output(filename, tmp)
    else
        call FTAPI_writeRestart(filename, tmp)
    end if

end subroutine FrontTracking_output
