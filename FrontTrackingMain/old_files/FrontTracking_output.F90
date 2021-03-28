subroutine FrontTracking_output(filename)
#include "FTAPI.h"
#include "constants.h"
    use FrontTracking_data, ONLY: fr_useFrontTracking
    implicit none

    character*MAX_STRING_LENGTH, intent(in) :: filename

    call ftapi_output(filename)

end subroutine FrontTracking_output
