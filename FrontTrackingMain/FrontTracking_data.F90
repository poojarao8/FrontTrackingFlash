!!****if* source/FrontTracking/FrontTrackingMain/FrontTracking_data
!!
!! NAME
!!    FrontTracking_data
!!
!! SYNOPSIS
!!    use FrontTracking_data, ONLY: fr_useFrontTracking, fr_restart
!!
!! DESCRIPTION
!!    Module to hold local variables and data types for FrontTracking unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!    useFrontTracking   BOOLEAN [TRUE]  Should front tracking be used in this simulation?
!!    fr_restart         BOOLEAN [TRUE]  Are we going to need restart the simulation?
!!***

module FrontTracking_data
  
  implicit none

  logical, save :: fr_useFrontTracking
  logical, save :: fr_restart
  
  !componenete ggid data struct
  integer, dimension(:), allocatable :: FrontTracking_compGrid 

end module FrontTracking_data
