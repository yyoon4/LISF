!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noah271_dasnow_Mod
!BOP
!
! !MODULE: noah271_dasnow_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! !USES:        

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah271_dasnow_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP

  SAVE
contains
!BOP
! 
! !ROUTINE: noah271_dasnow_init
! \label{noah271_dasnow_init}
! 
! !INTERFACE:
  subroutine noah271_dasnow_init()
! !USES:
! !DESCRIPTION:        
!
!EOP
    implicit none
  end subroutine noah271_dasnow_init
end module noah271_dasnow_Mod
