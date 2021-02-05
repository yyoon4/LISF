!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_ASOS
! \label{finalize_ASOS}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 06Sep2017: Minwook Kim; based on new TRMM 3B42V7 code that has changes
! 
! !INTERFACE:
subroutine finalize_ASOS(findex)

! !USES:
  use ASOS_forcingMod, only : ASOS_struc
! !DESCRIPTION: 
!  Routine to cleanup ASOS forcing related memory allocations.   
!
!EOP
  implicit none
  integer :: findex
  
  deallocate(ASOS_struc)

end subroutine finalize_ASOS
