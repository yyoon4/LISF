!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_ASOS
! \label{readcrd_ASOS}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 25 Aug 2006; Yudong Tian, Modification for 3B42, LIS 4.2 release 
! 21 Jun 2013: Soni Yatheendradas; change from earlier 3B42V6 code,
!              ported to 3B42V7
! 06 Sep 2017: Minwook Kim; Modification for ASOS
!
! !INTERFACE:
subroutine readcrd_ASOS()
! !USES:
  use ESMF
  use ASOS_forcingMod, only : ASOS_struc
  use LIS_coreMod, only : LIS_config,LIS_rc
  use LIS_logMod, only : LIS_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to ASOS forcing from
!  the LIS configuration file.
!
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config, "ASOS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ASOS_struc(n)%ASOSdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using ASOS forcing'
     write(LIS_logunit,*) 'ASOS forcing directory :',ASOS_struc(n)%ASOSDIR
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     !ASOS_struc(n)%ASOStime = 0.0 ! SY
  enddo

end subroutine readcrd_ASOS


