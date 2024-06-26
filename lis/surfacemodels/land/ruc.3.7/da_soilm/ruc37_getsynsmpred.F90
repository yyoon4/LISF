!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: RUC37_getsynsmpred
! \label{RUC37_getsynsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine RUC37_getsynsmpred(n, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use RUC37_lsmMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the Soil moisture obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: count1(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  integer                :: i,t,m,gid

  obs_pred = 0.0
  count1 = 0.0
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        obs_pred(gid,m)= obs_pred(gid,m)+ RUC37_struc(n)%ruc37(t)%smc(1)
        count1(gid,m) = count1(gid,m) +1
     enddo
  enddo

  do i=1,LIS_rc%ngrid(n)
     do m=1,LIS_rc%nensem(n)
        obs_pred(i,m) = obs_pred(i,m)/(count1(i,m))
     enddo
  enddo

end subroutine RUC37_getsynsmpred

