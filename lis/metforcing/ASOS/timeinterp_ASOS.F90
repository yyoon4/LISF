!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_ASOS
! \label{timeinterp_ASOS}
!
! !REVISION HISTORY: 
!  20 Jan 2006:   Yudong Tian: Initial Implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  19 Apr 2013: Jonathan Case; Included comparison against undefined value,
!                              set small -ve values to 0
!  21 Jun 2013: Soni Yatheendradas; Based on new TRMM3B42V6 code, which is
!               based on changes from earlier code to avoid (a) alternate file skip,
!               (b) jump to previous day TRMM, and (c) absence of rain rate weighting
!  06 Sep 2017: Minwook Kim; based on new TRMM 3B42V7 code that has changes
!
! !INTERFACE: 
 subroutine timeinterp_ASOS(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod, only           : LIS_rc, LIS_domain
    use LIS_FORC_AttributesMod
    use LIS_metforcingMod,  only    : LIS_FORC_Base_State, LIS_forc
    use LIS_logMod,         only    : LIS_verify, LIS_logunit 
    use ASOS_forcingMod, only : ASOS_struc
    use LIS_forecastMod, only : LIS_get_iteration_index

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP

    integer :: t, index1

    integer          :: status
    type(ESMF_Field) :: pcpField, cpcpField
    real, pointer    :: pcp(:), cpcp(:)
    real :: wt1,wt2 ! SY
    real*8 :: BreakpointTime ! SY
    real,  allocatable :: ratio(:)
    integer          :: m, mfactor, kk, k

!   ! SY: Begin calculating weights, note that weights are NOT linear interpolation-based
!   if ((ASOS_struc(n)%ASOSyr_TStepEnd .EQ. &
!        ASOS_struc(n)%ASOSyr_TStepStart) .AND. &
!       (ASOS_struc(n)%ASOSmo_TStepEnd .EQ. &
!        ASOS_struc(n)%ASOSmo_TStepStart) .AND. &
!       (ASOS_struc(n)%ASOSda_TStepEnd .EQ. &
!        ASOS_struc(n)%ASOSda_TStepStart) .AND. &
!       (ASOS_struc(n)%ASOShr_TStepEnd .EQ. &
!        ASOS_struc(n)%ASOShr_TStepStart)) then
!     wt1 = 0.5
!     wt2 = 0.5
!   else
!     BreakpointTime = (ASOS_struc(n)%ASOStime_TStepStart+&
!                       ASOS_struc(n)%ASOStime_TStepEnd)/2.0
!     wt1 = (abs(BreakpointTime-ASOS_struc(n)%LIS_timeAtTStepStart)) / &
!           (ASOS_struc(n)%LIS_timeAtTStepEnd- &
!            ASOS_struc(n)%LIS_timeAtTStepStart) ! SY
!     wt2 = 1.0 - wt1
!   endif
    wt1=(ASOS_struc(n)%ASOStime_TStepEnd-LIS_rc%time)/ &
      (ASOS_struc(n)%ASOStime_TStepEnd-ASOS_struc(n)%ASOStime_TStepStart)
    wt2=1.0-wt1
    ! SY: End calculating weights, note that weights are NOT linear interpolation-based
     
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status, 'Error: enable Rainf in forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)

!- Convert precipitation sum to rate:
    mfactor = LIS_rc%nensem(n)/ASOS_struc(n)%nIter
    do k=1,LIS_rc%ntiles(n)/mfactor
      do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)
        ! SY: Start modification
        if ( ( ASOS_struc(n)%metdata2(kk,1,index1) .ne. LIS_rc%udef .and. &
               ASOS_struc(n)%metdata2(kk,1,index1) >= 0.0 ) .OR. &
             ( ASOS_struc(n)%metdata1(kk,1,index1) .ne. LIS_rc%udef .and. &
               ASOS_struc(n)%metdata1(kk,1,index1) >= 0.0 ) ) then
          if ( ASOS_struc(n)%metdata2(kk,1,index1) .eq. LIS_rc%udef .or. &
                ASOS_struc(n)%metdata2(kk,1,index1) < 0.0 ) then
            pcp(t) = ASOS_struc(n)%metdata1(kk,1,index1) / 3600.0      ! mm/s
          elseif ( ASOS_struc(n)%metdata1(kk,1,index1) .eq. LIS_rc%udef .or. &
                ASOS_struc(n)%metdata1(kk,1,index1) < 0.0 ) then
            pcp(t) = ASOS_struc(n)%metdata2(kk,1,index1) / 3600.0      ! mm/s
          else
            pcp(t) = (wt1 * ASOS_struc(n)%metdata1(kk,1,index1) + wt2 * &
                       ASOS_struc(n)%metdata2(kk,1,index1)) / 3600.0    ! mm/s
          endif
         ! SY: End modification
        else
           pcp(t) = LIS_rc%udef
        endif
      ! J.Case -- ensure that small negative values are reset to 0.
        if ( pcp(t) /= LIS_rc%udef ) then
           if (pcp(t)  < 0)  pcp(t) = 0.0
        endif
      enddo
    enddo
!-- Applying a convective rainfall amount to observed precip field: 
!  elseif( LIS_FORC_CRainf%selectOpt == 1 ) then

!    allocate(ratio(LIS_rc%ntiles(n)))

!    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
!         rc=status)
!    call LIS_verify(status, 'Error: enable CRainf in forcing variables list')
!    
!    call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
!    call LIS_verify(status)

!!------------------------------------------------------------------------
!! Compute ratio between convective model precip and total model precip
!! so that it can be applied to the observed global precip
!!------------------------------------------------------------------------
!    do t = 1,LIS_rc%ntiles(n) 
!       if (pcp(t) .ne. 0.0 .and.  &
!            pcp(t) .ne. LIS_rc%udef .and.  &
!            cpcp(t) .ne. LIS_rc%udef) then
!          ratio(t) = cpcp(t) / pcp(t)
!          if (ratio(t) .gt. 1.0) ratio(t) = 1.0
!          if (ratio(t) .lt. 0.0) ratio(t) = 0.0
!       else
!          ratio(t) = 0.0
!       endif
!    enddo
!
!    do t = 1,LIS_rc%ntiles(n) 
!       index1 = LIS_domain(n)%tile(t)%index
!! J.Case (4/22/2013) -- make consistent with Stg4/NMQ forcing examples.
!! SY: Further modified to allow for averaging over 2 ASOS datasets
!!       if (LIS_forc(n,findex)%metdata2(1,index1) .ge. 0.0) then ! J.Case
!!       if ( LIS_forc(n,findex)%metdata2(1,index1) .ne. LIS_rc%udef .and. & ! J.Case ! SY
!!            LIS_forc(n,findex)%metdata2(1,index1) >= 0.0 ) then ! J.Case ! SY
!! SY: Start modification
!       if ( ( ASOS_struc(n)%metdata2(1,index1) .ne. LIS_rc%udef .and. &
!              ASOS_struc(n)%metdata2(1,index1) >= 0.0 ) .OR. &
!            ( ASOS_struc(n)%metdata1(1,index1) .ne. LIS_rc%udef .and. &
!              ASOS_struc(n)%metdata1(1,index1) >= 0.0 ) ) then
!          if ( ASOS_struc(n)%metdata2(1,index1) .eq. LIS_rc%udef .or. &
!               ASOS_struc(n)%metdata2(1,index1) < 0.0 ) then
!            pcp(t) = ASOS_struc(n)%metdata1(1,index1) / 3600.0       ! mm/s
!          elseif ( ASOS_struc(n)%metdata1(1,index1) .eq. LIS_rc%udef .or. &
!               ASOS_struc(n)%metdata1(1,index1) < 0.0 ) then
!            pcp(t) = ASOS_struc(n)%metdata2(1,index1) / 3600.0   ! mm/s
!          else
!            pcp(t) = (wt1 * ASOS_struc(n)%metdata1(1,index1) + wt2 * ASOS_struc(n)%metdata2(1,index1)) / 3600.0       ! mm/s
!          endif
!! SY: End modification
!          !pcp(t) = LIS_forc(n,findex)%metdata2(1,index1) / 3600.0       ! mm/s ! SY
!          cpcp(t) = ratio(t) * pcp(t)
!       else
!          pcp(t) = LIS_rc%udef
!          cpcp(t) = LIS_rc%udef
!       endif
!! J.Case -- ensure that small negative values are reset to 0.
!       if ( pcp(t) /= LIS_rc%udef ) then
!          if (pcp(t)  < 0)  pcp(t) = 0.0
!       endif
!       if ( cpcp(t) /= LIS_rc%udef ) then
!          if (cpcp(t) < 0) cpcp(t) = 0.0
!       endif
!    enddo
!    deallocate(ratio)

!  endif

  end subroutine timeinterp_ASOS
