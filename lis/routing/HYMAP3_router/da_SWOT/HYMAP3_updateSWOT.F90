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
! !ROUTINE: HYMAP3_updateSWOT
!  \label{HYMAP3_updateSWOT}
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
! 24 Mar 26: Yeosang Yoon; Updated the code to fit HyMAP3
! 10 Sep 26: Yeosang Yoon; Update on memory issue for localization
! 12 Sep 26: Added selectable MAX-AREA and BLEND-AREA methods
!
! !INTERFACE:
subroutine HYMAP3_updateSWOT(n, Routing_State, Routing_Incr_State)
! !USES:
  use ESMF
  use HYMAP3_daSWOT_Mod, only: HYMAP3_daSWOT_struc
  use HYMAP3_routingMod, only: HYMAP3_routing_struc
  use LIS_coreMod, only: LIS_rc, LIS_routing, LIS_localpet
  use LIS_logMod, only: LIS_verify, LIS_logunit, LIS_endrun

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
  type(ESMF_State)       :: Routing_Incr_State
!
! !DESCRIPTION:
!
!  This routine updates the water level prognostic variables.
!
!  Localization blending options:
!    1 = MAX-AREA   : retain the increment associated with the largest
!                     drainage-area-ratio weight at each target cell.
!    2 = BLEND-AREA : use all valid source increments at each indirect
!                     target cell according to
!
!       increment = max(weight) * sum(weight*increment) / sum(weight)
!
!  Direct DA increments are retained without blending for option 2.
!
!EOP

  type(ESMF_Field)       :: sfcelevField
  type(ESMF_Field)       :: sfcelevIncrField
  real, pointer          :: sfcelev(:)
  real, pointer          :: sfcelevIncr(:)
  real, allocatable      :: sfcelevIncr_tmp(:)
  real, allocatable      :: weightedIncrSum(:)
  real, allocatable      :: weightSum(:)
  real, allocatable      :: maxWeight(:)

  integer                :: t,i,i1,m
  integer                :: ix,iy
  integer                :: c,r,c1,c2,r1,r2,t1
  integer                :: siteid
  integer                :: status
  integer                :: nstate

  real                   :: localWeight(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                   :: pairWeight
  real                   :: thresh

  call ESMF_StateGet(Routing_State,"Surface elevation",&
       sfcelevField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Surface elevation failed in HYMAP3_updateSWOT")

  call ESMF_FieldGet(sfcelevField,localDE=0,&
       farrayPtr=sfcelev,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Surface elevation failed in HYMAP3_updateSWOT")

  call ESMF_StateGet(Routing_Incr_State,"Surface elevation",&
       sfcelevIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Surface elevation failed in HYMAP3_updateSWOT")

  call ESMF_FieldGet(sfcelevIncrField,localDE=0,&
       farrayPtr=sfcelevIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Surface elevation failed in HYMAP3_updateSWOT")

  nstate = HYMAP3_routing_struc(n)%nseqall*LIS_rc%nensem(n)

  allocate(sfcelevIncr_tmp(nstate))
  sfcelevIncr_tmp = 0.0

  if(HYMAP3_daSWOT_struc(n)%useLocalUpd.eq.1) then

     select case(HYMAP3_daSWOT_struc(n)%localBlendOption)

     case(1)
        ! Only the largest weight seen so far is needed.  This avoids
        ! storing every candidate weight in a two-dimensional array.
        allocate(maxWeight(nstate))
        maxWeight = -huge(1.0)

     case(2)
        allocate(weightedIncrSum(nstate))
        allocate(weightSum(nstate))
        allocate(maxWeight(nstate))
        weightedIncrSum = 0.0
        weightSum = 0.0
        maxWeight = 0.0

     case default
        write(LIS_logunit,*) &
             '[ERR] Invalid HYMAP3 localization blending option: ',&
             HYMAP3_daSWOT_struc(n)%localBlendOption
        write(LIS_logunit,*) &
             '[ERR] Valid options are 1=MAX-AREA, 2=BLEND-AREA'
        call LIS_endrun

     end select
  endif

  thresh = 5.0  ! upper update limit

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m

        ! An anomalous value differential was detected,
        ! and the update was limited to the established threshold.
        if(abs(sfcelevIncr(t)).ge.10.0) then
           sfcelevIncr(t) = 0.0
           cycle
        elseif(abs(sfcelevIncr(t)).gt.thresh) then
           sfcelevIncr(t) = &
                max(-thresh,min(sfcelevIncr(t),thresh))
        endif

        if(HYMAP3_daSWOT_struc(n)%useLocalUpd.eq.1) then
           if(abs(sfcelevIncr(t)).gt.0.0) then
              localWeight = -9999.0

              ix = HYMAP3_routing_struc(n)%seqx(i)
              iy = HYMAP3_routing_struc(n)%seqy(i)

              siteid = int(HYMAP3_daSWOT_struc(n)%sites(ix,iy))

              if(siteid.ge.1 .and. &
                   siteid.le.HYMAP3_daSWOT_struc(n)%nsites) then

                 localWeight(:,:) = &
                      HYMAP3_daSWOT_struc(n)%localWeight(:,:,siteid)

                 c1=max(1,ix-HYMAP3_daSWOT_struc(n)%localupdDX)
                 c2=min(LIS_rc%lnc(n),&
                      ix+HYMAP3_daSWOT_struc(n)%localupdDX)
                 r1=max(1,iy-HYMAP3_daSWOT_struc(n)%localupdDX)
                 r2=min(LIS_rc%lnr(n),&
                      iy+HYMAP3_daSWOT_struc(n)%localupdDX)

                 do r=r1,r2
                    do c=c1,c2
                       i1 = LIS_routing(n)%gindex(c,r)
                       if(i1.gt.0) then
                          t1 = (i1-1)*LIS_rc%nensem(n)+m

                          if(sfcelev(t1).ne.-9999.0 .and. &
                               localWeight(c,r).ne.-9999.0 .and. &
                               localWeight(ix,iy).ne.-9999.0) then

                             pairWeight = localWeight(c,r)

                             select case(&
                                  HYMAP3_daSWOT_struc(n)% &
                                  localBlendOption)

                             case(1)
                                if(pairWeight.gt.maxWeight(t1)) then
                                   maxWeight(t1) = pairWeight
                                   sfcelevIncr_tmp(t1) = &
                                        sfcelevIncr(t)*pairWeight
                                endif

                             case(2)
                                weightedIncrSum(t1) = &
                                     weightedIncrSum(t1) + &
                                     pairWeight*sfcelevIncr(t)

                                weightSum(t1) = weightSum(t1) + &
                                     pairWeight

                                maxWeight(t1) = max(maxWeight(t1),&
                                     pairWeight)

                             end select
                          endif
                       endif
                    enddo
                 enddo

              elseif(siteid.eq.-9999) then
                 ! Preserve a direct increment without a localization entry.
                 sfcelevIncr_tmp(t) = sfcelevIncr(t)

              else
                 write(LIS_logunit,*) &
                      '[ERR] Invalid localization site index: ',siteid
                 write(LIS_logunit,*) &
                      '[ERR] Expected -9999 or a value from 1 to ',&
                      HYMAP3_daSWOT_struc(n)%nsites
                 write(LIS_logunit,*) &
                      '[ERR] Local grid indices ix, iy: ',ix,iy
                 call LIS_endrun
              endif
           endif

        else
           sfcelevIncr_tmp(t) = sfcelevIncr(t)
        endif
     enddo
  enddo

  ! Finalize BLEND-AREA after all source observations have contributed.
  ! A nonzero direct DA increment always takes precedence at its own cell.
  if(HYMAP3_daSWOT_struc(n)%useLocalUpd.eq.1 .and. &
       HYMAP3_daSWOT_struc(n)%localBlendOption.eq.2) then

     do i=1,HYMAP3_routing_struc(n)%nseqall
        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m

           if(abs(sfcelevIncr(t)).gt.0.0) then
              sfcelevIncr_tmp(t) = sfcelevIncr(t)
           elseif(weightSum(t).gt.0.0) then
              sfcelevIncr_tmp(t) = maxWeight(t) * &
                   weightedIncrSum(t)/weightSum(t)
           endif
        enddo
     enddo
  endif

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        if(abs(sfcelevIncr_tmp(t)).gt.0.0) then
           sfcelev(t) = sfcelev(t)+sfcelevIncr_tmp(t)
        endif
     enddo
  enddo

  deallocate(sfcelevIncr_tmp)

  if(allocated(weightedIncrSum)) then
     deallocate(weightedIncrSum)
  endif
  if(allocated(weightSum)) then
     deallocate(weightSum)
  endif
  if(allocated(maxWeight)) then
     deallocate(maxWeight)
  endif

end subroutine HYMAP3_updateSWOT

