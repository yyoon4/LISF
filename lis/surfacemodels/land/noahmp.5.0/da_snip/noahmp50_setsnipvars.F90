!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp50_setsnipvars
! \label{noahmp50_setsnipvars}
!
! !REVISION HISTORY:
! 19 Oct 2026: Port SNIP DA to NoahMP 5.0 using existing DA interfaces.
! 18 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
subroutine noahmp50_setsnipvars(n, LSM_State)

! !USES:
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_logMod, only : LIS_verify, LIS_endrun, LIS_logunit
  use LIS_snowMod, only : LIS_snow_struc
  use noahmp50_lsmMod

  ! Defaults
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  This routine assigns the snow prognostic variables to NoahMP 5.0
!  model space. The state vector consists of total SWE and snow depth.
!  This routine also updates other model prognostics (snice, snliq,
!  snow thickness, snow temperature) based on the update.
!
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh
  integer                :: t
  integer                :: status
  integer                :: ncount(LIS_rc%ngrid(n))
  integer                :: tid, gid
  integer, save          :: check_call = 0
  integer                :: changed, reduced, layer_changes, failed, old_isnow
  real                   :: old_h, old_s, dh, ds, dm, dz, htol, stol
  real                   :: max_dh, max_ds, max_dm, max_dz

  integer                :: trace_normalized
  real                   :: trace_depth, max_trace_depth
  logical                :: normalized

  external :: noahmp50_snip_update, noahmp50_snip_qc_target

  call ESMF_StateGet(LSM_State, "SWE", sweField, rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State, "Snowdepth", snodField, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField, localDE=0, farrayPtr=swe, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField, localDE=0, farrayPtr=snod, rc=status)
  call LIS_verify(status)

  check_call = check_call + 1
  trace_normalized = 0
  max_trace_depth = 0.0
  changed = 0
  reduced = 0
  layer_changes = 0
  failed = 0
  max_dh = 0.0
  max_ds = 0.0
  max_dm = 0.0
  max_dz = 0.0
  if (check_call == 1) then
     write(LIS_logunit,'(a)') 'SNIP_TARGET_CHECK_COLUMNS nest call patches changed reduced '// &
          'layers_changed failed max_depth_error_m max_SWE_error_kgm2 '// &
          'max_layer_mass_error_kgm2 max_layer_depth_error_m'
  endif

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     old_h = noahmp50_struc(n)%noahmp50(t)%snowh
     old_s = noahmp50_struc(n)%noahmp50(t)%sneqv
     old_isnow = noahmp50_struc(n)%noahmp50(t)%isnow

     ! QC the state targets before computing increments. Do not change
     ! positive SWE, existing snow, or depths at/above the trace threshold.
     trace_depth = snod(t)
     call noahmp50_snip_qc_target(old_h, old_s, snod(t), swe(t), normalized)
     if (normalized) then
        trace_normalized = trace_normalized + 1
        max_trace_depth = max(max_trace_depth, trace_depth)
        if (trace_normalized <= 3) write(LIS_logunit,*) &
             'SNIP_TRACE_SNOW_QC_SAMPLE nest/tile/col/row=', n, t, &
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col, &
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row, &
             ' raw_depth_m=', trace_depth, ' normalized_depth_m=0 SWE=0'
     endif

     dsneqv = swe(t) - noahmp50_struc(n)%noahmp50(t)%sneqv   !in mm
     dsnowh = snod(t) - noahmp50_struc(n)%noahmp50(t)%snowh  !in m

     ! update
     call noahmp50_snip_update(n, t, dsneqv, dsnowh)

     if (dsneqv /= 0.0 .or. dsnowh /= 0.0) changed = changed + 1
     if (dsnowh < 0.0) reduced = reduced + 1
     if (old_isnow /= noahmp50_struc(n)%noahmp50(t)%isnow) layer_changes = layer_changes + 1
     dh = abs(noahmp50_struc(n)%noahmp50(t)%snowh - snod(t))
     ds = abs(noahmp50_struc(n)%noahmp50(t)%sneqv - swe(t))
     dm = 0.0
     dz = 0.0
     if (noahmp50_struc(n)%noahmp50(t)%isnow < 0) then
        dm = abs(sum(noahmp50_struc(n)%noahmp50(t)%snowice) + &
             sum(noahmp50_struc(n)%noahmp50(t)%snowliq) - swe(t))
        dz = abs(-noahmp50_struc(n)%noahmp50(t)%zss(noahmp50_struc(n)%nsnow) - snod(t))
     endif
     ! The increment interface subtracts and re-adds single precision totals.
     ! Include the old magnitude in tolerances to allow cancellation roundoff.
     htol = max(1.e-6, 8.0*epsilon(dh)*max(1.0, abs(old_h), abs(snod(t))))
     stol = max(1.e-4, 8.0*epsilon(ds)*max(1.0, abs(old_s), abs(swe(t))))
     max_dh = max(max_dh, dh)
     max_ds = max(max_ds, ds)
     max_dm = max(max_dm, dm)
     max_dz = max(max_dz, dz)
     if (.not.all(ieee_is_finite((/dh, ds, dm, dz/))) .or. &
         dh > htol .or. ds > stol .or. dm > stol .or. dz > htol) then
        failed = failed + 1
        if (failed <= 5) write(LIS_logunit,*) '[ERR] SNIP_TARGET_FAIL nest/tile/col/row=', &
             n, t, LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col, &
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row, &
             ' target depth/SWE=', snod(t), swe(t), ' errors=', dh, ds, dm, dz
     endif

  enddo

  write(LIS_logunit,'(a,3(1x,i0),1x,es16.8)') 'SNIP_TRACE_SNOW_QC', &
       n, check_call, trace_normalized, max_trace_depth
  write(LIS_logunit,'(a,7(1x,i0),4(1x,es16.8))') 'SNIP_TARGET_CHECK', &
       n, check_call, LIS_rc%npatch(n,LIS_rc%lsm_index), changed, reduced, layer_changes, failed, &
       max_dh, max_ds, max_dm, max_dz
  if (failed > 0) call LIS_endrun()

  if (LIS_rc%snowsrc(n) .gt. 0) then

     ncount = 0 ! Number of tiles per grid id (over land)
     LIS_snow_struc(n)%snowdepth = 0 ! At grid points
     LIS_snow_struc(n)%sneqv = 0     ! At tiles

     ! Collect SWE at tiles
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
        LIS_snow_struc(n)%sneqv(tid) = LIS_snow_struc(n)%sneqv(tid) + &
             noahmp50_struc(n)%noahmp50(t)%sneqv
     end do

     ! Collect mean snow depth at grid points
     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        LIS_snow_struc(n)%snowdepth(gid) = &
             LIS_snow_struc(n)%snowdepth(gid) + &
             noahmp50_struc(n)%noahmp50(t)%snowh
        ncount(gid) = ncount(gid) + 1
     end do
     do t = 1, LIS_rc%ngrid(n)
        if (ncount(t).gt.0) then
           LIS_snow_struc(n)%snowdepth(t) = &
                LIS_snow_struc(n)%snowdepth(t) / ncount(t)
        else
           LIS_snow_struc(n)%snowdepth(t) = 0.0
        endif
     end do
  end if

end subroutine noahmp50_setsnipvars

