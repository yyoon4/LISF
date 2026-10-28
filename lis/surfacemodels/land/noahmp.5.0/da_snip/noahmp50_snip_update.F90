!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp50_snip_update
!
! !DESCRIPTION:
!  Apply SNIP analysis increments directly to Noah-MP 5.0 snow states.
!  Scale layer thickness and water mass independently to the requested
!  depth and SWE. Preserve relative layer profiles and ice/liquid fractions;
!  merge/split layers conservatively using Noah-MP geometry thresholds.
!  This is an analysis mapping, not a model time step: do not compact snow,
!  drain liquid water, or apply a second glacier cap here. Those processes
!  belong to the normal model advance (QC has already bounded the targets).
!
! !REVISION HISTORY:
! 19 Oct 2026: Initial SNIP interface.
! 28 Oct 2026: Independent, target-conserving SNIP state update.
!EOP
subroutine noahmp50_snip_update(n, t, dsneqv, dsnowh)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use NoahMP50_lsmMod, only: Noahmp50_struc
  use NoahMP50_snowphys_updateMod, only: Divide, Combo
  use LIS_logMod, only: LIS_logunit, LIS_endrun
  implicit none

  integer, intent(in) :: n, t
  real, intent(in) :: dsneqv, dsnowh
  real, parameter :: tfrz = 273.16, thin_depth = 0.025
  real, parameter :: dzmin(3) = (/0.025, 0.025, 0.1/)
  integer :: nsnow, nsoil, old_isnow, isnow, nl, i, j, k, neighbor, lower, upper
  real :: old_h, old_s, target_h, target_s, sum_h, sum_s, depth, htolerance, stolerance
  real :: dz(3), ice(3), liquid(3), temperature(3)
  real :: fdz(-2:Noahmp50_struc(n)%nsoil), fstc(-2:Noahmp50_struc(n)%nsoil)
  real :: fice(-2:0), fliq(-2:0)
  logical :: merged

  nsnow = Noahmp50_struc(n)%nsnow
  nsoil = Noahmp50_struc(n)%nsoil
  if (nsnow /= 3) call fail('SNIP layer mapping requires nsnow=3')
  old_isnow = Noahmp50_struc(n)%noahmp50(t)%isnow
  old_h = Noahmp50_struc(n)%noahmp50(t)%snowh
  old_s = Noahmp50_struc(n)%noahmp50(t)%sneqv
  if (.not.all(ieee_is_finite((/old_h, old_s, dsnowh, dsneqv/)))) &
       call fail('non-finite snow state or increment')
  target_h = old_h + dsnowh
  target_s = old_s + dsneqv
  if (.not.all(ieee_is_finite((/target_h, target_s/)))) call fail('non-finite analysis target')
  if (old_isnow > 0 .or. old_isnow < -nsnow) call fail('invalid input snow layer count')
  if (old_h < 0.0 .or. old_s < 0.0) call fail('negative input snow total')
  if (target_h < 0.0 .or. target_s < 0.0) call fail('negative analysis target')
  ! A positive thickness with no mass (or vice versa) is not a realizable
  ! snow analysis. Do not silently replace a QC target or invent a density.
  if ((target_h == 0.0) .neqv. (target_s == 0.0)) call fail('inconsistent depth/SWE targets')

  ! A zero update must not perform extra physics or change the vertical state.
  if (dsnowh == 0.0 .and. dsneqv == 0.0) return

  dz = 0.0
  ice = 0.0
  liquid = 0.0
  temperature = 0.0
  fdz = 0.0
  fstc = 0.0
  fice = 0.0
  fliq = 0.0
  isnow = 0

  if (target_h >= thin_depth) then
     nl = -old_isnow
     if (nl > 0) then
        do i = 1, nl
           k = nsnow + old_isnow + i
           if (i == 1) then
              dz(i) = -Noahmp50_struc(n)%noahmp50(t)%zss(k)
           else
              dz(i) = Noahmp50_struc(n)%noahmp50(t)%zss(k-1) - &
                      Noahmp50_struc(n)%noahmp50(t)%zss(k)
           endif
           ice(i) = Noahmp50_struc(n)%noahmp50(t)%snowice(k)
           liquid(i) = Noahmp50_struc(n)%noahmp50(t)%snowliq(k)
           temperature(i) = Noahmp50_struc(n)%noahmp50(t)%tsno(k)
        enddo
        if (.not.all(ieee_is_finite(dz(1:nl))) .or. &
            .not.all(ieee_is_finite(ice(1:nl))) .or. &
            .not.all(ieee_is_finite(liquid(1:nl))) .or. &
            .not.all(ieee_is_finite(temperature(1:nl)))) call fail('non-finite active layer')
        if (any(dz(1:nl) <= 0.0) .or. any(ice(1:nl) < 0.0) .or. &
            any(liquid(1:nl) < 0.0) .or. any(temperature(1:nl) <= 0.0)) &
            call fail('invalid active layer')
        sum_h = sum(dz(1:nl))
        sum_s = sum(ice(1:nl)) + sum(liquid(1:nl))
        if (sum_s <= 0.0 .or. any(ice(1:nl)+liquid(1:nl) <= 0.0)) &
             call fail('active snow layer has no water mass')
        dz(1:nl) = (dz(1:nl)/sum_h)*target_h
        ice(1:nl) = (ice(1:nl)/sum_s)*target_s
        liquid(1:nl) = (liquid(1:nl)/sum_s)*target_s
     else
        ! No explicit profile is available for snow-free or thin-snow input.
        ! As in Noah-MP snowfall initialization, create ice at surface T.
        nl = 1
        dz(1) = target_h
        ice(1) = target_s
        if (.not.ieee_is_finite(Noahmp50_struc(n)%noahmp50(t)%sfctmp)) &
             call fail('non-finite surface temperature for new snow layer')
        temperature(1) = min(tfrz, Noahmp50_struc(n)%noahmp50(t)%sfctmp)
        if (.not.ieee_is_finite(temperature(1)) .or. temperature(1) <= 0.0) &
             call fail('invalid surface temperature for new snow layer')
     endif

     ! Merge thin / nearly ice-free layers without Combine's drainage or
     ! deletion of liquid SWE. The subsequent model step handles phase changes.
     do while (nl > 1)
        merged = .false.
        do i = 1, nl
           if (dz(i) >= dzmin(i) .and. ice(i) > 0.1) cycle
           if (i == 1) then
              neighbor = i + 1
           elseif (i == nl) then
              neighbor = i - 1
           else
              neighbor = i + 1
              if (dz(i-1) < dz(i+1)) neighbor = i - 1
           endif
           upper = min(i, neighbor)
           lower = max(i, neighbor)
           call Combo(Noahmp50_struc(n)%noahmp50(t)%param, &
                dz(lower), liquid(lower), ice(lower), temperature(lower), &
                dz(upper), liquid(upper), ice(upper), temperature(upper))
           do j = upper, nl-1
              dz(j) = dz(j+1)
              ice(j) = ice(j+1)
              liquid(j) = liquid(j+1)
              temperature(j) = temperature(j+1)
           enddo
           nl = nl - 1
           merged = .true.
           exit
        enddo
        if (.not.merged) exit
     enddo

     isnow = -nl
     do i = 1, nl
        j = isnow + i
        fdz(j) = dz(i)
        fice(j) = ice(i)
        fliq(j) = liquid(i)
        fstc(j) = temperature(i)
     enddo
     ! Divide changes geometry but conserves the total ice and liquid masses.
     call Divide(Noahmp50_struc(n)%noahmp50(t)%param, nsnow, nsoil, &
          isnow, fstc, fice, fliq, fdz)

     if (isnow >= 0 .or. isnow < -nsnow) call fail('invalid reconstructed layer count')
     i = isnow + 1
     if (.not.all(ieee_is_finite(fdz(i:0))) .or. &
         .not.all(ieee_is_finite(fice(i:0))) .or. &
         .not.all(ieee_is_finite(fliq(i:0))) .or. &
         .not.all(ieee_is_finite(fstc(i:0)))) call fail('non-finite reconstructed layer')
     if (any(fdz(i:0) <= 0.0) .or. any(fice(i:0) < 0.0) .or. &
         any(fliq(i:0) < 0.0) .or. any(fstc(i:0) <= 0.0)) &
         call fail('invalid reconstructed layer')
     sum_h = sum(fdz(i:0))
     sum_s = sum(fice(i:0)) + sum(fliq(i:0))
     if (sum_s <= 0.0) call fail('reconstructed snow has no water mass')
     ! Only roundoff may remain after conservative merge/split operations.
     htolerance = 8.0*epsilon(target_h)*max(1.0, target_h)
     stolerance = 8.0*epsilon(target_s)*max(1.0, target_s)
     if (abs(sum_h-target_h) > htolerance .or. abs(sum_s-target_s) > stolerance) &
          call fail('layer remapping did not conserve analysis targets')
     fdz(i:0) = (fdz(i:0)/sum_h)*target_h
     fice(i:0) = (fice(i:0)/sum_s)*target_s
     fliq(i:0) = (fliq(i:0)/sum_s)*target_s
  endif

  ! Commit the result only after validation. Inactive output slots stay zero.
  ! Thin snow has aggregate SWE/depth and no explicit layer or liquid profile.
  Noahmp50_struc(n)%noahmp50(t)%isnow = isnow
  Noahmp50_struc(n)%noahmp50(t)%snowh = target_h
  Noahmp50_struc(n)%noahmp50(t)%sneqv = target_s
  Noahmp50_struc(n)%noahmp50(t)%snowice = 0.0
  Noahmp50_struc(n)%noahmp50(t)%snowliq = 0.0
  Noahmp50_struc(n)%noahmp50(t)%tsno = 0.0
  Noahmp50_struc(n)%noahmp50(t)%zss(1:nsnow) = 0.0
  depth = 0.0
  do i = isnow+1, 0
     k = i + nsnow
     depth = depth + fdz(i)
     Noahmp50_struc(n)%noahmp50(t)%snowice(k) = fice(i)
     Noahmp50_struc(n)%noahmp50(t)%snowliq(k) = fliq(i)
     Noahmp50_struc(n)%noahmp50(t)%tsno(k) = fstc(i)
     Noahmp50_struc(n)%noahmp50(t)%zss(k) = -depth
  enddo
  ! Soil thicknesses and temperatures do not change; only their origin shifts.
  if (isnow < 0) depth = target_h
  do i = 1, nsoil
     depth = depth + Noahmp50_struc(n)%sldpth(i)
     Noahmp50_struc(n)%noahmp50(t)%zss(nsnow+i) = -depth
  enddo

contains
  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(LIS_logunit,*) '[ERR] SNIP_UPDATE ', trim(message), ' nest/tile=', n, t
    write(LIS_logunit,*) '[ERR] SNIP_UPDATE old depth/SWE, increments=', &
         Noahmp50_struc(n)%noahmp50(t)%snowh, Noahmp50_struc(n)%noahmp50(t)%sneqv, dsnowh, dsneqv
    call LIS_endrun()
  end subroutine fail
end subroutine noahmp50_snip_update
