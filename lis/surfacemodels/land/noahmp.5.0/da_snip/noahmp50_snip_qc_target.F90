! Normalize only trace, massless targets over initially snow-free ground.
! The strict 2.54E-3 m boundary matches noahmp50_map_snip's trace threshold.
! This routine runs before increments and target-conservation diagnostics.
subroutine noahmp50_snip_qc_target(old_h, old_s, target_h, target_s, normalized)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real, intent(in) :: old_h, old_s, target_s
  real, intent(inout) :: target_h
  logical, intent(out) :: normalized
  real, parameter :: trace_depth = 2.54e-3

  normalized = .false.
  ! Leave invalid values to the updater's existing error checks.
  if (.not.all(ieee_is_finite((/old_h, old_s, target_h, target_s/)))) return
  if (old_h /= 0.0 .or. old_s /= 0.0 .or. target_s /= 0.0) return
  if (target_h > 0.0 .and. target_h < trace_depth) then
     target_h = 0.0
     normalized = .true.
  endif
end subroutine noahmp50_snip_qc_target
