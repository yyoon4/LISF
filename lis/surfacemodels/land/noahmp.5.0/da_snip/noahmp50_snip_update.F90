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
! \label{noahmp50_snip_update}
!
! !REVISION HISTORY:
! 19 Oct 2026: Reuse the NoahMP 5.0 snow DA update for SNIP.
!
! !INTERFACE:
subroutine noahmp50_snip_update(n, t, dsneqv, dsnowh)

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: t
  real :: dsneqv ! SWE increment (mm)
  real :: dsnowh ! snow-depth increment (m)
!
! !DESCRIPTION:
!  SNIP uses the same snow state and units as the NoahMP 5.0 USAFSI
!  interface. Reuse its update so snow-layer physics and per-tile
!  parameters remain consistent with the existing NoahMP 5.0 snow DA.
!  The da_usafsi build path is included with NoahMP 5.0 virtual_da,
!  independently of whether the USAFSI observation reader is enabled.
!EOP
  external :: noahmp50_usafsi_update

  call noahmp50_usafsi_update(n, t, dsneqv, dsnowh)

end subroutine noahmp50_snip_update
