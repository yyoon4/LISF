program check_interfaces
  use ESMF
  use LIS_coreMod
  use LIS_snowMod
  use NoahMP50_lsmMod
  use noahmp50_dasnip_Mod
  implicit none
  type(ESMF_State) :: state, increments, observations
  real :: prediction(1,5), reference(31), result(31)
  real, parameter :: changes(2,5) = reshape([ &
       0.0, 0.0, 10.0, 0.05, -10.0, -0.05, &
       -100.0, -0.5, 100.0, 0.5], [2,5])
  integer :: t, i

  allocate(Noahmp50_struc(1))
  allocate(Noahmp50_struc(1)%noahmp50(5))
  allocate(Noahmp50_struc(1)%sldpth(4))
  Noahmp50_struc(1)%sldpth = [0.1, 0.3, 0.6, 1.0]
  Noahmp50_struc(1)%nsnow = 3
  Noahmp50_struc(1)%nsoil = 4
  Noahmp50_struc(1)%ts = 3600.0
  do t = 1, 5
     allocate(Noahmp50_struc(1)%noahmp50(t)%snowice(3))
     allocate(Noahmp50_struc(1)%noahmp50(t)%snowliq(3))
     allocate(Noahmp50_struc(1)%noahmp50(t)%tsno(3))
     allocate(Noahmp50_struc(1)%noahmp50(t)%tslb(4))
     allocate(Noahmp50_struc(1)%noahmp50(t)%smc(4))
     allocate(Noahmp50_struc(1)%noahmp50(t)%sh2o(4))
     allocate(Noahmp50_struc(1)%noahmp50(t)%zss(7))
  end do
  allocate(state%swe%values(5), state%snowdepth%values(5))
  allocate(increments%swe%values(5), increments%snowdepth%values(5))
  allocate(observations%observation%values(1))
  call noahmp50_dasnip_init()

  ! Compare all snow-update outputs with the existing v5.0 update.
  ! This checks reuse, not the independent physical accuracy of that routine.
  do i = 1, size(changes, 2)
     call reset_model(100.0, 0.5)
     call noahmp50_usafsi_update(1, 1, changes(1,i), changes(2,i))
     call snapshot(reference)
     call reset_model(100.0, 0.5)
     call noahmp50_snip_update(1, 1, changes(1,i), changes(2,i))
     call snapshot(result)
     call check(all(abs(result-reference) < 1.e-5), 'snow update reuse')
  end do

  call reset_model(0.0, 0.0)
  call noahmp50_snip_update(1, 1, 10.0, 0.05)
  call check(Noahmp50_struc(1)%noahmp50(1)%isnow < 0, 'new snow layer')
  call check(abs(Noahmp50_struc(1)%noahmp50(1)%sneqv-10.0) < 1.e-5, 'new snow SWE')

  call reset_model(100.0, 0.5)
  call noahmp50_snip_update(1, 1, -100.0, -0.5)
  call check(Noahmp50_struc(1)%noahmp50(1)%isnow == 0, 'complete snow removal')
  call check(Noahmp50_struc(1)%noahmp50(1)%sneqv == 0.0, 'zero SWE after removal')
  call check(Noahmp50_struc(1)%noahmp50(1)%snowh == 0.0, 'zero depth after removal')

  call reset_model(100.0, 0.5)
  call noahmp50_getsnipvars(1, state)
  call check(all(state%swe%values == 100.0), 'SWE in mm')
  call check(all(state%snowdepth%values == 0.5), 'snow depth in m')
  call noahmp50_getsnippred(1, 1, prediction)
  call check(all(prediction == 0.5), 'observation operator in m')
  call noahmp50_scale_snip(1, state)
  call noahmp50_descale_snip(1, state, increments)
  observations%observation%values = 0.55
  call noahmp50_transform_snip(1, observations)
  call check(all(observations%observation%values == 0.55), 'observation units')
  call noahmp50_qc_snipobs(1, 1, observations)
  call noahmp50_map_snip(1, 1, observations, increments)
  call check(all(abs(increments%swe%values-10.0) < 1.e-4), 'mapped SWE increment')
  call noahmp50_updatesnipvars(1, state, increments)
  call noahmp50_qcsnip(1, state)
  call check(all(abs(state%swe%values-110.0) < 1.e-4), 'analysis SWE')
  LIS_rc%snowsrc = 1
  call noahmp50_setsnipvars(1, state)
  call check(abs(LIS_snow_struc(1)%sneqv(1)-550.0) < 1.e-3, 'snow-source tile collection')
  do t = 1, 5
     call check(abs(Noahmp50_struc(1)%noahmp50(t)%sneqv-110.0) < 1.e-4, 'model SWE update')
     call check(abs(sum(Noahmp50_struc(1)%noahmp50(t)%snowice) + &
          sum(Noahmp50_struc(1)%noahmp50(t)%snowliq)-110.0) < 1.e-4, 'layer mass')
  end do
  print *, 'PASS: all 12 SNIP sources compile; interface, units and snow-update checks pass.'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not. condition) then
       print *, 'FAIL: ', label
       error stop 1
    end if
  end subroutine
  subroutine reset_model(swe, depth)
    real, intent(in) :: swe, depth
    integer :: j
    do j = 1, 5
       associate(p => Noahmp50_struc(1)%noahmp50(j))
         p%sneqv = swe
         p%snowh = depth
         p%isnow = 0
         if (depth >= 0.025) p%isnow = -1
         p%snowice = [0.0, 0.0, swe]
         p%snowliq = 0.0
         p%tsno = [0.0, 0.0, 268.0]
         p%tslb = 275.0
         p%sfctmp = 268.0
         p%tv = 268.0
         p%fveg = 0.5
         p%smc = 0.3
         p%sh2o = 0.3
         p%soiltype = 1
         p%param%BEXP = 4.0
         p%param%PSISAT = 0.1
         p%param%SMCMAX = 0.45
         p%zss = [0.0, 0.0, -depth, -depth-0.1, -depth-0.4, -depth-1.0, -depth-2.0]
       end associate
    end do
  end subroutine
  subroutine snapshot(values)
    real, intent(out) :: values(31)
    associate(p => Noahmp50_struc(1)%noahmp50(1))
      values = [real(p%isnow), p%sneqv, p%snowh, p%snowice, &
           p%snowliq, p%tsno, p%tslb, p%zss, p%smc, p%sh2o]
    end associate
  end subroutine
end program
