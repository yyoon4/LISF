program test_snip_update
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use NoahMP50_lsmMod
  implicit none
  external :: noahmp50_snip_update
  integer :: i, j, k, cases
  real :: h, s, oldh, olds, mass, depth, temp0(3), ice0(3), liq0(3), z0(7)
  real :: soil_t(4), soil_w(4), soil_s(4)
  real, parameter :: targets(21) = (/0., .001, .024, .025, .026, .049, .05, .051, &
       .099, .1, .249, .25, .251, .4, .8, .95, .98, .999, 1., 1.1, 2./)
  character(len=64) :: arg
  cases=0
  call get_command_argument(1,arg)
  if (len_trim(arg)>0) then
     call initial(3, .false.)
     select case(trim(arg))
     case('negative')
        call noahmp50_snip_update(1,1,-301.,-1.1)
     case('inconsistent')
        call noahmp50_snip_update(1,1,-300.,-.6)
     case('nan')
        h=ieee_value(0.,ieee_quiet_nan)
        call noahmp50_snip_update(1,1,0.,h)
     case default
        error stop 2
     end select
     error stop 'Invalid target was accepted'
  endif

  ! Original counterexamples, thresholds, increase and complete snow removal.
  do i=1,size(targets)
     do j=0,1
        call initial(3,j==1)
        call apply_and_check(targets(i),300.*targets(i))
     enddo
  enddo
  ! One/two-layer inputs and layers created from thin/no snow.
  do j=0,2
     do i=1,size(targets)
        call initial(j,.false.)
        call apply_and_check(targets(i),200.*targets(i))
     enddo
  enddo
  call initial(0,.false.)
  Noahmp50_struc(1)%noahmp50(1)%snowh=.01
  Noahmp50_struc(1)%noahmp50(1)%sneqv=2.
  Noahmp50_struc(1)%noahmp50(1)%sfctmp=280.
  call apply_and_check(.3,60.)
  call require(maxval(Noahmp50_struc(1)%noahmp50(1)%tsno)<=273.161,'new snow temperature')

  ! Independent depth/SWE targets (opposite signs and zero increments).
  call initial(3,.true.)
  call apply_and_check(.8,330.)
  call initial(3,.true.)
  call apply_and_check(1.2,240.)
  call initial(3,.true.)
  call apply_and_check(1.,240.)
  call initial(3,.true.)
  call apply_and_check(.8,300.)

  ! Identity operation must not compact, rebuild or alter soil/snow states.
  call initial(3,.true.)
  z0=Noahmp50_struc(1)%noahmp50(1)%zss
  temp0=Noahmp50_struc(1)%noahmp50(1)%tsno
  ice0=Noahmp50_struc(1)%noahmp50(1)%snowice
  liq0=Noahmp50_struc(1)%noahmp50(1)%snowliq
  call apply_and_check(1.,300.)
  call require(all(z0==Noahmp50_struc(1)%noahmp50(1)%zss),'no-op geometry')
  call require(all(temp0==Noahmp50_struc(1)%noahmp50(1)%tsno),'no-op temperature')
  call require(all(ice0==Noahmp50_struc(1)%noahmp50(1)%snowice),'no-op ice')
  call require(all(liq0==Noahmp50_struc(1)%noahmp50(1)%snowliq),'no-op liquid')

  ! A small reduction keeps all three valid layers and their temperatures.
  call initial(3,.true.)
  temp0=Noahmp50_struc(1)%noahmp50(1)%tsno
  ice0=Noahmp50_struc(1)%noahmp50(1)%snowice
  liq0=Noahmp50_struc(1)%noahmp50(1)%snowliq
  call apply_and_check(.98,294.)
  call require(all(temp0==Noahmp50_struc(1)%noahmp50(1)%tsno),'retained temperature profile')
  call require(maxval(abs(Noahmp50_struc(1)%noahmp50(1)%snowice-.98*ice0))<1.e-4,'ice scaling')
  call require(maxval(abs(Noahmp50_struc(1)%noahmp50(1)%snowliq-.98*liq0))<1.e-5,'liquid scaling')

  ! No hidden 5000 mm glacier cap in the analysis mapping.
  call initial(3,.false.)
  call apply_and_check(20.,6000.)

  ! A nearly ice-free layer must merge without draining analysis liquid water.
  call initial(3,.true.)
  Noahmp50_struc(1)%noahmp50(1)%snowice(1)=0.01
  Noahmp50_struc(1)%noahmp50(1)%snowliq(1)=14.99
  call apply_and_check(.4,120.)
  call require(abs(sum(Noahmp50_struc(1)%noahmp50(1)%snowliq)-.4*19.99)<1.e-4,'wet merge liquid')

  ! Deterministic spread of depth/mass ratios and initial layer counts.
  do i=1,1000
     j=mod(i,4)
     call initial(j,mod(i,2)==0)
     h=.0005+real(mod(37*i,2001))/1000.
     s=h*(50.+real(mod(53*i,451)))
     call apply_and_check(h,s)
  enddo
  write(*,'(a,1x,i0)') 'SNIP_UNIT_PASS',cases

contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
       write(*,*) 'SNIP_UNIT_FAIL case=',cases,' ',trim(label)
       error stop 1
    endif
  end subroutine

  subroutine initial(nlayers,wet)
    integer, intent(in) :: nlayers
    logical, intent(in) :: wet
    real :: d(3), m(3), l(3), tt(3), z
    integer :: ii
    d=0.;m=0.;l=0.;tt=0.
    select case(nlayers)
    case(1)
       d(3)=.04;m(3)=12.;tt(3)=260.
    case(2)
       d(2:3)=(/.04,.06/);m(2:3)=(/12.,18./);tt(2:3)=(/255.,260./)
    case(3)
       d=(/.05,.2,.75/);m=(/15.,60.,225./);tt=(/250.,255.,260./)
    end select
    if (wet) then
       l=.1*m
       m=.9*m
       if(nlayers==3) then
          l=(/1.,2.,3./);m=(/14.,58.,222./)
       endif
    endif
    Noahmp50_struc(1)%noahmp50(1)%isnow=-nlayers
    Noahmp50_struc(1)%noahmp50(1)%snowh=sum(d)
    Noahmp50_struc(1)%noahmp50(1)%sneqv=sum(m)+sum(l)
    Noahmp50_struc(1)%noahmp50(1)%snowice=m
    Noahmp50_struc(1)%noahmp50(1)%snowliq=l
    Noahmp50_struc(1)%noahmp50(1)%tsno=tt
    Noahmp50_struc(1)%noahmp50(1)%sfctmp=260.
    Noahmp50_struc(1)%noahmp50(1)%tslb=(/260.,265.,270.,275./)
    Noahmp50_struc(1)%noahmp50(1)%sh2o=.2
    Noahmp50_struc(1)%noahmp50(1)%smc=.3
    z=0.
    do ii=1,3
       z=z+d(ii)
       Noahmp50_struc(1)%noahmp50(1)%zss(ii)=-z
    enddo
    do ii=1,4
       z=z+Noahmp50_struc(1)%sldpth(ii)
       Noahmp50_struc(1)%noahmp50(1)%zss(3+ii)=-z
    enddo
  end subroutine

  subroutine apply_and_check(goal_h,goal_s)
    real, intent(in) :: goal_h, goal_s
    real :: dh, ds, ht, st, expected_z, liq_before, total_before, next_h
    integer :: count_layers, first, ii
    cases=cases+1
    oldh=Noahmp50_struc(1)%noahmp50(1)%snowh
    olds=Noahmp50_struc(1)%noahmp50(1)%sneqv
    soil_t=Noahmp50_struc(1)%noahmp50(1)%tslb
    soil_w=Noahmp50_struc(1)%noahmp50(1)%sh2o
    soil_s=Noahmp50_struc(1)%noahmp50(1)%smc
    liq_before=sum(Noahmp50_struc(1)%noahmp50(1)%snowliq)
    total_before=sum(Noahmp50_struc(1)%noahmp50(1)%snowice)+liq_before
    dh=goal_h-oldh;ds=goal_s-olds
    ! Threshold classification follows the representable increment target.
    next_h=oldh+dh
    call noahmp50_snip_update(1,1,ds,dh)
    ht=max(1.e-6,8.*epsilon(h)*max(1.,oldh,goal_h))
    st=max(1.e-4,8.*epsilon(s)*max(1.,olds,goal_s))
    call require(abs(Noahmp50_struc(1)%noahmp50(1)%snowh-goal_h)<=ht,'depth target')
    call require(abs(Noahmp50_struc(1)%noahmp50(1)%sneqv-goal_s)<=st,'SWE target')
    count_layers=-Noahmp50_struc(1)%noahmp50(1)%isnow
    call require(count_layers>=0.and.count_layers<=3,'layer count bounds')
    if(next_h>=.025) call require(count_layers>0,'deep snow has active layers')
    if(next_h<.025) call require(count_layers==0,'thin snow aggregate representation')
    first=4-count_layers
    if(count_layers>0) then
       mass=sum(Noahmp50_struc(1)%noahmp50(1)%snowice)+sum(Noahmp50_struc(1)%noahmp50(1)%snowliq)
       call require(abs(mass-goal_s)<=st,'layer mass target')
       call require(abs(-Noahmp50_struc(1)%noahmp50(1)%zss(3)-goal_h)<=ht,'layer depth target')
       if(total_before>0.) call require(abs(sum(Noahmp50_struc(1)%noahmp50(1)%snowliq)- &
            liq_before/total_before*goal_s)<=st,'liquid fraction conserved')
    endif
    depth=0.
    do ii=1,3
       call require(ieee_is_finite(Noahmp50_struc(1)%noahmp50(1)%tsno(ii)),'finite temperature')
       if(ii<first) then
          call require(Noahmp50_struc(1)%noahmp50(1)%snowice(ii)==0.,'inactive ice')
          call require(Noahmp50_struc(1)%noahmp50(1)%snowliq(ii)==0.,'inactive liquid')
          call require(Noahmp50_struc(1)%noahmp50(1)%tsno(ii)==0.,'inactive temperature')
          call require(Noahmp50_struc(1)%noahmp50(1)%zss(ii)==0.,'inactive depth')
       else
          call require(Noahmp50_struc(1)%noahmp50(1)%snowice(ii)>=0.,'nonnegative ice')
          call require(Noahmp50_struc(1)%noahmp50(1)%snowliq(ii)>=0.,'nonnegative liquid')
          call require(Noahmp50_struc(1)%noahmp50(1)%tsno(ii)>150.,'active temperature')
          call require(-Noahmp50_struc(1)%noahmp50(1)%zss(ii)>depth,'positive thickness')
          depth=-Noahmp50_struc(1)%noahmp50(1)%zss(ii)
       endif
    enddo
    expected_z=0.
    if(count_layers>0) expected_z=Noahmp50_struc(1)%noahmp50(1)%snowh
    do ii=1,4
       expected_z=expected_z+Noahmp50_struc(1)%sldpth(ii)
       call require(abs(-Noahmp50_struc(1)%noahmp50(1)%zss(3+ii)-expected_z)<=ht,'soil geometry')
    enddo
    call require(all(soil_t==Noahmp50_struc(1)%noahmp50(1)%tslb),'soil temperature unchanged')
    call require(all(soil_w==Noahmp50_struc(1)%noahmp50(1)%sh2o),'soil liquid unchanged')
    call require(all(soil_s==Noahmp50_struc(1)%noahmp50(1)%smc),'soil moisture unchanged')
  end subroutine
end program
