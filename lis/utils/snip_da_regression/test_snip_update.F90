program test_snip_update
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use NoahMP50_lsmMod
  implicit none
  external :: noahmp50_snip_update, noahmp50_usafsi_update
  integer :: i,j,k,cases,legacy_cases
  real :: h,s,a,b
  real, parameter :: fractions(9) = (/0.,.1,.25,.4,.75,.8,.95,.98,.999/)
  real, parameter :: reported(3) = (/2.0394144e-4,1.87592127e-3,7.3130843e-3/)
  real :: temp0(3),ice0(3),liq0(3)
  cases=0;legacy_cases=0
  ! Compare every non-reduction path with the actual original updater.
  do k=0,1
     Noahmp50_struc(1)%ts=1800.*real(k)
     do j=0,1
        do i=0,3
           call initial(i,j==1)
           call compare_legacy(0.,0.)
           call compare_legacy(3.,.01)
           call compare_legacy(0.,.01)
           call compare_legacy(3.,0.)
           call compare_legacy(-3.,.01)
           call compare_legacy(3.,-.01)
           call compare_legacy(0.,-.01)
           call compare_legacy(-3.,0.)
           call compare_legacy(-400.,-2.)
        enddo
     enddo
  enddo
  ! Exact reported inputs: no abort, no new SWE or depth.
  do i=1,size(reported)
     call initial(0,.false.)
     call compare_legacy(0.,reported(i))
     a=0.;b=reported(i)
     call noahmp50_snip_update(1,1,a,b)
     call require(Noahmp50_struc(1)%noahmp50(1)%snowh==0.,'reported depth ignored as before')
     call require(Noahmp50_struc(1)%noahmp50(1)%sneqv==0.,'no invented SWE')
  enddo
  ! Dry dt=0 isolates mapping; dt=1800 retains the old compaction.
  do k=0,1
     Noahmp50_struc(1)%ts=1800.*real(k)
     do j=1,3
        do i=1,size(fractions)
           call initial(j,.false.)
           h=Noahmp50_struc(1)%noahmp50(1)%snowh*fractions(i)
           s=Noahmp50_struc(1)%noahmp50(1)%sneqv*fractions(i)
           call reduction(h,s,k==0,.false.)
        enddo
     enddo
  enddo
  Noahmp50_struc(1)%ts=0.
  call initial(3,.false.)
  call reduction(.8,210.,.true.,.false.)
  call initial(3,.true.)
  temp0=Noahmp50_struc(1)%noahmp50(1)%tsno
  ice0=Noahmp50_struc(1)%noahmp50(1)%snowice
  liq0=Noahmp50_struc(1)%noahmp50(1)%snowliq
  call reduction(.98,294.,.true.,.false.)
  call require(all(temp0==Noahmp50_struc(1)%noahmp50(1)%tsno),'retained temperature')
  call require(maxval(abs(.98*ice0-Noahmp50_struc(1)%noahmp50(1)%snowice))<1.e-4,'ice scaling')
  call require(maxval(abs(.98*liq0-Noahmp50_struc(1)%noahmp50(1)%snowliq))<1.e-5,'liquid scaling')
  call initial(0,.false.)
  Noahmp50_struc(1)%noahmp50(1)%snowh=.01
  Noahmp50_struc(1)%noahmp50(1)%sneqv=3.
  call reduction(.005,1.5,.true.,.false.)
  ! Wet thin snow retains Combine's existing liquid drainage convention.
  Noahmp50_struc(1)%ts=1800.
  call initial(3,.true.)
  call reduction(.01,3.,.false.,.true.)
  call require(abs(Noahmp50_struc(1)%noahmp50(1)%sneqv-2.94)<1.e-4,'legacy thin-snow ice total')
  write(*,'(a,2(1x,i0))') 'SNIP_UNIT_PASS',legacy_cases,cases
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
       write(*,*) 'SNIP_UNIT_FAIL ',trim(label)
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

  subroutine compare_legacy(delta_s,delta_h)
    real, intent(in) :: delta_s,delta_h
    real :: aa,bb
    type(tile_type) :: before,expected,actual
    before=Noahmp50_struc(1)%noahmp50(1)
    aa=delta_s;bb=delta_h
    call noahmp50_usafsi_update(1,1,aa,bb)
    expected=Noahmp50_struc(1)%noahmp50(1)
    Noahmp50_struc(1)%noahmp50(1)=before
    aa=delta_s;bb=delta_h
    call noahmp50_snip_update(1,1,aa,bb)
    actual=Noahmp50_struc(1)%noahmp50(1)
    call require(actual%isnow==expected%isnow,'legacy layer count')
    call require(actual%snowh==expected%snowh,'legacy depth')
    call require(actual%sneqv==expected%sneqv,'legacy SWE')
    call require(all(actual%snowice==expected%snowice),'legacy ice')
    call require(all(actual%snowliq==expected%snowliq),'legacy liquid')
    call require(all(actual%tsno==expected%tsno),'legacy snow temperature')
    call require(all(actual%zss==expected%zss),'legacy coordinates')
    call require(all(actual%tslb==expected%tslb),'legacy soil temperature')
    call require(all(actual%sh2o==expected%sh2o),'legacy soil water')
    call require(all(actual%smc==expected%smc),'legacy soil moisture')
    Noahmp50_struc(1)%noahmp50(1)=before
    legacy_cases=legacy_cases+1
  end subroutine
  subroutine reduction(goal_h,goal_s,check_depth,allow_drainage)
    real, intent(in) :: goal_h,goal_s
    logical, intent(in) :: check_depth,allow_drainage
    real :: aa,bb,ht,st,depth
    integer :: first,ii
    ht=max(1.e-6,8.*epsilon(goal_h)*max(1.,Noahmp50_struc(1)%noahmp50(1)%snowh))
    st=max(1.e-4,8.*epsilon(goal_s)*max(1.,Noahmp50_struc(1)%noahmp50(1)%sneqv))
    aa=goal_s-Noahmp50_struc(1)%noahmp50(1)%sneqv
    bb=goal_h-Noahmp50_struc(1)%noahmp50(1)%snowh
    call noahmp50_snip_update(1,1,aa,bb)
    associate(q=>Noahmp50_struc(1)%noahmp50(1))
      call require(ieee_is_finite(q%snowh).and.ieee_is_finite(q%sneqv),'finite totals')
      call require(q%snowh>=0..and.q%sneqv>=0.,'nonnegative totals')
      if(check_depth) call require(abs(q%snowh-goal_h)<=ht,'dry dt=0 target depth')
      if(.not.allow_drainage) call require(abs(q%sneqv-goal_s)<=st,'target SWE')
      call require(q%isnow>=-3.and.q%isnow<=0,'layer count bounds')
      first=4+q%isnow
      depth=0.
      do ii=1,3
         if(ii<first) then
            call require(q%snowice(ii)==0..and.q%snowliq(ii)==0.,'inactive mass')
            call require(q%tsno(ii)==0..and.q%zss(ii)==0.,'inactive temperature/depth')
         else
            call require(q%snowice(ii)>=0..and.q%snowliq(ii)>=0.,'active mass')
            call require(ieee_is_finite(q%tsno(ii)),'finite temperature')
            call require(q%tsno(ii)>150..and.q%tsno(ii)<=273.17,'active temperature range')
            call require(-q%zss(ii)>depth,'positive layer thickness')
            depth=-q%zss(ii)
         endif
      enddo
      if(q%isnow<0) then
         call require(abs(depth-q%snowh)<=ht,'depth agrees with layers')
         call require(abs(sum(q%snowice)+sum(q%snowliq)-q%sneqv)<=st,'SWE agrees with layers')
      endif
    end associate
    cases=cases+1
  end subroutine
end program
