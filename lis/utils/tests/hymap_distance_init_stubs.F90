! Synthetic API stubs. No operational NetCDF file or ESMF library is used.
module fixture
  use,intrinsic :: ieee_arithmetic
  implicit none
  integer :: scenario=0,nest_count=1,cursor=0
  character(160) :: label=''
  integer :: opens=0,closes=0,weight_reads=0,distance_reads=0,distance_lookups=0
contains
  real function base_weight(x,y,s)
    integer,intent(in) :: x,y,s
    base_weight=0.8
    if(x==1.and.y==1.and.s==1) base_weight=-9999.
    if(x==3.and.y==2.and.s==1) then
       select case(scenario)
       case(0,1,10)
          base_weight=ieee_value(0.,ieee_quiet_nan)
       case(8)
          base_weight=-0.1
       case(9)
          base_weight=1.1
       case(26)
          base_weight=ieee_value(0.,ieee_positive_inf)
       end select
    endif
  end function
  real(kind=8) function network_distance(x,y,s)
    integer,intent(in) :: x,y,s
    network_distance=25.0_8*mod(x+2*y+3*s,6)
    if(x==1.and.y==1.and.s==1) network_distance=-9999.0_8
    if(x==3.and.y==2.and.s==1) network_distance=ieee_value(0.0_8,ieee_quiet_nan)
    if(x==4.and.y==2.and.s==1) network_distance=-1.0_8
    if(scenario==23.and.x==3.and.y==2.and.s==1) network_distance=ieee_value(0.0_8,ieee_positive_inf)
  end function
end module

module ESMF
  use fixture
  implicit none
  interface ESMF_ConfigGetAttribute
     module procedure get_integer,get_real,get_character
  end interface
contains
  subroutine ESMF_ConfigFindLabel(config,text,rc)
    integer,intent(in) :: config
    character(*),intent(in) :: text
    integer,intent(out) :: rc
    label=text;cursor=0;rc=0
    if(scenario==0.and.(index(text,'weighting method:')>0.or.index(text,'distance cutoff')>0)) rc=1
    if(scenario==25.and.index(text,'distance cutoff')>0) rc=1
  end subroutine
  subroutine get_integer(config,value,rc)
    integer,intent(in) :: config
    integer,intent(out) :: value,rc
    cursor=cursor+1;rc=0
    if(index(label,'use localization')>0) then
       value=1
       if(scenario==12.or.(scenario==13.and.cursor==2)) value=0
    else
       value=2
    endif
  end subroutine
  subroutine get_real(config,value,rc)
    integer,intent(in) :: config
    real,intent(out) :: value
    integer,intent(out) :: rc
    cursor=cursor+1;rc=0;value=100.
    if(scenario==4.or.scenario==12.or.(scenario==13.and.cursor==2)) value=0.
    if(scenario==14) value=ieee_value(0.,ieee_quiet_nan)
    if(scenario==15) value=ieee_value(0.,ieee_positive_inf)
    if(scenario==22) value=-1.
    if(scenario==24.and.cursor==1) value=0.
  end subroutine
  subroutine get_character(config,value,rc)
    integer,intent(in) :: config
    character(*),intent(out) :: value
    integer,intent(out) :: rc
    cursor=cursor+1;rc=0
    if(index(label,'weighting method:')>0) then
       value='area-distance-max'
       if(scenario==1) value='area-max'
       if(scenario==3) value='bad-method'
       if(scenario==24.and.cursor==1) value='area-max'
    else
       value='synthetic_localization_map.nc'
       if(scenario==12.or.(scenario==13.and.cursor==2)) value='no-such-disabled-map.nc'
    endif
  end subroutine
end module

module LIS_constantsMod
  integer,parameter :: LIS_CONST_PATH_LEN=512
end module
module LIS_coreMod
  implicit none
  type rc_type
     integer :: nnest=1
     integer :: gnc(2)=5,gnr(2)=4,lnc(2)=2,lnr(2)=3
  end type
  type(rc_type) :: LIS_rc
  integer :: LIS_config=0,LIS_localPet=0
  integer :: LIS_ews_halo_ind(2,1)=3,LIS_ewe_halo_ind(2,1)=4
  integer :: LIS_nss_halo_ind(2,1)=2,LIS_nse_halo_ind(2,1)=4
end module
module LIS_logMod
  implicit none
  integer,parameter :: LIS_logunit=6
contains
  subroutine LIS_verify(status,message)
    integer,intent(in) :: status
    character(*),intent(in) :: message
    if(status/=0) then
       print *,'EXPECTED INIT ERROR: ',message
       stop 91
    endif
  end subroutine
  subroutine LIS_endrun()
    stop 92
  end subroutine
end module
module LIS_dataAssimMod
end module
module HYMAP3_routingMod
end module
module HYMAP3_initMod
end module

module netcdf
  use fixture
  implicit none
  integer,parameter :: NF90_NOWRITE=0
  interface nf90_get_var
     module procedure get_real3,get_real2,get_double2
  end interface
  interface nf90_get_att
     module procedure get_attr_char,get_attr_double
  end interface
contains
  integer function nf90_open(path,mode,ncid) result(status)
    character(*),intent(in) :: path
    integer,intent(in) :: mode
    integer,intent(out) :: ncid
    opens=opens+1;ncid=1;status=0
  end function
  integer function nf90_close(ncid) result(status)
    integer,intent(in) :: ncid
    closes=closes+1;status=0
  end function
  integer function nf90_inq_dimid(ncid,name,id) result(status)
    integer,intent(in) :: ncid
    character(*),intent(in) :: name
    integer,intent(out) :: id
    select case(name)
    case('east_west');id=1
    case('north_south');id=2
    case('count');id=3
    case default;stop 80
    end select
    status=0
  end function
  integer function nf90_inquire_dimension(ncid,id,len) result(status)
    integer,intent(in) :: ncid,id
    integer,intent(out) :: len
    select case(id)
    case(1);len=5
    case(2);len=4
    case(3);len=2
    case default;stop 81
    end select
    if(scenario==18.and.id==1) len=6
    status=0
  end function
  integer function nf90_inq_varid(ncid,name,id) result(status)
    integer,intent(in) :: ncid
    character(*),intent(in) :: name
    integer,intent(out) :: id
    status=0
    select case(name)
    case('weight');id=1
    case('sites');id=2
    case('network_distance_km')
       id=3;distance_lookups=distance_lookups+1
       if(scenario==5) status=1
    case default;stop 82
    end select
  end function
  integer function nf90_inquire_variable(ncid,id,ndims,dimids) result(status)
    integer,intent(in) :: ncid,id
    integer,intent(out),optional :: ndims,dimids(:)
    if(present(ndims)) then
       ndims=3
       if(scenario==17.and.id==3) ndims=2
    endif
    if(present(dimids)) then
       dimids=[1,2,3]
       if(scenario==7.and.id==3) dimids=[2,1,3]
       if(scenario==19.and.id==1) dimids=[2,1,3]
    endif
    status=0
  end function
  integer function get_attr_char(ncid,id,name,value) result(status)
    integer,intent(in) :: ncid,id
    character(*),intent(in) :: name
    character(*),intent(out) :: value
    value='km';status=0
    if(scenario==6) value='m'
    if(scenario==20) status=1
  end function
  integer function get_attr_double(ncid,id,name,value) result(status)
    integer,intent(in) :: ncid,id
    character(*),intent(in) :: name
    real(kind=8),intent(out) :: value
    value=-9999.0_8;status=0
    if(scenario==16) value=-999.0_8
    if(scenario==21) status=1
  end function
  integer function get_real3(ncid,id,values,start,count) result(status)
    integer,intent(in) :: ncid,id
    real,intent(out) :: values(:,:,:)
    integer,intent(in),optional :: start(:),count(:)
    integer :: c,r,s,x0,y0,s0
    if(id/=1) stop 83
    weight_reads=weight_reads+1
    x0=1;y0=1;s0=1
    if(present(start)) then
       x0=start(1);y0=start(2);s0=start(3)
       if(any(count/=[size(values,1),size(values,2),size(values,3)])) stop 84
    endif
    do s=1,size(values,3)
       do r=1,size(values,2)
          do c=1,size(values,1)
             values(c,r,s)=base_weight(x0+c-1,y0+r-1,s0+s-1)
          enddo
       enddo
    enddo
    status=0
  end function
  integer function get_real2(ncid,id,values,start,count) result(status)
    integer,intent(in) :: ncid,id,start(:),count(:)
    real,intent(out) :: values(:,:)
    if(id/=2.or.size(start)/=2.or.any(count/=shape(values))) stop 85
    values=1.0;status=0
  end function
  integer function get_double2(ncid,id,values,start,count) result(status)
    integer,intent(in) :: ncid,id,start(:),count(:)
    real(kind=8),intent(out) :: values(:,:)
    integer :: c,r
    if(id/=3.or.size(start)/=3.or.count(3)/=1) stop 86
    if(any(count(1:2)/=shape(values))) stop 87
    distance_reads=distance_reads+1
    do r=1,size(values,2)
       do c=1,size(values,1)
          values(c,r)=network_distance(start(1)+c-1,start(2)+r-1,start(3))
       enddo
    enddo
    status=0
  end function
end module

! INSERT_MODULES_HERE

program init_regression
  use fixture
  use LIS_coreMod
  use HYMAP3_daSWOT_Mod
  use HYMAP3_daSWOT_Mod_baseline, only: baseline_init=>HYMAP3_daSWOT_init,baseline=>HYMAP3_daSWOT_struc
  implicit none
  character(20) :: argument
  integer :: c,r,s,x0,y0,x,y,native_checks,old_reads,old_opens,tn,expected_opens
  real :: expected,w
  real(kind=8) :: d
  call get_command_argument(1,argument)
  read(argument,*) scenario
  if(scenario==13.or.scenario==24) LIS_rc%nnest=2
  if(scenario==0.or.scenario==1) then
     call baseline_init(1)
     if(distance_reads/=0.or.distance_lookups/=0) stop 70
     opens=0;closes=0;weight_reads=0;cursor=0
  endif
  call HYMAP3_daSWOT_init(1)
  if(scenario==12) then
     if(opens/=0.or.distance_reads/=0.or.allocated(HYMAP3_daSWOT_struc(1)%localWeight)) stop 71
     print *,'PASS INIT: disabled distance mode performs no file I/O'
     stop
  endif
  expected_opens=1
  if(scenario==24) expected_opens=2
  if(opens/=expected_opens.or.closes/=expected_opens.or.weight_reads/=expected_opens) stop 72
  if(scenario==13) then
     if(allocated(HYMAP3_daSWOT_struc(2)%localWeight)) stop 73
  endif
  if(scenario==0.or.scenario==1) then
     if(distance_reads/=0.or.distance_lookups/=0) stop 74
     if(any(shape(HYMAP3_daSWOT_struc(1)%localWeight)/=shape(baseline(1)%localWeight))) stop 75
     if(any(transfer(HYMAP3_daSWOT_struc(1)%localWeight,[0],size(baseline(1)%localWeight))/= &
            transfer(baseline(1)%localWeight,[0],size(baseline(1)%localWeight)))) stop 76
     print *,'PASS INIT: legacy weights bitwise equal; no distance I/O'
     stop
  endif
  if(distance_reads/=2.or.distance_lookups/=1) stop 77
  x0=SLICE_X_START;y0=SLICE_Y_START
  tn=1
  if(scenario==24) then
     tn=2
     if(trim(HYMAP3_daSWOT_struc(1)%weightingMethod)/='area-max') stop 68
     do s=1,2
        do r=1,size(HYMAP3_daSWOT_struc(1)%localWeight,2)
           do c=1,size(HYMAP3_daSWOT_struc(1)%localWeight,1)
              if(HYMAP3_daSWOT_struc(1)%localWeight(c,r,s)/=base_weight(x0+c-1,y0+r-1,s)) stop 69
           enddo
        enddo
     enddo
  endif
  native_checks=0
  do s=1,2
     do r=1,size(HYMAP3_daSWOT_struc(tn)%localWeight,2)
        do c=1,size(HYMAP3_daSWOT_struc(tn)%localWeight,1)
           x=x0+c-1;y=y0+r-1
           w=base_weight(x,y,s);d=network_distance(x,y,s)
           if(w==-9999.) then
              expected=-9999.
           elseif(.not.ieee_is_finite(d)) then
              expected=-9999.
           elseif(d<0._8.or.d>=100._8) then
              expected=-9999.
           else
              expected=real(dble(w)*0.5_8*(1._8+cos(acos(-1._8)*d/100._8)))
           endif
           if(HYMAP3_daSWOT_struc(tn)%localWeight(c,r,s)/=expected) stop 78
           native_checks=native_checks+1
        enddo
     enddo
  enddo
  old_reads=distance_reads;old_opens=opens
  call HYMAP3_daSWOT_init(1)
  if(old_reads/=distance_reads.or.old_opens/=opens) stop 79
  print *,'PASS INIT: tapered offset/edge values, one-slice I/O, single application; checks=',native_checks
end program
