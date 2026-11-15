! Synthetic only. The runner inserts the actual baseline and candidate routine.
module ESMF
  implicit none
  type ESMF_State
     real, pointer :: values(:) => null()
  end type
  type ESMF_Field
     real, pointer :: values(:) => null()
  end type
contains
  subroutine ESMF_StateGet(state, name, field, rc)
    type(ESMF_State), intent(in) :: state
    character(*), intent(in) :: name
    type(ESMF_Field), intent(out) :: field
    integer, intent(out) :: rc
    field%values => state%values
    rc=0
  end subroutine
  subroutine ESMF_FieldGet(field, localDE, farrayPtr, rc)
    type(ESMF_Field), intent(in) :: field
    integer, intent(in) :: localDE
    real, pointer :: farrayPtr(:)
    integer, intent(out) :: rc
    farrayPtr => field%values
    rc=0
  end subroutine
end module

module LIS_coreMod
  implicit none
  type rc_type
     integer, allocatable :: lnc(:),lnr(:),nensem(:)
  end type
  type routing_type
     integer, allocatable :: gindex(:,:)
  end type
  type(rc_type) :: LIS_rc
  type(routing_type), allocatable :: LIS_routing(:)
  integer :: LIS_localpet=0
  integer :: LIS_ews_halo_ind(1,1)=1,LIS_ewe_halo_ind(1,1)=1
  integer :: LIS_nss_halo_ind(1,1)=1,LIS_nse_halo_ind(1,1)=1
end module

module LIS_routingMod
  implicit none
end module

module LIS_logMod
  implicit none
contains
  subroutine LIS_verify(status,message)
    integer, intent(in) :: status
    character(*), intent(in) :: message
    if(status/=0) then
       print *,message,status
       stop 1
    endif
  end subroutine
end module

module HYMAP3_routingMod
  implicit none
  type hymap_type
     integer :: nseqall
     integer, allocatable :: seqx(:),seqy(:)
  end type
  type(hymap_type),allocatable :: HYMAP3_routing_struc(:)
end module

module HYMAP3_daSWOT_Mod
  implicit none
  type da_type
     integer :: useLocalUpd,localupdDX
     real, allocatable :: sites(:,:),localWeight(:,:,:)
  end type
  type(da_type),allocatable :: HYMAP3_daSWOT_struc(:)
end module

! INSERT_ACTUAL_ROUTINES_HERE

program regression
  use ESMF
  use LIS_coreMod
  use HYMAP3_routingMod
  use HYMAP3_daSWOT_Mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  real, allocatable, target :: state(:),incr(:),refstate(:),refincr(:)
  type(ESMF_State) :: actual_state,actual_incr,reference_state,reference_incr
  integer :: i,j,m,t,c,r,fixture
  real :: expected
  character(32) :: mode
  integer :: checks=0
  call get_command_argument(1,mode)

  if(trim(mode)=='--legacy-dense') then
     call setup(21,21,3,10,1)
     incr=0.1
     HYMAP3_daSWOT_struc(1)%localWeight=1.0
     ! Expected bounds-check failure: the 101st contributor exceeds DX**2.
     call HYMAP3_updateSWOT_legacy(1,actual_state,actual_incr)
     print *,'ERROR: legacy dense call unexpectedly completed'
     stop 2
  endif

  do fixture=1,4
     call setup(4,3,3,4,1)
     do i=1,HYMAP3_routing_struc(1)%nseqall
        do m=1,3
           t=(i-1)*3+m
           incr(t)=real(mod(7*i+3*m,21)-10)
        enddo
        do r=1,3
           do c=1,4
              HYMAP3_daSWOT_struc(1)%localWeight(c,r,i)=real(mod(5*i+c+2*r,17))/16.0
           enddo
        enddo
     enddo
     select case(fixture)
     case(1)
        ! Rejection/capping boundaries and member indexing.
        incr(1:9)=[-10.0,-9.999,-5.1,-5.0,0.0,5.0,5.1,9.999,10.0]
        state(10)=-9999.0
     case(2)
        ! Equal weights: retain first contributor including source-cell updates.
        HYMAP3_daSWOT_struc(1)%localWeight=0.5
     case(3)
        ! Preserve legacy missing-site direct writes and missing self weights.
        HYMAP3_daSWOT_struc(1)%sites(1,2)=-9999.0
        HYMAP3_daSWOT_struc(1)%sites(3,2)=-9999.0
        HYMAP3_daSWOT_struc(1)%localWeight(2,1,4)=-9999.0
        HYMAP3_daSWOT_struc(1)%localWeight(4,1,10)=-9999.0
     case(4)
        ! Existing finite zero/negative/>1 weight behavior remains unchanged.
        do i=1,12
           do r=1,3
              do c=1,4
                 HYMAP3_daSWOT_struc(1)%localWeight(c,r,i)=real(mod(i+c+3*r,29)-10)/10.0
              enddo
           enddo
        enddo
     end select
     call compare_legacy()
  enddo

  call setup(4,3,3,4,1)
  incr=0.0
  call compare_legacy()

  call setup(21,21,3,10,1)
  HYMAP3_daSWOT_struc(1)%localWeight=1.0
  do i=1,441
     do m=1,3
        incr((i-1)*3+m)=real(m)*0.1
     enddo
  enddo
  call HYMAP3_updateSWOT(1,actual_state,actual_incr)
  do i=1,441
     do m=1,3
        expected=25.0+real(m)*0.1
        if(state((i-1)*3+m)/=expected) stop 3
     enddo
  enddo
  checks=checks+1

  ! Disabled LZN: no sites/map allocation, poisoned radius must not be used.
  call setup(4,3,3,-huge(0),0)
  incr=9.0
  call HYMAP3_updateSWOT(1,actual_state,actual_incr)
  if(any(state/=30.0).or.any(incr/=5.0)) stop 4
  checks=checks+1

  ! Radius zero has a well-defined self-only update without any history extent.
  call setup(4,3,3,0,1)
  incr=1.0
  HYMAP3_daSWOT_struc(1)%localWeight=1.0
  call HYMAP3_updateSWOT(1,actual_state,actual_incr)
  if(any(state/=26.0)) stop 5
  checks=checks+1

  ! Nonfinite raw increments have the existing behavior: NaN ignored, inf rejected.
  call setup(1,1,3,1,1)
  incr=[ieee_value(0.0,ieee_quiet_nan),ieee_value(0.0,ieee_positive_inf),0.0]
  HYMAP3_daSWOT_struc(1)%localWeight=1.0
  call HYMAP3_updateSWOT(1,actual_state,actual_incr)
  if(any(state/=25.0)) stop 6
  checks=checks+1
  print '(a,i0)','PASS synthetic Fortran fixtures: ',checks
contains
  subroutine setup(nx,ny,nens,dx,localized)
    integer,intent(in) :: nx,ny,nens,dx,localized
    integer :: x,y,index,sz
    if(allocated(LIS_rc%lnc)) deallocate(LIS_rc%lnc,LIS_rc%lnr,LIS_rc%nensem)
    if(allocated(LIS_routing)) deallocate(LIS_routing)
    if(allocated(HYMAP3_routing_struc)) deallocate(HYMAP3_routing_struc)
    if(allocated(HYMAP3_daSWOT_struc)) deallocate(HYMAP3_daSWOT_struc)
    if(allocated(state)) deallocate(state,incr,refstate,refincr)
    allocate(LIS_rc%lnc(1),LIS_rc%lnr(1),LIS_rc%nensem(1),LIS_routing(1))
    allocate(HYMAP3_routing_struc(1),HYMAP3_daSWOT_struc(1))
    LIS_rc%lnc=nx; LIS_rc%lnr=ny; LIS_rc%nensem=nens
    LIS_ewe_halo_ind=nx; LIS_nse_halo_ind=ny
    HYMAP3_routing_struc(1)%nseqall=nx*ny
    allocate(LIS_routing(1)%gindex(nx,ny))
    allocate(HYMAP3_routing_struc(1)%seqx(nx*ny),HYMAP3_routing_struc(1)%seqy(nx*ny))
    HYMAP3_daSWOT_struc(1)%useLocalUpd=localized
    HYMAP3_daSWOT_struc(1)%localupdDX=dx
    if(localized==1) then
       allocate(HYMAP3_daSWOT_struc(1)%sites(nx,ny))
       allocate(HYMAP3_daSWOT_struc(1)%localWeight(nx,ny,nx*ny))
    endif
    index=0
    do x=1,nx
       do y=1,ny
          index=index+1
          HYMAP3_routing_struc(1)%seqx(index)=x
          HYMAP3_routing_struc(1)%seqy(index)=y
          LIS_routing(1)%gindex(x,y)=index
          if(localized==1) HYMAP3_daSWOT_struc(1)%sites(x,y)=real(index)
       enddo
    enddo
    sz=nx*ny*nens
    allocate(state(sz),incr(sz),refstate(sz),refincr(sz))
    state=25.0; incr=0.0
    actual_state%values=>state; actual_incr%values=>incr
    reference_state%values=>refstate; reference_incr%values=>refincr
  end subroutine

  subroutine compare_legacy()
    refstate=state; refincr=incr
    call HYMAP3_updateSWOT_legacy(1,reference_state,reference_incr)
    call HYMAP3_updateSWOT(1,actual_state,actual_incr)
    if(any(state/=refstate).or.any(incr/=refincr)) then
       print *,'Mismatch fixture ',fixture,maxval(abs(state-refstate))
       stop 7
    endif
    checks=checks+1
  end subroutine
end program
