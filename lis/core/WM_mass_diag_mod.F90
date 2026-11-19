! Optional routing mass diagnostics. Runtime selection stays outside source.
! Test build: HYMAP3_MASS_DIAG=1 ./compile, using the existing configuration.
! Logging is off unless WM_DIAG_DIR is set to a fresh existing directory.
! WM_DIAG_CELLS: one positive row-major global cell ID per line (one-based).
! WM_DIAG_OUTLET_CELL must be in that list; WM_DIAG_EXPECTED_CELLS is optional.
! WM_DIAG_START and WM_DIAG_END: inclusive YYYYMMDDHHMM model clock values.
! Output: tick, context, phase, member, owned-cell count, phase-specific values.
! Each rank writes its own file; sum ranks separately for each member/phase.
! No MPI reduction or random-number call is introduced.
module WM_mass_diag_mod
 use iso_fortran_env,only:int64,real64
 implicit none
 character(16)::wm_context='NONE'
 integer::wm_member=0,wm_outlet_cell=0
 logical,save::initialized=.false.,enabled=.false.
 logical,allocatable,save::selected(:)
 integer(int64),save::start_tick=0,end_tick=0
 character(1024),save::outdir=''
 integer,save::unit=-1
contains
 logical function wm_active(tick)
  integer(int64),intent(in)::tick
  character(1024)::cellfile,buf
  integer::stat,u,cell,n,maxcell,expected
  if(.not.initialized)then
   initialized=.true.
   call get_environment_variable('WM_DIAG_DIR',outdir,status=stat)
   enabled=stat==0.and.len_trim(outdir)>0
   if(enabled)then
    call get_environment_variable('WM_DIAG_CELLS',cellfile,status=stat)
    if(stat/=0)error stop 'WM_DIAG_CELLS required'
    open(newunit=u,file=trim(cellfile),status='old',action='read')
    n=0;maxcell=0
    do
     read(u,*,iostat=stat)cell
     if(stat<0)exit
     if(stat>0)error stop 'WM invalid cell list'
     if(cell<1)error stop 'WM invalid cell list'
     n=n+1;maxcell=max(maxcell,cell)
    enddo
    if(n<1)error stop 'WM empty cell list'
    allocate(selected(maxcell));selected=.false.
    rewind(u)
    do
     read(u,*,iostat=stat)cell
     if(stat<0)exit
     if(stat>0)error stop 'WM changed cell list'
     if(cell<1.or.cell>maxcell)error stop 'WM changed cell list'
     if(selected(cell))error stop 'WM duplicate cell'
     selected(cell)=.true.
    enddo
    close(u)
    call get_environment_variable('WM_DIAG_EXPECTED_CELLS',buf,status=stat)
    if(stat==0)then
     read(buf,*,iostat=stat)expected
     if(stat/=0)error stop 'WM invalid expected cell count'
     if(expected/=count(selected))error stop 'WM selected cell count mismatch'
    endif
    call get_environment_variable('WM_DIAG_OUTLET_CELL',buf,status=stat)
    if(stat/=0)error stop 'WM_DIAG_OUTLET_CELL required'
    read(buf,*,iostat=stat)wm_outlet_cell
    if(stat/=0)error stop 'WM invalid outlet cell'
    if(wm_outlet_cell<1.or.wm_outlet_cell>maxcell)error stop 'WM outlet outside selection'
    if(.not.selected(wm_outlet_cell))error stop 'WM outlet not selected'
    call get_environment_variable('WM_DIAG_START',buf,status=stat)
    if(stat/=0)error stop 'WM_DIAG_START required'
    read(buf,*,iostat=stat)start_tick
    if(stat/=0)error stop 'WM invalid start tick'
    call get_environment_variable('WM_DIAG_END',buf,status=stat)
    if(stat/=0)error stop 'WM_DIAG_END required'
    read(buf,*,iostat=stat)end_tick
    if(stat/=0)error stop 'WM invalid end tick'
    if(start_tick>end_tick)error stop 'WM reversed diagnostic interval'
   endif
  endif
  wm_active=enabled.and.tick>=start_tick.and.tick<=end_tick
 end function
 logical function wm_contains(cell)
  integer,intent(in)::cell
  wm_contains=.false.
  if(.not.allocated(selected))return
  if(cell<1.or.cell>size(selected))return
  wm_contains=selected(cell)
 end function
 subroutine wm_emit(rank,tick,phase,member,ncell,values)
  integer,intent(in)::rank,member,ncell
  integer(int64),intent(in)::tick
  character(*),intent(in)::phase
  real(real64),intent(in)::values(:)
  character(1024)::file
  character(12)::rankstr
  if(.not.wm_active(tick))return
  if(unit==-1)then
   write(rankstr,'(I6.6)')rank
   file=trim(outdir)//'/rank_'//trim(rankstr)//'.txt'
   open(newunit=unit,file=trim(file),status='new',action='write')
  endif
  write(unit,'(I0,1X,A,1X,A,1X,I0,1X,I0,1X,*(ES25.16,1X))')tick,trim(wm_context),trim(phase),member,ncell,values
 end subroutine
end module
