!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

module RAPID_damStateMod
!BOP
!
! !MODULE: RAPID_damStateMod
!
! !DESCRIPTION:
!  Preserve reservoir and routing state across LIS coupling calls, ensemble
!  trajectories, and NetCDF restarts. Estimate storage from initial discharge
!  when no full reservoir restart state is available.
!
!  Source: Victor Roland supplied the Qinit-to-storage initialization law in
!  standalone RAPID version B, src/rapid_init.F90. Its original patch date was
!  not provided. This attribution covers the inverse law, not the LIS state
!  management, NetCDF restart schema, or validation routines below.
!
! !REVISION HISTORY:
! 03 Nov 2026: Yeosang Yoon: Initial LIS dam state management and restart I/O;
!                          reuse Victor Roland's supplied initialization law.
! 03 Nov 2026: Yeosang Yoon: Match LIS MPI include/module bindings for broadcast.
!
!EOP
  use ESMF, only: ESMF_KIND_R8
  use RAPID_routingMod, only: RAPID_routing_struc
  use LIS_coreMod, only: LIS_rc
  use LIS_logMod, only: LIS_logunit, LIS_verify, LIS_endrun
  use netcdf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef PETSc
#include <petsc/finclude/petscmat.h>
  use petscmat
#ifdef SPMD
#ifdef USE_INCLUDE_MPI
  use LIS_mpiMod, only: LIS_MPI_INTEGER => MPI_INTEGER
#else
  use LIS_mpiMod, only: MPI_Bcast, LIS_MPI_INTEGER => MPI_INTEGER
#endif
#endif
  use rapid_var, only: rank, ierr, IS_one, ZS_one, IS_riv_bas, IS_dam_tot, &
       IV_riv_loc1, IV_riv_bas_id, IV_dam_tot_id, IV_dam_pos, &
       ZV_S_dam, ZV_Smin_dam, ZV_Smax_dam, ZV_k_dam, ZV_p_dam, &
       ZV_Qin_dam_prev, ZV_Qout_dam_prev, ZV_Qout_dam0, &
       ZV_QoutinitR, ZV_QoutR, ZV_QoutbarR, ZV_Qext, ZV_QinbarR, ZM_Net, &
       vecscat, ZV_SeqZero, ZV_pointer, ZS_TauR, ZS_dtR, &
       BS_opt_Qinit, BS_lis_dam_inflow_ready
#endif
  implicit none
  private
  public :: RAPID_dam_read_restart, RAPID_dam_define_restart, RAPID_dam_write_restart
#ifdef PETSc
  public :: RAPID_dam_warm_start, RAPID_dam_restore, RAPID_dam_capture
#endif
  integer, parameter :: dam_state_version = 1
  character(len=32), parameter :: dam_names(3) = [character(len=32) :: &
       'rapid_dam_storage', 'rapid_dam_inflow_next', 'rapid_dam_outflow_prev']

contains

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(LIS_logunit,*) '[ERR] RAPID dam state: ',trim(message)
    call LIS_endrun()
  end subroutine fail

  integer function trajectories(n)
    integer, intent(in) :: n
    trajectories = 1
    if (RAPID_routing_struc(n)%useens == 2) trajectories = LIS_rc%nensem(n)
  end function trajectories

  subroutine allocate_state(n, ndam)
    integer, intent(in) :: n, ndam
    integer :: nriv, nt
    nriv = RAPID_routing_struc(n)%n_riv_bas
    nt = trajectories(n)
    associate(r => RAPID_routing_struc(n))
      if (allocated(r%dam_storage)) call fail('state allocated more than once')
      allocate(r%dam_ids(ndam), r%dam_river_ids(nriv), r%dam_parameters(ndam,4))
      allocate(r%dam_storage(ndam,nt), r%dam_inflow_next(ndam,nt), r%dam_outflow_prev(ndam,nt))
      allocate(r%dam_river_qout(nriv,nt), r%dam_state_ready(nt))
      r%dam_state_ready = .false.
    end associate
  end subroutine allocate_state

  subroutine read_real(ncid, name, values, dims)
    integer, intent(in) :: ncid, dims(:)
    character(len=*), intent(in) :: name
    real(ESMF_KIND_R8), intent(out) :: values(:,:)
    integer :: varid, ndims, actual_dims(NF90_MAX_VAR_DIMS)
    call LIS_verify(nf90_inq_varid(ncid,name,varid),'Missing RAPID dam restart variable '//name)
    call LIS_verify(nf90_inquire_variable(ncid,varid,ndims=ndims,dimids=actual_dims), &
         'Invalid RAPID dam restart variable '//name)
    if (ndims /= size(dims)) call fail('wrong rank for '//name)
    if (any(actual_dims(1:ndims) /= dims)) call fail('wrong dimensions for '//name)
    call LIS_verify(nf90_get_var(ncid,varid,values),'Cannot read '//name)
    if (.not. all(ieee_is_finite(values))) call fail('nonfinite values in '//name)
  end subroutine read_real

  subroutine read_ids(ncid, name, values, dimid)
    integer, intent(in) :: ncid, dimid
    character(len=*), intent(in) :: name
    integer, intent(out) :: values(:)
    integer :: varid, ndims, actual_dims(NF90_MAX_VAR_DIMS)
    call LIS_verify(nf90_inq_varid(ncid,name,varid),'Missing '//name)
    call LIS_verify(nf90_inquire_variable(ncid,varid,ndims=ndims,dimids=actual_dims),'Invalid '//name)
    if (ndims /= 1) call fail('wrong rank for '//name)
    if (actual_dims(1) /= dimid) call fail('wrong dimension for '//name)
    call LIS_verify(nf90_get_var(ncid,varid,values),'Cannot read '//name)
  end subroutine read_ids

  ! Called on the LIS master before the first PETSc model invocation.
  subroutine RAPID_dam_read_restart(n, ncid)
    integer, intent(in) :: n, ncid
    integer :: status, version, dd, dm, dr, dp, ndam, nt, nriv, np, varid
    if (.not. RAPID_routing_struc(n)%bdam) return
    status = nf90_get_att(ncid,NF90_GLOBAL,'rapid_dam_state_version',version)
    if (status == NF90_ENOTATT) then
      if (nf90_inq_varid(ncid,'rapid_dam_storage',varid) == NF90_NOERR) &
           call fail('dam variables are present but schema version is missing')
      write(LIS_logunit,*) '[INFO] RAPID dam: legacy Qout-only restart; storage will be estimated.'
      return
    endif
    call LIS_verify(status,'Cannot read RAPID dam restart version')
    if (version /= dam_state_version) call fail('unsupported restart schema')
    call LIS_verify(nf90_inq_dimid(ncid,'rapid_dam',dd),'Missing rapid_dam dimension')
    call LIS_verify(nf90_inq_dimid(ncid,'rapid_dam_member',dm),'Missing rapid_dam_member dimension')
    call LIS_verify(nf90_inq_dimid(ncid,'rivid',dr),'Missing rivid dimension')
    call LIS_verify(nf90_inq_dimid(ncid,'rapid_dam_parameter',dp),'Missing rapid_dam_parameter dimension')
    call LIS_verify(nf90_inquire_dimension(ncid,dd,len=ndam),'Invalid rapid_dam dimension')
    call LIS_verify(nf90_inquire_dimension(ncid,dm,len=nt),'Invalid rapid_dam_member dimension')
    call LIS_verify(nf90_inquire_dimension(ncid,dr,len=nriv),'Invalid rivid dimension')
    call LIS_verify(nf90_inquire_dimension(ncid,dp,len=np),'Invalid rapid_dam_parameter dimension')
    if (ndam <= 0 .or. np /= 4) call fail('invalid dam parameter dimensions')
    if (nt /= trajectories(n)) call fail('restart trajectory count differs from configuration')
    if (nriv /= RAPID_routing_struc(n)%n_riv_bas) call fail('restart river count differs from configuration')
    call allocate_state(n,ndam)
    associate(r => RAPID_routing_struc(n))
      call read_ids(ncid,'rapid_dam_id',r%dam_ids,dd)
      call read_ids(ncid,'rapid_dam_river_id',r%dam_river_ids,dr)
      call read_real(ncid,'rapid_dam_parameters',r%dam_parameters,[dd,dp])
      call read_real(ncid,dam_names(1),r%dam_storage,[dd,dm])
      call read_real(ncid,dam_names(2),r%dam_inflow_next,[dd,dm])
      call read_real(ncid,dam_names(3),r%dam_outflow_prev,[dd,dm])
      call read_real(ncid,'rapid_dam_river_qout',r%dam_river_qout,[dr,dm])
      call LIS_verify(nf90_get_att(ncid,NF90_GLOBAL,'rapid_dam_tau',r%dam_tau),'Missing rapid_dam_tau')
      call LIS_verify(nf90_get_att(ncid,NF90_GLOBAL,'rapid_dam_dt',r%dam_dt),'Missing rapid_dam_dt')
      if (.not. ieee_is_finite(r%dam_tau) .or. .not. ieee_is_finite(r%dam_dt)) &
           call fail('invalid restart time steps')
      r%dam_state_ready = .true.
      r%dam_restart_loaded = .true.
    end associate
  end subroutine RAPID_dam_read_restart

  subroutine RAPID_dam_define_restart(n, ncid)
    integer, intent(in) :: n, ncid
    integer :: dd, dm, dr, dp, varid, j
    if (.not. RAPID_routing_struc(n)%bdam) return
    associate(r => RAPID_routing_struc(n))
      if (.not. allocated(r%dam_storage)) return
      if (.not. all(r%dam_state_ready)) call fail('not all dam trajectories have been initialized')
      call LIS_verify(nf90_put_att(ncid,NF90_GLOBAL,'rapid_dam_state_version',dam_state_version), &
           'Cannot define RAPID dam schema')
      call LIS_verify(nf90_put_att(ncid,NF90_GLOBAL,'rapid_dam_tau',r%dam_tau),'Cannot define dam tau')
      call LIS_verify(nf90_put_att(ncid,NF90_GLOBAL,'rapid_dam_dt',r%dam_dt),'Cannot define dam dt')
      call LIS_verify(nf90_def_dim(ncid,'rapid_dam',size(r%dam_ids),dd),'Cannot define dam dimension')
      call LIS_verify(nf90_def_dim(ncid,'rapid_dam_member',size(r%dam_state_ready),dm), &
           'Cannot define dam trajectory dimension')
      call LIS_verify(nf90_def_dim(ncid,'rapid_dam_parameter',4,dp),'Cannot define dam parameter dimension')
      call LIS_verify(nf90_inq_dimid(ncid,'rivid',dr),'Missing rivid dimension')
      call LIS_verify(nf90_def_var(ncid,'rapid_dam_id',NF90_INT,[dd],varid),'Cannot define dam IDs')
      call LIS_verify(nf90_def_var(ncid,'rapid_dam_river_id',NF90_INT,[dr],varid),'Cannot define river IDs')
      call LIS_verify(nf90_def_var(ncid,'rapid_dam_parameters',NF90_DOUBLE,[dd,dp],varid), &
           'Cannot define dam parameters')
      call LIS_verify(nf90_put_att(ncid,varid,'parameter_order','Smax Smin p k'),'Cannot define parameter order')
      do j=1,size(dam_names)
        call LIS_verify(nf90_def_var(ncid,trim(dam_names(j)),NF90_DOUBLE,[dd,dm],varid), &
             'Cannot define '//trim(dam_names(j)))
      enddo
      call LIS_verify(nf90_def_var(ncid,'rapid_dam_river_qout',NF90_DOUBLE,[dr,dm],varid), &
           'Cannot define full-precision routing state')
    end associate
  end subroutine RAPID_dam_define_restart

  subroutine put_real(ncid, name, values)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(ESMF_KIND_R8), intent(in) :: values(:,:)
    integer :: varid
    call LIS_verify(nf90_inq_varid(ncid,name,varid),'Missing '//name)
    call LIS_verify(nf90_put_var(ncid,varid,values),'Cannot write '//name)
  end subroutine put_real

  subroutine RAPID_dam_write_restart(n, ncid)
    integer, intent(in) :: n, ncid
    integer :: varid
    if (.not. RAPID_routing_struc(n)%bdam) return
    associate(r => RAPID_routing_struc(n))
      if (.not. allocated(r%dam_storage)) return
      call LIS_verify(nf90_inq_varid(ncid,'rapid_dam_id',varid),'Missing dam IDs')
      call LIS_verify(nf90_put_var(ncid,varid,r%dam_ids),'Cannot write dam IDs')
      call LIS_verify(nf90_inq_varid(ncid,'rapid_dam_river_id',varid),'Missing river IDs')
      call LIS_verify(nf90_put_var(ncid,varid,r%dam_river_ids),'Cannot write river IDs')
      call put_real(ncid,'rapid_dam_parameters',r%dam_parameters)
      call put_real(ncid,dam_names(1),r%dam_storage)
      call put_real(ncid,dam_names(2),r%dam_inflow_next)
      call put_real(ncid,dam_names(3),r%dam_outflow_prev)
      call put_real(ncid,'rapid_dam_river_qout',r%dam_river_qout)
    end associate
  end subroutine RAPID_dam_write_restart

#ifdef PETSc
  ! Victor Roland's supplied inverse release law. Its use remains an estimate
  ! when the supplied discharge is not an actual, time-aligned dam release.
  real(ESMF_KIND_R8) function storage_from_flow(q, j) result(storage)
    real(ESMF_KIND_R8), intent(in) :: q
    integer, intent(in) :: j
    real(ESMF_KIND_R8) :: capacity
    if (.not. ieee_is_finite(q)) call fail('nonfinite initial discharge')
    capacity = ZV_Smax_dam(j)-ZV_Smin_dam(j)
    if (.not. ieee_is_finite(capacity) .or. capacity <= 0) call fail('Smax must exceed Smin')
    if (.not. ieee_is_finite(ZV_k_dam(j)) .or. ZV_k_dam(j) <= 0) call fail('k must be positive')
    if (.not. ieee_is_finite(ZV_p_dam(j)) .or. ZV_p_dam(j) < 0) call fail('p must be nonnegative')
    if (.not. ieee_is_finite(ZS_TauR) .or. ZS_TauR <= 0) call fail('TauR must be positive')
    storage = ZV_Smin_dam(j)
    if (q > 0) then
      storage = ZV_Smin_dam(j) + &
           (q*ZS_TauR*capacity**ZV_p_dam(j)/ZV_k_dam(j))**(1.0_8/(1.0_8+ZV_p_dam(j)))
      if (.not. ieee_is_finite(storage)) call fail('initial storage inversion overflow')
      storage = min(max(storage,real(ZV_Smin_dam(j),ESMF_KIND_R8)), &
           real(ZV_Smax_dam(j),ESMF_KIND_R8))
    endif
  end function storage_from_flow

  subroutine gather_vector(vector, values)
    Vec, intent(in) :: vector
    real(ESMF_KIND_R8), allocatable, intent(out) :: values(:)
    call VecScatterBegin(vecscat,vector,ZV_SeqZero,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(vecscat,vector,ZV_SeqZero,INSERT_VALUES,SCATTER_FORWARD,ierr)
    if (rank == 0) then
      call VecGetArrayF90(ZV_SeqZero,ZV_pointer,ierr)
      allocate(values(IS_riv_bas))
      values = ZV_pointer
      call VecRestoreArrayF90(ZV_SeqZero,ZV_pointer,ierr)
    endif
  end subroutine gather_vector

  subroutine RAPID_dam_warm_start(vector)
    Vec, intent(in) :: vector
    real(ESMF_KIND_R8), allocatable :: flow(:)
    integer :: j
    call gather_vector(vector,flow)
    if (rank == 0) then
      ZV_S_dam = ZV_Smin_dam
      do j=1,IS_dam_tot
        if (IV_dam_pos(j) > 0) then
          ZV_Qout_dam0(j) = flow(IV_dam_pos(j))
          ZV_S_dam(j) = storage_from_flow(flow(IV_dam_pos(j)),j)
        endif
      enddo
    endif
  end subroutine RAPID_dam_warm_start

  subroutine check_geometry(n)
    integer, intent(in) :: n
    integer :: j
    associate(r => RAPID_routing_struc(n))
      do j=1,IS_dam_tot
        if (.not. ieee_is_finite(storage_from_flow(0.0_8,j))) call fail('invalid initial storage')
      enddo
      if (size(r%dam_ids) /= IS_dam_tot) call fail('restart dam count changed')
      if (any(r%dam_ids /= IV_dam_tot_id)) call fail('restart dam ID order changed')
      if (any(r%dam_river_ids /= IV_riv_bas_id)) call fail('restart basin ID order changed')
      if (r%dam_tau /= ZS_TauR .or. r%dam_dt /= ZS_dtR) call fail('restart routing time steps changed')
      if (any(r%dam_parameters(:,1) /= ZV_Smax_dam) .or. &
          any(r%dam_parameters(:,2) /= ZV_Smin_dam) .or. &
          any(r%dam_parameters(:,3) /= ZV_p_dam) .or. &
          any(r%dam_parameters(:,4) /= ZV_k_dam)) call fail('restart dam parameters changed')
      if (any(r%dam_storage < spread(real(ZV_Smin_dam,ESMF_KIND_R8),2,size(r%dam_storage,2)))) &
           call fail('restart storage is below Smin')
      ! Do not clip restored storage at Smax: the existing forward model can
      ! temporarily exceed Smax, and a restart must retain its actual state.
      if (r%dam_restart_loaded) write(LIS_logunit,*) '[INFO] RAPID dam: restoring full restart state.'
      r%dam_state_checked = .true.
    end associate
  end subroutine check_geometry

  subroutine next_inflow(vector, include_external, values)
    Vec, intent(in) :: vector
    logical, intent(in) :: include_external
    real(ESMF_KIND_R8), allocatable, intent(out) :: values(:)
    real(ESMF_KIND_R8), allocatable :: river_inflow(:)
    integer :: j
    call MatMult(ZM_Net,vector,ZV_QinbarR,ierr)
    if (include_external) call VecAXPY(ZV_QinbarR,ZS_one,ZV_Qext,ierr)
    call gather_vector(ZV_QinbarR,river_inflow)
    if (rank == 0) then
      allocate(values(IS_dam_tot))
      values = 0
      do j=1,IS_dam_tot
        if (IV_dam_pos(j) > 0) values(j) = river_inflow(IV_dam_pos(j))
      enddo
    endif
  end subroutine next_inflow

  subroutine RAPID_dam_restore(n, member)
    integer, intent(in) :: n, member
    integer :: trajectory, initialize, j
#if defined(SPMD) && defined(USE_INCLUDE_MPI)
    ! LIS uses mpif.h in this configuration; MPI_Bcast is an external routine.
    external :: MPI_Bcast
#endif
    PetscScalar, allocatable :: routed_flow(:)
    real(ESMF_KIND_R8), allocatable :: flow(:), inflow(:)
    trajectory = max(1,member)
    initialize = 0
    if (rank == 0) then
      associate(r => RAPID_routing_struc(n))
        if (.not. allocated(r%dam_storage)) then
          call allocate_state(n,IS_dam_tot)
          r%dam_ids = IV_dam_tot_id
          r%dam_river_ids = IV_riv_bas_id
          r%dam_parameters(:,1) = ZV_Smax_dam
          r%dam_parameters(:,2) = ZV_Smin_dam
          r%dam_parameters(:,3) = ZV_p_dam
          r%dam_parameters(:,4) = ZV_k_dam
          r%dam_tau = ZS_TauR
          r%dam_dt = ZS_dtR
          r%dam_storage = spread(real(ZV_Smin_dam,ESMF_KIND_R8),2,trajectories(n))
        endif
        if (.not. r%dam_state_checked) call check_geometry(n)
        if (.not. r%dam_state_ready(trajectory)) initialize = 1
      end associate
    endif
#ifdef SPMD
    call MPI_Bcast(initialize,1,LIS_MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
#endif
    if (initialize == 1) then
      call gather_vector(ZV_QoutinitR,flow)
      ! Legacy restarts have no previous mean flow/external inflow. Use their
      ! initial river discharge as a mean-flow estimate and zero external flow.
      call next_inflow(ZV_QoutinitR,.false.,inflow)
      if (rank == 0) then
        associate(r => RAPID_routing_struc(n))
          r%dam_river_qout(:,trajectory) = flow
          r%dam_inflow_next(:,trajectory) = inflow
          r%dam_outflow_prev(:,trajectory) = 0
          if (r%startmode == 'restart' .or. BS_opt_Qinit) then
            do j=1,IS_dam_tot
              if (IV_dam_pos(j) > 0) then
                r%dam_storage(j,trajectory) = storage_from_flow(flow(IV_dam_pos(j)),j)
                r%dam_outflow_prev(j,trajectory) = flow(IV_dam_pos(j))
              endif
            enddo
            write(LIS_logunit,*) '[INFO] RAPID dam: estimated initial storage for trajectory ',trajectory
          endif
          r%dam_state_ready(trajectory) = .true.
        end associate
      endif
    endif
    if (rank == 0) then
      associate(r => RAPID_routing_struc(n))
        ZV_S_dam = r%dam_storage(:,trajectory)
        ZV_Qin_dam_prev = r%dam_inflow_next(:,trajectory)
        ZV_Qout_dam_prev = r%dam_outflow_prev(:,trajectory)
        allocate(routed_flow(IS_riv_bas))
        routed_flow = r%dam_river_qout(:,trajectory)
        call VecSetValues(ZV_QoutinitR,IS_riv_bas,IV_riv_loc1,routed_flow,INSERT_VALUES,ierr)
      end associate
    endif
    call VecAssemblyBegin(ZV_QoutinitR,ierr)
    call VecAssemblyEnd(ZV_QoutinitR,ierr)
    BS_lis_dam_inflow_ready = .true.
  end subroutine RAPID_dam_restore

  subroutine RAPID_dam_capture(n, member)
    integer, intent(in) :: n, member
    integer :: trajectory
    real(ESMF_KIND_R8), allocatable :: flow(:), inflow(:)
    trajectory = max(1,member)
    call gather_vector(ZV_QoutR,flow)
    ! This is exactly Net*previous_QoutbarR + previous_Qext used by the next
    ! standalone dam call. Store only its dam entries, not two full river vectors.
    call next_inflow(ZV_QoutbarR,.true.,inflow)
    if (rank == 0) then
      associate(r => RAPID_routing_struc(n))
        if (.not. all(ieee_is_finite(ZV_S_dam))) call fail('nonfinite simulated storage')
        if (.not. all(ieee_is_finite(flow)) .or. .not. all(ieee_is_finite(inflow))) &
             call fail('nonfinite simulated flow')
        r%dam_storage(:,trajectory) = ZV_S_dam
        r%dam_inflow_next(:,trajectory) = inflow
        r%dam_outflow_prev(:,trajectory) = ZV_Qout_dam_prev
        r%dam_river_qout(:,trajectory) = flow
      end associate
    endif
  end subroutine RAPID_dam_capture
#endif
end module RAPID_damStateMod
