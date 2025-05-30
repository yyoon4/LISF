!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: USAFSIobs_Mod
!
! REVISION HISTORY:
! 09 Apr 2019  Eric Kemp  Initial specification
!
! DESCRIPTION:
! Source code for using USAF Snow and Ice analysis.
!------------------------------------------------------------------------------

#include "LIS_misc.h"

module USAFSIobs_Mod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: USAFSIobs_finalize
   public :: USAFSIobs_setup

   ! Public variables
   public :: USAFSI_obs, USAFSI_obs_t

   ! Declare type
   type :: USAFSI_obs_t
      real, allocatable :: rlat1(:)
      real, allocatable :: rlon1(:)
      integer, allocatable :: n111(:)
      integer :: mo1 ! On LIS grid
      integer :: nc_ldt ! On LDT grid
      integer :: nr_ldt ! On LDT grid
      integer :: nc_lis ! On LIS grid
      integer :: nr_lis ! On LIS grid
      integer :: mi ! On LDT grid
      real :: gridDesco(50)
   end type USAFSI_obs_t

   type(USAFSI_obs_t), allocatable :: USAFSI_obs(:)

contains

   ! Complete rnntime initializations for USAFSI product
   subroutine USAFSIobs_setup(k, OBS_State, OBS_Pert_State)

      ! Imports
      use ESMF
      use LIS_constantsMod, only: LIS_CONST_PATH_LEN
      use LIS_coreMod, only: LIS_rc, LIS_config
      use LIS_DAobservationsMod, only: LIS_obsVecGrid, LIS_obsEnsOnGrid
      use LIS_logMod, only: LIS_verify, LIS_logunit, LIS_getNextUnitNumber, &
           LIS_releaseUnitNumber, LIS_endrun
      use LIS_perturbMod, only: LIS_readPertAttributes, pert_dec_type
      use LIS_timeMgrMod, only: LIS_registerAlarm

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: k
      type(ESMF_State), intent(inout) :: OBS_State(LIS_rc%nnest)
      type(ESMF_State), intent(inout) :: OBS_Pert_State(LIS_rc%nnest)

      ! Local variables
      integer                ::  ftn
      integer                ::  i
      type(ESMF_ArraySpec)   ::  intarrspec
      character(len=LIS_CONST_PATH_LEN) ::  USAFSIobsdir
      character(40)          ::  USAFSI_infile_name
      integer                ::  n
      type(pert_dec_type)    ::  obs_pert
      real, pointer          ::  obs_temp(:,:)
      type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
      type(ESMF_ArraySpec)   ::  pertArrSpec
      type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
      type(ESMF_ArraySpec)   ::  realarrspec
      real,  allocatable     ::  ssdev(:)
      integer                ::  status
      character(100)          ::  temp
      real,   allocatable    ::  varmin(:)
      real,   allocatable    ::  varmax(:)
      character(1)            ::  vid(2)
      character(40), allocatable  ::  vname(:)
      
      allocate(USAFSI_obs(LIS_rc%nnest))
      call ESMF_ArraySpecSet(intarrspec, rank=1, typekind=ESMF_TYPEKIND_I4,&
           rc=status)
      call LIS_verify(status)

      call ESMF_ArraySpecSet(realarrspec, rank=1, typekind=ESMF_TYPEKIND_R4,&
           rc=status)
      call LIS_verify(status)

      call ESMF_ArraySpecSet(pertArrSpec, rank=2, typekind=ESMF_TYPEKIND_R4,&
           rc=status)
      call LIS_verify(status)

      call ESMF_ConfigFindLabel(LIS_config, "USAFSI data directory:",&
           rc=status)
      do n = 1, LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config, USAFSIobsdir,&
              rc=status)
         if (status .ne. ESMF_SUCCESS)then
            write(LIS_logunit,*) "[ERR] USAFSI data directory is missing"
         end if
         call LIS_verify(status)
         call ESMF_AttributeSet(OBS_State(n), "Data Directory",&
              USAFSIobsdir, rc=status)
         call LIS_verify(status)
      end do ! n

      call ESMF_ConfigFindLabel(LIS_config, "USAFSI netcdf filename prefix:",&
           rc=status)
      do n = 1, LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config, USAFSI_infile_name,&
              rc=status)
         if (status .ne. ESMF_SUCCESS)then
            write(LIS_logunit,*) "[ERR] USAFSI netcdf filename prefix is missing"
            call LIS_endrun()
         end if
         call LIS_verify(status)
         call ESMF_AttributeSet(OBS_State(n), "Netcdf filename prefix",&
              USAFSI_infile_name, rc=status)
         call LIS_verify(status)
      end do ! n

      do n=1,LIS_rc%nnest
         call ESMF_AttributeSet(OBS_State(n), "Data Update Status",&
              .false., rc=status)
         call LIS_verify(status)

         call ESMF_AttributeSet(OBS_State(n), "Data Update Time",&
              -99.0, rc=status)
         call LIS_verify(status)

         call ESMF_AttributeSet(OBS_State(n), "Data Assimilate Status",&
              .false., rc=status)
         call LIS_verify(status)
         
         call ESMF_AttributeSet(OBS_State(n), "Number Of Observations",&
              LIS_rc%obs_ngrid(k), rc=status)
         call LIS_verify(status)
      end do ! n

      write(LIS_logunit,*)'[INFO] Read USAFSI data specifications'

      ! Create the array containers that will contain the observations
      ! and the perturbations
      do n = 1, LIS_rc%nnest
         write(unit=temp, fmt='(i2.2)') 1
         read(unit=temp, fmt='(2a1)') vid

         obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_obsVecGrid(n,k),&
            name="Observation"//vid(1)//vid(2), rc=status)
         call LIS_verify(status)

         ! Perturbations state
         write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
              trim(LIS_rc%obsattribfile(k))
         ftn = LIS_getNextUnitNumber()
         open(ftn, file=trim(LIS_rc%obsattribfile(k)), status='old')
         read(ftn,*)
         read(ftn,*) LIS_rc%nobtypes(k)
         read(ftn,*)

         allocate(vname(LIS_rc%nobtypes(k)))
         allocate(varmax(LIS_rc%nobtypes(k)))
         allocate(varmin(LIS_rc%nobtypes(k)))

         do i = 1, LIS_rc%nobtypes(k)
            read(ftn, fmt='(a40)') vname(i)
            read(ftn,*) varmin(i),varmax(i)
            write(LIS_logunit,*) '[INFO] ',vname(i), varmin(i), varmax(i)
         enddo ! i
         call LIS_releaseUnitNumber(ftn)  
         
         allocate(ssdev(LIS_rc%obs_ngrid(k)))

         if (trim(LIS_rc%perturb_obs(k)) .ne. "none") then 
            allocate(obs_pert%vname(1))
            allocate(obs_pert%perttype(1))
            allocate(obs_pert%ssdev(1))
            allocate(obs_pert%stdmax(1))
            allocate(obs_pert%zeromean(1))
            allocate(obs_pert%tcorr(1))
            allocate(obs_pert%xcorr(1))
            allocate(obs_pert%ycorr(1))
            allocate(obs_pert%ccorr(1,1))

            call LIS_readPertAttributes(1, LIS_rc%obspertAttribfile(k),&
                 obs_pert)

            ! Set obs err to be uniform (will be rescaled later for each
            ! grid point)
            ssdev = obs_pert%ssdev(1)

            pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                 grid=LIS_obsEnsOnGrid(n,k), &
                 name="Observation"//vid(1)//vid(2),&
                 rc=status)
            call LIS_verify(status)
            
            ! Initialize perturbations to be zero
            call ESMF_FieldGet(pertField(n), localDE=0, farrayPtr=obs_temp,&
                 rc=status)
            call LIS_verify(status)
            obs_temp(:,:) = 0 

            call ESMF_AttributeSet(pertField(n), "Perturbation Type",&
                 obs_pert%perttype(1), rc=status)
            call LIS_verify(status)

            if (LIS_rc%obs_ngrid(k) .gt. 0) then 
               call ESMF_AttributeSet(pertField(n), &
                    "Standard Deviation",&
                    ssdev, itemCount=LIS_rc%obs_ngrid(k), rc=status)
               call LIS_verify(status)
            endif

            call ESMF_AttributeSet(pertField(n), "Std Normal Max",&
                 obs_pert%stdmax(1), rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(pertField(n), "Ensure Zero Mean",&
                 obs_pert%zeromean(1), rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(pertField(n), &
                 "Temporal Correlation Scale",&
                 obs_pert%tcorr(1), rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(pertField(n), "X Correlation Scale",&
                 obs_pert%xcorr(1), rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(pertField(n), "Y Correlation Scale",&
               obs_pert%ycorr(1), rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(pertField(n), &
               "Cross Correlation Strength",&
               obs_pert%ccorr(1,:), itemCount=1, rc=status)
            call LIS_verify(status)

            call ESMF_StateAdd(OBS_Pert_State(n), (/pertField(n)/), rc=status)
            call LIS_verify(status)

         end if

         deallocate(vname)
         deallocate(varmax)
         deallocate(varmin)
         deallocate(ssdev)

         ! The USAFSI grid information will be read in from the LDT
         ! netCDF file.  For now, put in dummy data.
         USAFSI_obs(n)%gridDesco(:) = 0         
         USAFSI_obs%mo1 = 0
         USAFSI_obs%nc_ldt = 0
         USAFSI_obs%nr_ldt = 0
         USAFSI_obs%nc_lis = 0
         USAFSI_obs%nr_lis = 0
         USAFSI_obs%mi = 0
         
         call LIS_registerAlarm("USAFSI read alarm", 21600.0, 21600.0)         
         call ESMF_StateAdd(OBS_State(n), (/obsField(n)/), rc=status)
         call LIS_verify(status)
         
      end do ! n

      write(LIS_logunit,*) &
           '[INFO] Created the States to hold the observations data'
   end subroutine USAFSIobs_setup

   ! Dummy subroutine to finalize work with USAFSI data
   subroutine USAFSIobs_finalize()
      implicit none
   end subroutine USAFSIobs_finalize

end module USAFSIobs_Mod
