!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noah33_lsmMod
!BOP
!
! !MODULE: noah33_lsmMod
!
! !DESCRIPTION:
!  
! This module provides the definition of derived data type used to 
! control the operation of Noah3.3 LSM. It also provides the entry method
! for the initialization of Noah3.3-specific variables. The derived
!  data type {\tt noah3.3\_struc} includes the variables that specify
!  the runtime options and other control variables as described below:
!
! \begin{description}
!  \item[rfile]
!    name of the noah3.3 restart file
!  \item[vfile]
!    name of the static vegetation parameter table
!  \item[sfile]
!    name of the soil parameter table
!  \item[gfile]
!    name of the general parameter table
!  \item[useptf]
!    whether PTFs should be used to map soil properties
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[soilscheme]
!    type of soil mapping scheme usd (1-zobler,2-statsgo)
!  \item[nslay]
!    number of soil layers
!  \item[nvegp]
!    number of static vegetation parameters in the table
!  \item[nstxts]
!    number of soil texture classes in the classification scheme
!  \item[varid]
!    NETCDF ids for variables (used for netcdf output)
!  \item[inittemp]
!    initial soil temperatures for a cold start run
!  \item[initsm]
!    initial total soil moisture for a cold start run
!  \item[initsmliq]
!    initial liquid soil moisture for a cold start run
!  \item[rstInterval]
!   restart writing interval
!  \item[zh]
!   reference height of T and q forcing
!  \item[zm]
!   reference height of u and v forcing
!  \item[forcing\_ch]
!   flag indicating whether to use aerodynamic conductance field from forcing
!  \item[lyrthk]
!   thickness of soil layers
!  \item[noah]
!   Noah3.3 LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
! Apr 2003; Sujay Kumar, Initial Code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!
! !USES:        
  use noah33_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah33_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noah33_struc
!EOP
  type, public ::  noah33_type_dec 

     character(len=LIS_CONST_PATH_LEN) :: rfile
     character(len=LIS_CONST_PATH_LEN) :: vfile
     character(len=LIS_CONST_PATH_LEN) :: sfile
     character(len=LIS_CONST_PATH_LEN) :: gfile
     integer                    :: useptf
     integer                    :: soilscheme
     integer                    :: nslay
     integer                    :: varid(21)
     integer                    :: usedsoilmap
     integer                    :: usedrootmap
     integer                    :: fixedvegtype
     integer                    :: fixedsoiltype
     integer                    :: fixedslopetype
     real                       :: fixedtbot
     real                       :: fixedmxsnalb
     real                       :: shdfac_monthly(12)
     real                       :: albedo_monthly(12)
     real                       :: z0brd_monthly(12)
     integer                    :: iz0tlnd
     integer                    :: sfcdifoption
     real                       :: rstInterval
     real                       :: zh
     real                       :: zm
     integer                    :: forcing_ch
     integer                    :: forcing_cm
     real, allocatable              :: lyrthk(:)
     real, allocatable              :: inittemp(:)
     real, allocatable              :: initsm(:)
     real, allocatable              :: initsmliq(:)
     real                       :: initskintemp
     real                       :: initcanopywater
     real                       :: initsnowdepth
     real                       :: initsnowequiv
     real :: ztmax2
     real :: dzeta2
     real :: psih2(10001)
     real :: psim2(10001)

! these flags check to see if these parameters are set from OPTUE
! (and subsequently should not be modified)
     integer :: z0brd_upd
     integer :: lai_upd
     integer :: embrd_upd
     integer :: alb_upd
     integer :: optStartFlag

     integer          :: lucats
     real, allocatable    :: shdtbl(:)
     real, allocatable    :: nrotbl(:)
     real, allocatable    :: rstbl(:)
     real, allocatable    :: rgltbl(:)
     real, allocatable    :: hstbl(:)
     real, allocatable    :: snuptbl(:)
     real, allocatable    :: maxalb(:)
     real, allocatable    :: emissmin(:)
     real, allocatable    :: emissmax(:)
     real, allocatable    :: laimin(:)
     real, allocatable    :: laimax(:)
     real, allocatable    :: albmin(:)
     real, allocatable    :: albmax(:)
     real, allocatable    :: z0min(:)
     real, allocatable    :: z0max(:)
     real, allocatable    :: cziltbl(:)
     real, allocatable    :: cmcmaxtbl(:)

     integer                  :: forc_count
     real                     :: ts
     type(noah33dec), allocatable :: noah(:)
  end type noah33_type_dec

  type(noah33_type_dec), allocatable :: noah33_struc(:)
  SAVE
contains
!BOP
! 
! !ROUTINE: noah33_lsm_ini
! \label{noah33_lsm_ini}
! 
! !INTERFACE:
  subroutine noah33_lsm_ini()
! !USES:
    use LIS_coreMod,      only : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_logMod,       only : LIS_verify
! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for Noah3.3-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for Noah3.3 from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[noah33\_readcrd](\ref{noah33_readcrd}) \newline
!    reads the runtime options for Noah3.3 LSM
!  \end{description}
!EOP
    implicit none
    integer                 :: i,n,t
    character*3             :: fnest
    integer                 :: status

    allocate(noah33_struc(LIS_rc%nnest))
    call noah33_readcrd()
    do n=1,LIS_rc%nnest
       allocate(noah33_struc(n)%noah(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          allocate(noah33_struc(n)%noah(i)%stc(noah33_struc(n)%nslay))
          allocate(noah33_struc(n)%noah(i)%smc(noah33_struc(n)%nslay))
          allocate(noah33_struc(n)%noah(i)%sh2o(noah33_struc(n)%nslay))
          allocate(noah33_struc(n)%noah(i)%relsmc(noah33_struc(n)%nslay))
#ifdef PARFLOW
          allocate(noah33_struc(n)%noah(i)%wtrflx(noah33_struc(n)%nslay))
#endif
       enddo
       noah33_struc(n)%forc_count = 0 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah33_struc(n)%noah(i)%tair = 0 
          noah33_struc(n)%noah(i)%qair = 0 
          noah33_struc(n)%noah(i)%swdown = 0 
          noah33_struc(n)%noah(i)%lwdown = 0
          noah33_struc(n)%noah(i)%uwind = 0
          noah33_struc(n)%noah(i)%vwind = 0 
          noah33_struc(n)%noah(i)%psurf = 0 
          noah33_struc(n)%noah(i)%rainf = 0 
          noah33_struc(n)%noah(i)%snowf = 0 
          noah33_struc(n)%noah(i)%rainf_c = 0 
          noah33_struc(n)%noah(i)%ch = 0
          
          noah33_struc(n)%noah(i)%shdfac = 0 
          noah33_struc(n)%noah(i)%alb = 0 
          noah33_struc(n)%noah(i)%z0 = 0 

          noah33_struc(n)%noah(i)%stc = LIS_rc%udef
          noah33_struc(n)%noah(i)%smc = LIS_rc%udef
          noah33_struc(n)%noah(i)%sh2o = LIS_rc%udef
          noah33_struc(n)%noah(i)%relsmc = LIS_rc%udef
#if WRF_HYDRO
          noah33_struc(n)%noah(i)%sfhead1rt = 0
          noah33_struc(n)%noah(i)%infxs1rt = 0
          noah33_struc(n)%noah(i)%soldrain1rt = 0
#endif
#if PARFLOW
          noah33_struc(n)%noah(i)%wtrflx = LIS_rc%udef
#endif
          
       enddo

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, noah33_struc(n)%ts)
       
       write(fnest,'(i3.3)') n
       call LIS_registerAlarm("Noah33 model alarm "//trim(fnest),&
            noah33_struc(n)%ts,&
            noah33_struc(n)%ts)

       call LIS_registerAlarm("Noah33 restart alarm "//trim(fnest),&
            noah33_struc(n)%ts,&
            noah33_struc(n)%rstInterval)

       ! Initialize min/max values to implausible values.
       noah33_struc(n)%noah(:)%tair_agl_min = 999.0


       noah33_struc(n)%z0brd_upd = 0 
       noah33_struc(n)%lai_upd = 0 
       noah33_struc(n)%embrd_upd = 0 
       noah33_struc(n)%alb_upd = 0 

       noah33_struc(n)%optStartFlag = 1 
!uninitialized variables: 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah33_struc(n)%noah(i)%cqs2 = LIS_rc%udef
          noah33_struc(n)%noah(i)%qsfc = LIS_rc%udef
          noah33_struc(n)%noah(i)%sca = 0.0
       enddo


       LIS_sfmodel_struc(n)%nsm_layers = noah33_struc(n)%nslay
       LIS_sfmodel_struc(n)%nst_layers = noah33_struc(n)%nslay
       allocate(LIS_sfmodel_struc(n)%lyrthk(noah33_struc(n)%nslay))
       do i = 1,noah33_struc(n)%nslay
          LIS_sfmodel_struc(n)%lyrthk(i) = noah33_struc(n)%lyrthk(i)*100.0
       enddo
       LIS_sfmodel_struc(n)%ts = noah33_struc(n)%ts

!!!!       !initialize state vars; will be overridden later if restart mode
!!!!       call noah33_coldstart(LIS_rc%lsm_index)
!!!!       ! Set initial q1 to a negative value - D. Mocko
!!!!       ! Will be set to q2 for the first timestep in noah33_main
!!!!       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!!!!          noah33_struc(n)%noah(i)%q1 = -0.5
!!!!       enddo
    enddo

!moved the coldstart setting to init to help the resetting of 
! initial conditions for OPT cases (if initial conditions are part of 
! the parameteters being optimized. 
    
    call noah33_coldstart(LIS_rc%lsm_index)
     ! Set initial q1 to a negative value - D. Mocko
! Will be set to q2 for the first timestep in noah33_main
    do n=1,LIS_rc%nnest
       do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah33_struc(n)%noah(t)%q1 = -0.5
          noah33_struc(n)%noah(t)%wchange_prev = 0 
       enddo
    enddo

  end subroutine noah33_lsm_ini

end module noah33_lsmMod
