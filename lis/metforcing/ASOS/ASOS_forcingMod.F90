!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ASOS_forcingMod
!BOP
! !MODULE: ASOS_forcingMod
!
! !REVISION HISTORY:
! 21 Jun 2013: Soni Yatheendradas; changes from earlier code to avoid
!              alternate file skip, jump to previous day ASOS and 
!              absence of rain rate weighting
! 06 Sep 2017: Minwook Kim; changes from earlier code
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  NASA TRMM 3B42 (Version 7) merged analysis of precipitation (TRMM 3B42V7).
!  TRMM 3B42V7 algorithm consists of two separate steps to produce a
!  quasi-global 0.25 degree 3-hourly precipitation analysis. The first step
!  uses the TRMM VIRS and TMI orbit data (TRMM products 1B01 and 2A12)
!  and the monthly TMI/TRMM Combined Instrument (TCI) calibration parameters
!  (from TRMM product 3B31) to produce monthly IR calibration parameters. The
!  second step uses these derived monthly IR calibration parameters to adjust
!  the merged-IR precipitation data, which consists of GMS, GOES-E, GOES-W,
!  Meteosat-7, Meteosat-5, and NOAA-12 data.
!
!  The implementation in LIS has the derived data type {\tt ASOS\_struc}
!  that includes the variables to specify the runtime options, and
!  the weights and neighbor information for spatial interpolation
!
!  They are desribed below:
! \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[ASOSdir]
!    Directory containing the input data
!  \item[ASOStime]
!    The nearest instance of the incoming
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES:
    implicit none
    PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_ASOS     !defines the native resolution of
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ASOS_struc

!EOP

  type, public :: ASOS_type_dec
     real                     :: ts 
     integer                  :: ncold
     integer                  :: nrold
     character*40             :: ASOSdir
     !real*8                   :: ASOStime ! SY
     real*8                   :: ASOStime_TStepStart ! SY
     integer                  :: ASOSyr_TStepStart ! SY
     integer                  :: ASOSmo_TStepStart ! SY
     integer                  :: ASOSda_TStepStart ! SY
     integer                  :: ASOShr_TStepStart ! SY
     integer                  :: ASOSyr_TStepStart_Previous ! SY
     integer                  :: ASOSmo_TStepStart_Previous ! SY
     integer                  :: ASOSda_TStepStart_Previous ! SY
     integer                  :: ASOShr_TStepStart_Previous ! SY
     real*8                   :: ASOStime_TStepEnd ! SY
     integer                  :: ASOSyr_TStepEnd ! SY
     integer                  :: ASOSmo_TStepEnd ! SY
     integer                  :: ASOSda_TStepEnd ! SY
     integer                  :: ASOShr_TStepEnd ! SY
     integer                  :: ASOSyr_TStepEnd_Previous ! SY
     integer                  :: ASOSmo_TStepEnd_Previous ! SY
     integer                  :: ASOSda_TStepEnd_Previous ! SY
     integer                  :: ASOShr_TStepEnd_Previous ! SY
     real*8                   :: LIS_timeAtTStepStart ! SY
     real*8                   :: LIS_timeAtTStepEnd ! SY
     !real*8                   :: griduptime1 ! SY
     !real*8                   :: griduptime2 ! SY
     !logical                  :: gridchange1 ! SY
     !logical                  :: gridchange2 ! SY
     integer                  :: mi
     integer, allocatable     :: n112(:,:)
     integer, allocatable     :: n122(:,:)
     integer, allocatable     :: n212(:,:)
     integer, allocatable     :: n222(:,:)
     real,    allocatable     :: w112(:,:),w122(:,:)
     real,    allocatable     :: w212(:,:),w222(:,:)

     integer           :: nIter, st_iterid,en_iterid  ! Forecast mode

     real, allocatable :: metdata1(:,:,:)
     real, allocatable :: metdata2(:,:,:)

  end type ASOS_type_dec

  type(ASOS_type_dec), allocatable :: ASOS_struc(:)
contains

!BOP
!
! !ROUTINE: init_ASOS
! \label{init_ASOS}
!
! !REVISION HISTORY:
! 11Dec2003: Sujay Kumar; Initial Specification 
! 25Aug2006: Yudong Tian; Implementation for TRMM 3B42 V7
! 06Sep2006: Minwook Kim; Implementation for ASOS 
!
! !INTERFACE:
  subroutine init_ASOS(findex)
! !USES:
    use LIS_coreMod, only: LIS_rc, LIS_config, LIS_domain
    !use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep ! SY
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep, &
                               LIS_tick, LIS_parseTimeString, &
                               LIS_registerAlarm ! SY
    use LIS_logMod,    only : LIS_logunit, LIS_endrun, &
                              LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use LIS_FORC_AttributesMod
    use LIS_forecastMod
    use ESMF

    implicit none
    
    integer,  intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for ASOS 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_ASOS](\ref{readcrd_ASOS}) \newline
!     reads the runtime options specified for ASOS data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    real :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n
    integer :: ts1,ts2 ! SY
    integer :: LIS_syr,LIS_smo,LIS_sda,LIS_shr,LIS_smn,LIS_sss ! SY

    real*8                   :: DataAvailabilityStartTime ! SY
    real*8                   :: LIS_StartTime ! SY

    integer             :: rc ! SY
    character (len=10):: time ! SY

    allocate(ASOS_struc(LIS_rc%nnest))

    ! SY: Begin LIS time check against ASOS data availability time
    yr1 = 1980
    mo1 = 01
    da1 = 01
    hr1 = 0
    mn1 = 0; ss1 = 0
    ts1 = (-1) * 12 * 60 * 60 ! For 12*60 min behind 1st available ASOS data
    call LIS_tick(DataAvailabilityStartTime,updoy,upgmt,yr1,mo1,&
         da1,hr1,mn1,ss1,real(ts1) )
    LIS_syr = LIS_rc%syr
    LIS_smo = LIS_rc%smo
    LIS_sda = LIS_rc%sda
    LIS_shr = LIS_rc%shr
    LIS_smn = LIS_rc%smn
    LIS_sss = LIS_rc%sss
    ts2 = 0
    call LIS_tick(LIS_StartTime,updoy,upgmt,LIS_syr,LIS_smo,&
         LIS_sda,LIS_shr,LIS_smn,LIS_sss,real(ts2))
    if (DataAvailabilityStartTime .gt. LIS_StartTime) then
       write(LIS_logunit,*) 'LIS start time is earlier than ASOS data availability time!'
       write(LIS_logunit,*) 'Program stopping ... '
       call LIS_endrun()
    endif
    ! SY: End LIS time check against ASOS data availability time

    call readcrd_ASOS()

    ! SY: Start for obtaining intended ASOS read time step
!    call ESMF_ConfigFindLabel(LIS_config,"ASOS timestep:",rc=rc)
    do n=1, LIS_rc%nnest
!       call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
!       call LIS_verify(rc,'ASOS timestep: not defined')

!       call LIS_parseTimeString(time,ASOS_struc(n)%ts)

       ASOS_struc(n)%ts = 12*60*60 ! SY: Because a LIS 12-hr timestep is typically valid from Z to Z, while ASOS 12-hr data validity is from 0.5 Z to 0.5 Z 
       call LIS_update_timestep(LIS_rc, n, ASOS_struc(n)%ts)

       call LIS_registerAlarm("ASOS alarm",&
            ASOS_struc(n)%ts,&
            ASOS_struc(n)%ts, alarm_offset=12*60*60)
    enddo

    LIS_rc%met_nf(findex) = 2 !number of met variables in ASOS forcing

    do n=1,LIS_rc%nnest
      if(LIS_rc%forecastMode.eq.1) then
        if(mod(LIS_rc%nensem(n),&
            LIS_forecast_struc(1)%niterations).ne.0) then
           write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
           write(LIS_logunit,*) '[ERR] of the number of iterations '
           write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
           write(LIS_logunit,*) '[ERR] niter  = ',LIS_forecast_struc(1)%niterations
           call LIS_endrun()
        endif

        ASOS_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
        ASOS_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
        ASOS_struc(n)%nIter = LIS_forecast_struc(1)%niterations

        allocate(ASOS_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))
        allocate(ASOS_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
             LIS_rc%met_nf(findex),&
             LIS_rc%ngrid(n)))

      else  ! Regular retrospective or non-forecast mode:

         ASOS_struc(n)%st_iterid = 1
         ASOS_struc(n)%en_iterId = 1
         ASOS_struc(n)%nIter = 1

         allocate(ASOS_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
         allocate(ASOs_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
      endif

      ASOS_struc(n)%metdata1 = 0
      ASOS_struc(n)%metdata2 = 0

      gridDesci = 0
      gridDesci(1) = 0
      gridDesci(2) = 33 
      gridDesci(3) = 45
      gridDesci(4) = 33. 
      gridDesci(5) = 124.
      gridDesci(6) = 128
      gridDesci(7) = 44. 
      gridDesci(8) = 132. 
      gridDesci(9) = 0.25
      gridDesci(10) = 0.25
      gridDesci(20) = 64

      ASOS_struc(n)%ncold = 33 
      ASOS_struc(n)%nrold = 45

      ! SY: Begin initialization of ASOS_struc time components to 0
      ASOS_struc(n)%ASOStime_TStepStart = 0.0
      ASOS_struc(n)%ASOSyr_TStepStart = 0
      ASOS_struc(n)%ASOSmo_TStepStart = 0
      ASOS_struc(n)%ASOSda_TStepStart = 0
      ASOS_struc(n)%ASOShr_TStepStart = 0
      ASOS_struc(n)%ASOSyr_TStepStart_Previous = 0
      ASOS_struc(n)%ASOSmo_TStepStart_Previous = 0
      ASOS_struc(n)%ASOSda_TStepStart_Previous = 0
      ASOS_struc(n)%ASOShr_TStepStart_Previous = 0
      ASOS_struc(n)%ASOStime_TStepEnd = 0.0
      ASOS_struc(n)%ASOSyr_TStepEnd = 0
      ASOS_struc(n)%ASOSmo_TStepEnd = 0
      ASOS_struc(n)%ASOSda_TStepEnd = 0
      ASOS_struc(n)%ASOShr_TStepEnd = 0
      ASOS_struc(n)%ASOSyr_TStepEnd_Previous = 0
      ASOS_struc(n)%ASOSmo_TStepEnd_Previous = 0
      ASOS_struc(n)%ASOSda_TStepEnd_Previous = 0
      ASOS_struc(n)%ASOShr_TStepEnd_Previous = 0
      ASOS_struc(n)%LIS_timeAtTStepStart = 0.0
      ASOS_struc(n)%LIS_timeAtTStepEnd = 0.0
      ! SY: End initialization of ASOS_struc time components to 0

       allocate(ASOS_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(ASOS_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

      if(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
         call conserv_interp_input(n,gridDesci,&
              ASOS_struc(n)%n112,ASOS_struc(n)%n122,&
              ASOS_struc(n)%n212,ASOS_struc(n)%n222,&
              ASOS_struc(n)%w112,ASOS_struc(n)%w122,&
              ASOS_struc(n)%w212,ASOS_struc(n)%w222)
      elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then ! SY
         ! SY: begin for nearest neighbor
         call neighbor_interp_input(n,gridDesci,&
              ASOS_struc(n)%n112)
      ! SY: end for nearest neighbor
      else
         write(LIS_logunit,*) 'This interpolation not defined for ASOS data'
         write(LIS_logunit,*) 'Program stopping ... '
         call LIS_endrun()
      endif

!       yr1 = 1998     !grid update time
!       mo1 = 01
!       da1 = 01
!       hr1 = 0
!       mn1 = 0; ss1 = 0
!       call LIS_date2time(ASOS_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,&
!            da1,hr1,mn1,ss1 )
!
!       yr1 = 2100     !grid update time
!       mo1 = 05
!       da1 = 31
!       hr1 = 0
!       mn1 = 0; ss1 = 0
!       call LIS_date2time(ASOS_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,&
!            da1,hr1,mn1,ss1 )
!
!       ASOS_struc(n)%gridchange1 = .true.
!       ASOS_struc(n)%gridchange2 = .false.

    enddo
  end subroutine init_ASOS

end module ASOS_forcingMod
                                                                        
