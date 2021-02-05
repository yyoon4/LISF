!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_ASOS
! \label{get_ASOS}
!
! !REVISION HISTORY:
! 21 Jun 2013: Soni Yatheendradas; based on new TRMM 3B42V6 code that has
!              changes from earlier code to avoid (a) alternate file skip,
!              (b) jump to previous day TRMM, and
!              (c) absence of rain rate weighting
! 06 Sep 2017: Minwook Kim;  based on new TRMM 3B42V7 code that has
!              changes
!
! !INTERFACE:
subroutine get_ASOS(n, findex)
! !USES:
  use LIS_coreMod, only           : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only        : LIS_time2date, LIS_tick, LIS_get_nstep, &
                                    LIS_isAlarmRinging ! SY
  use LIS_logMod, only            : LIS_logunit, LIS_endrun
  use ASOS_forcingMod, only : ASOS_struc
  use LIS_metforcingMod, only     : LIS_forc ! SY

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Opens, reads, and interpolates 12-hrly ASOS forcing.
!  At the beginning of a simulation, the code
!  reads the most recent past data (nearest 6 hour interval), and
!  the nearest future data. These two datasets are used to
!  temporally interpolate the data to the current model timestep.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the ASOS data times
!  \item[ASOSfile](\ref{ASOSfile}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_ASOS](\ref{read_ASOS}) \newline
!      Interpolates ASOS data to LIS grid
!  \end{description}
!EOP

!==== Local Variables=======================
  integer :: ferror_ASOS      ! Error flag for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1,ts1                ! SY: Time parameters for start and end of current LDAS time step
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2               ! SY: Time parameters for ASOS data time nearest to start of model time step
  integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3, ts3               ! SY: Time parameters for ASOS data time nearest to end of model time step
  real    :: gmt1, gmt2, gmt3 ! SY ,kgmt3, mgmt3
  character(len=80) :: name ! Filename variables for precip data sources
  real*8 :: LIS_timeAtTStepStart ! SY
  real*8 :: LIS_timeAtTStepEnd ! SY
  integer :: order, kk
  logical             :: alarmCheck ! SY

!=== End Variable Definition =======================

!  Disable alarm.  It throws off the timing of reading the
!  next data set, and it throws off copying the bookends.
!  alarmCheck = LIS_isAlarmRinging(LIS_rc,"ASOS alarm") ! SY
!  if(alarmCheck) then ! SY

 !------------------------------------------------------------------------
 ! SY : Determine start and end time of model time step
 !------------------------------------------------------------------------
   yr1 = LIS_rc%yr  !current time, ! SY: i.e., end of model time step
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr
   mn1 = LIS_rc%mn
   ss1 = LIS_rc%ss
   ts1 = 0
   call LIS_tick( ASOS_struc(n)%LIS_timeAtTStepEnd, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
   ts1 = (-1)*LIS_rc%nts(n) ! SY: for start of model time step
   call LIS_tick( ASOS_struc(n)%LIS_timeAtTStepStart, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
 
 !------------------------------------------------------------------------
 ! SY: ASOS product times.
 !     I get these times based on a very simple principle: modify the
 !     time step start/end by adding 90 min. Whichever ASOS data stamp the
 !     modified start/end coincides with or immediately comes after, will be
 !     the ASOS data nearest to the original time step start/end. Hence,
 !     values of that ASOS data are assigned to that original start/end as
 !     preparation for weighting data in the time_interp* subroutine.
 !------------------------------------------------------------------------
 
   yr2 = LIS_rc%yr
   mo2 = LIS_rc%mo
   da2 = LIS_rc%da
   hr2 = LIS_rc%hr
   mn2 = LIS_rc%mn
   ss2 = LIS_rc%ss
   ts2 = (-1)*LIS_rc%nts(n)
   call LIS_tick( LIS_timeAtTStepStart, doy2, gmt2, &
        yr2, mo2, da2, hr2, mn2, ss2, real(ts2))
   ! SY: Now start calculations for ASOS data time nearest to start of model time step
   ASOS_struc(n)%ASOSyr_TStepStart = yr2
   ASOS_struc(n)%ASOSmo_TStepStart = mo2
   ASOS_struc(n)%ASOSda_TStepStart = da2
   ASOS_struc(n)%ASOShr_TStepStart = 12*(hr2/12)
   mn2 = 0
   ss2 = 0
   ts2 = 0
   call LIS_tick( ASOS_struc(n)%ASOStime_TStepStart, doy2, gmt2, &
        ASOS_struc(n)%ASOSyr_TStepStart, &
        ASOS_struc(n)%ASOSmo_TStepStart, &
        ASOS_struc(n)%ASOSda_TStepStart, &
        ASOS_struc(n)%ASOShr_TStepStart, mn2, ss2, real(ts2))
 
   yr3 = LIS_rc%yr
   mo3 = LIS_rc%mo
   da3 = LIS_rc%da
   hr3 = LIS_rc%hr
   mn3 = LIS_rc%mn
   ss3 = LIS_rc%ss
   ts3 = 0
   call LIS_tick( LIS_timeAtTStepEnd, doy3, gmt3, &
        yr3, mo3, da3, hr3, mn3, ss3, real(ts3))
   ! SY: Now start calculations for ASOS data time nearest to end of model time step
   ASOS_struc(n)%ASOSyr_TStepEnd = yr3
   ASOS_struc(n)%ASOSmo_TStepEnd = mo3
   ASOS_struc(n)%ASOSda_TStepEnd = da3
   ASOS_struc(n)%ASOShr_TStepEnd = 12*(hr3/12)
   mn3 = 0
   ss3 = 0
   ts3 = 12*60*60
   call LIS_tick( ASOS_struc(n)%ASOStime_TStepEnd, doy3, gmt3, &
        ASOS_struc(n)%ASOSyr_TStepEnd, &
        ASOS_struc(n)%ASOSmo_TStepEnd, &
        ASOS_struc(n)%ASOSda_TStepEnd, &
        ASOS_struc(n)%ASOShr_TStepEnd, mn3, ss3, real(ts3))
 
 !------------------------------------------------------------------------
 ! Check for and get ASOS observed Precipitation data
 !------------------------------------------------------------------------
 
   ! SY: Now update ASOS data corresponding to start of model time step if required
   if  (LIS_get_nstep(LIS_rc, n).eq. 1 .or. LIS_rc%rstflag(n) .eq. 1) then

      if (LIS_rc%nts(n) .ge. (12*60-1)*60) then ! almost 12 hrs
         write(LIS_logunit,*) 'LIS time step should be < 12 hrs!'
         write(LIS_logunit,*) 'ASOS reader functionality for >= 12 hrs not written yet'
         write(LIS_logunit,*) 'will involve weights combining > 2 ASOS data files'
         write(LIS_logunit,*) 'Program stopping ... '
         call LIS_endrun()
      endif

     ASOS_struc(n)%metdata1 = LIS_rc%udef
     ASOS_struc(n)%metdata2 = LIS_rc%udef
     ferror_ASOS = 0
     write(LIS_logunit, *) 'ASOS yr,mo,da,hr for time step start:', &
                           ASOS_struc(n)%ASOSyr_TStepStart, &
                           ASOS_struc(n)%ASOSmo_TStepStart, &
                           ASOS_struc(n)%ASOSda_TStepStart, &
                           ASOS_struc(n)%ASOShr_TStepStart
     order = 1
     do kk= ASOS_struc(n)%st_iterid, ASOS_struc(n)%en_iterid
       call ASOSfile( name, n, ASOS_struc(n)%ASOSyr_TStepStart, &
                            ASOS_struc(n)%ASOSmo_TStepStart, &
                            ASOS_struc(n)%ASOSda_TStepStart, &
                            ASOS_struc(n)%ASOShr_TStepStart )
       write(LIS_logunit, *)'Getting new ASOS observation precip data:', name
       call read_ASOS(n, kk, name, findex, order, ferror_ASOS)
     end do
   elseif (.NOT. ((ASOS_struc(n)%ASOSyr_TStepStart .EQ. &
                 ASOS_struc(n)%ASOSyr_TStepStart_Previous) .AND. &
                (ASOS_struc(n)%ASOSmo_TStepStart .EQ. &
                 ASOS_struc(n)%ASOSmo_TStepStart_Previous) .AND. &
                (ASOS_struc(n)%ASOSda_TStepStart .EQ. &
                 ASOS_struc(n)%ASOSda_TStepStart_Previous) .AND. &
                (ASOS_struc(n)%ASOShr_TStepStart .EQ. &
                 ASOS_struc(n)%ASOShr_TStepStart_Previous))) then
       write(LIS_logunit, *) 'ASOS yr,mo,da,hr for time step start:', &
                             ASOS_struc(n)%ASOSyr_TStepStart, &
                             ASOS_struc(n)%ASOSmo_TStepStart, &
                             ASOS_struc(n)%ASOSda_TStepStart, &
                             ASOS_struc(n)%ASOShr_TStepStart
       write(LIS_logunit, *)'Values from time step end assigned'
       ASOS_struc(n)%metdata1 = ASOS_struc(n)%metdata2 ! SY: i.e., equal to LIS_forc(n,findex)%metdata2 before possible change in LIS_forc(n,findex)%metdata2 below
   endif
 
   ! SY: Now update ASOS data corresponding to end of model time step if required.
   if (.NOT. ((ASOS_struc(n)%ASOSyr_TStepEnd .EQ. &
               ASOS_struc(n)%ASOSyr_TStepEnd_Previous) .AND. &
              (ASOS_struc(n)%ASOSmo_TStepEnd .EQ. &
               ASOS_struc(n)%ASOSmo_TStepEnd_Previous) .AND. &
              (ASOS_struc(n)%ASOSda_TStepEnd .EQ. &
               ASOS_struc(n)%ASOSda_TStepEnd_Previous) .AND. &
              (ASOS_struc(n)%ASOShr_TStepEnd .EQ. &
               ASOS_struc(n)%ASOShr_TStepEnd_Previous))) then
     ferror_ASOS = 0
     write(LIS_logunit, *) 'ASOS yr,mo,da,hr for time step end:', &
                           ASOS_struc(n)%ASOSyr_TStepEnd, &
                           ASOS_struc(n)%ASOSmo_TStepEnd, &
                           ASOS_struc(n)%ASOSda_TStepEnd, &
                           ASOS_struc(n)%ASOShr_TStepEnd
     order = 2
     do kk= ASOS_struc(n)%st_iterid, ASOS_struc(n)%en_iterid
       call ASOSfile( name, n, ASOS_struc(n)%ASOSyr_TStepEnd, &
                            ASOS_struc(n)%ASOSmo_TStepEnd, &
                            ASOS_struc(n)%ASOSda_TStepEnd, &
                            ASOS_struc(n)%ASOShr_TStepEnd )
       write(LIS_logunit, *)'Getting new ASOS observation precip data:', name
       call read_ASOS(n, kk, name, findex, order, ferror_ASOS)
     end do
   endif
 
   ! SY: Begin reassigning *_Previous values for use in next get* subroutine call
   ASOS_struc(n)%ASOSyr_TStepStart_Previous = &
                              ASOS_struc(n)%ASOSyr_TStepStart
   ASOS_struc(n)%ASOSmo_TStepStart_Previous = &
                              ASOS_struc(n)%ASOSmo_TStepStart
   ASOS_struc(n)%ASOSda_TStepStart_Previous = &
                              ASOS_struc(n)%ASOSda_TStepStart
   ASOS_struc(n)%ASOShr_TStepStart_Previous = &
                              ASOS_struc(n)%ASOShr_TStepStart
 
   ASOS_struc(n)%ASOSyr_TStepEnd_Previous = &
                              ASOS_struc(n)%ASOSyr_TStepEnd
   ASOS_struc(n)%ASOSmo_TStepEnd_Previous = &
                              ASOS_struc(n)%ASOSmo_TStepEnd
   ASOS_struc(n)%ASOSda_TStepEnd_Previous = &
                              ASOS_struc(n)%ASOSda_TStepEnd
   ASOS_struc(n)%ASOShr_TStepEnd_Previous = &
                              ASOS_struc(n)%ASOShr_TStepEnd
   ! SY: End reassigning *_Previous values for use in next get* subroutine call
 
 ! J.Case (4/19/2013) -- test print of the suppdata2 array.
 ! write (96,*) suppdata2(1,:)

!  end if ! SY: if(alarmCheck) then
end subroutine get_ASOS

!BOP
! !ROUTINE: ASOSfile
! \label{ASOSfile}
!
! !INTERFACE:
subroutine ASOSfile( name, n, yr, mo, da, hr)

! !DESCRIPTION: This subroutine puts together ASOS file name for
!               12 hour file intervals

! !USES:
  use ASOS_forcingMod, only : ASOS_struc

!EOP
  implicit none

  integer, intent(in) :: n

!==== Local Variables=======================

  character(len=80) :: name, ASOSdir
  character*160 temp
  integer :: yr, mo, da, hr
  integer :: i, j
  integer :: uyr, umo, uda, uhr, umn, uss
  integer :: original

  uyr = mod(yr, 100)          ! 1980 -> 80
  umo = mo
  uda = da
  uhr = 12*(hr/12)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

   ASOSdir = ASOS_struc(n)%ASOSdir

   name=''        !clean it
   temp=''

   original = 2
   if (original .eq. 1) then     !  1. original: /abc/3B42.980131.12.6.precipitation
     write(temp, '(a, a, I4, I2.2, a, 3I2.2, a, I2, a)') ASOSdir, '/', yr, mo, '/ASOS.', uyr, umo, uda, '.', &
          uhr,  '.12.precipitation'
   else                          !  2. renamed: ASOS.2005110809
     write(temp, '(a, a, I4, I2.2, a, I4, 3I2.2, a)') ASOSdir, '/', yr, mo, '/', yr, umo, uda, uhr, '_interp_rain.dat'
   end if

  !strip off the spaces
  j = 1
  Do i=1, len(temp)
   if( temp(i:i) .ne. ' ' ) then
      name(j:j)=temp(i:i)
      j=j+1
   end if
  End Do
end subroutine ASOSfile


