!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE:  get_ecmwf
!  \label{get_ecmwf}
!
! !REVISION HISTORY:
!  18 Jun 2003: Urszula Jambor; original code based on getreanlecmwf.F90
!  20 Feb 2006: Sujay Kumar; Modified with nesting options
! 
! !INTERFACE:
subroutine get_ecmwf(n,findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_logMod,         only : LIS_logunit
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_tick
  use ecmwf_forcingMod,   only : ecmwf_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, 1/4 degree 
!  ECMWF forcing. At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
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
!    call to advance or retract time
!  \item[read\_ecmwf](\ref{read_ecmwf}) \newline
!    call to read the ECMWF data and perform spatial interpolation 
!   \item[read\_ecmwf\_elev](\ref{read_ecmwf_elev}) \newline
!    reads the native elevation of the ECMWF data to be used
!    for topographic adjustments to the forcing
!  \end{description}
!EOP
  integer, parameter :: ndays = 10  ! # days to look back for forcing data
  integer :: try, ferror
  integer :: c,f,order
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  real*8  :: time1,time2,dumbtime1,dumbtime2
  character(len=LIS_CONST_PATH_LEN) :: avgfilename1, instfilename, avgfilename2
  real    :: gmt1,gmt2,ts1,ts2
  integer :: movetime      ! 1=move time 2 data into time 1
  integer :: nforce  ! # forcing variables
  integer :: nstep

  ecmwf_struc(n)%findtime1=0
  ecmwf_struc(n)%findtime2=0
  movetime=0

!  nforce = ecmwf_struc(n)%nmif
  nforce = LIS_rc%met_nf(findex)
  nstep = LIS_get_nstep(LIS_rc,n)

  !=== Determine Required Data Times (The previous hour & the future hour)
  yr1=LIS_rc%yr    !Previous Hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=3*((LIS_rc%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

  yr2=LIS_rc%yr    !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=3*((LIS_rc%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2) 

!=== Get elevation file
  if(trim(LIS_rc%met_ecor(findex)).ne."none") then
      if ( time1 > ecmwf_struc(n)%griduptime1 .and. &
           time1 < ecmwf_struc(n)%griduptime2 .and. & 
          ecmwf_struc(n)%gridchange1) then 

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 2003 -- 2/2006'
        call read_ecmwf_elev(n, findex, 1)
        ecmwf_struc(n)%gridchange1 = .false.

      elseif ( time1 > ecmwf_struc(n)%griduptime2 .and. & 
               time1 < ecmwf_struc(n)%griduptime3 .and. &
               ecmwf_struc(n)%gridchange2) then 

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 2/2006 - 6/2008'
        call read_ecmwf_elev(n,findex, 2)
        ecmwf_struc(n)%gridchange2 = .false.

      elseif ( time1 > ecmwf_struc(n)%griduptime3 .and. & 
               time1 < ecmwf_struc(n)%griduptime4 .and. &
               ecmwf_struc(n)%gridchange3) then 

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 6/2008 - 3/2009'
        call read_ecmwf_elev(n,findex, 3)
        ecmwf_struc(n)%gridchange3 = .false.

      elseif ( time1 > ecmwf_struc(n)%griduptime4 .and. & 
               time1 < ecmwf_struc(n)%griduptime5 .and. &
               ecmwf_struc(n)%gridchange4) then 

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 3/2009 - 9/2009'
        call read_ecmwf_elev(n,findex, 4)
        ecmwf_struc(n)%gridchange4 = .false.

      elseif ( time1 > ecmwf_struc(n)%griduptime5 .and. & 
               time1 < ecmwf_struc(n)%griduptime6 .and. &
               ecmwf_struc(n)%gridchange5) then 

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 9/2009 - 1/2010'
        call read_ecmwf_elev(n,findex, 5)
        ecmwf_struc(n)%gridchange5 = .false.

      elseif ( time1 > ecmwf_struc(n)%griduptime6 .and. & 
               time1 < ecmwf_struc(n)%griduptime7 .and. &
               ecmwf_struc(n)%gridchange6) then

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 1/2010 - 5/2011'
        call read_ecmwf_elev(n,findex, 6)
        ecmwf_struc(n)%gridchange6 = .false.

      elseif ( time1 > ecmwf_struc(n)%griduptime7 .and. & 
               ecmwf_struc(n)%gridchange7) then 

        write(LIS_logunit,*) 'MSG: get_ecmwf -- changing elev to 5/2011 -'
        call read_ecmwf_elev(n,findex, 7)
        ecmwf_struc(n)%gridchange7 = .false.

      endif
     endif
  !-----------------------------------------------------------------
  ! Use these if need to roll back time.
  !-----------------------------------------------------------------
  dumbtime1 = time1
  dumbtime2 = time2
  write(112,*)'nstep:',nstep,ecmwf_struc(n)%ecmwftime1,ecmwf_struc(n)%ecmwftime2

  !=== Check if time interval boundary was crossed
  if(LIS_rc%time.gt.ecmwf_struc(n)%ecmwftime2) then
     movetime=1
     ecmwf_struc(n)%findtime2=1
  endif

  if (LIS_get_nstep(LIS_rc,n) .eq. 1.or.LIS_rc%rstflag(n).eq.1) then ! beginning of the run
     ecmwf_struc(n)%findtime1=1
     ecmwf_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif

  !=== Establish fmodeltime1
  if (ecmwf_struc(n)%findtime1==1) then  !need to get new time1 from the past
     order=1   !Get data for glbdata1
     ferror = 0
     try = 0
     ts1 = -24*60*60
     do
        if ( ferror /= 0 ) then
           exit
        end if
        try = try+1
        call create_ecmwf_filename(ecmwf_struc(n)%ecmwfdir,&
             avgfilename1, avgfilename2, instfilename,&
             yr1,mo1,da1,hr1)
        write(LIS_logunit,*) 'opening Bookend1 (I) ',trim(instfilename)
        write(LIS_logunit,*) 'opening Bookend1 file1 (A) ',trim(avgfilename1)
        write(LIS_logunit,*) 'opening Bookend1 file2 (A) ',trim(avgfilename2)
     
        call read_ecmwf(order,n, findex, &
             avgfilename1, avgfilename2, instfilename, &
             yr1,mo1,da1,hr1,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ecmwf_struc(n)%ecmwftime1=time1
        else  !ferror still=0, so roll back one day
           call LIS_tick(dumbtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        end if
        if ( try > ndays ) then 
           write(LIS_logunit,*) 'ERROR: ECMWF data gap exceeds 10 days'
           STOP
        end if
     end do
  endif

  !=== Find new time2 value and tranfer time2 data to time1 
  if(movetime.eq.1) then
     ecmwf_struc(n)%ecmwftime1=ecmwf_struc(n)%ecmwftime2
     ecmwf_struc(n)%findtime2=1 !to ensure getting new time2 data
     do f=1,nforce
        do c=1,LIS_rc%ngrid(n)
           ecmwf_struc(n)%metdata1(f,c)=ecmwf_struc(n)%metdata2(f,c)
        enddo
     enddo
  endif  ! if movetime=1
  
  if(ecmwf_struc(n)%findtime2.eq.1) then ! need new time2 data
     order=2   !Get data for glbdata2
     ferror = 0
     try = 0
     ts2 = -24*60*60
     do
        if ( ferror /= 0 ) exit
        try = try+1
        call create_ecmwf_filename(ecmwf_struc(n)%ecmwfdir,&
             avgfilename1, avgfilename2, instfilename,&
             yr2,mo2,da2,hr2)
        write(LIS_logunit,*) 'opening Bookend2 (I) ',trim(instfilename)
        write(LIS_logunit,*) 'opening Bookend2 file1 (A) ',trim(avgfilename1)
        write(LIS_logunit,*) 'opening Bookend2 file2 (A) ',trim(avgfilename2)
        call read_ecmwf(order,n,findex, &
             avgfilename1, avgfilename2, instfilename, &
             yr2,mo2,da2,hr2,ferror)
        if ( ferror == 1 ) then !successfully retrieved forcing data
           ecmwf_struc(n)%ecmwftime2=time2
        else  !ferror still=0, so roll back one day
           call LIS_tick(dumbtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        end if
        if ( try > ndays ) then
           write(LIS_logunit,*)'ERROR: ECMWF data gap exceeds 10 days'
           STOP
        end if
     end do
  endif
  
end subroutine get_ecmwf
   
!BOP
! 
! !ROUTINE: create_ecmwf_filename
! \label{create_ecmwf_filename}
! 
! !INTERFACE: 
subroutine create_ecmwf_filename(dir,avgfilename1, avgfilename2, instfilename,&
     yr,mo,da,hr)
! !USES: 
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  character(len=*)   :: dir
  character(len=*)   :: avgfilename1
  character(len=*)   :: avgfilename2
  character(len=*)   :: instfilename
  integer            :: yr, mo, da, hr
! 
! !DESCRIPTION: 
!  creates the ECMWF filenames
! 
!EOP

  integer     :: remainder
  character :: cyr*4,cmo*2,cda*2,chr*2,fhr*2, chr3*2
  character :: ciyr*4, cimo*2, cida*2
  real*8 :: itime
  real :: igmt
  integer :: iyr,imo,ida,ihr,imn,iss,ts,idoy
  character(len=LIS_CONST_PATH_LEN) :: filename
  character(len=LIS_CONST_PATH_LEN) :: file1, file2


  !instantaneous files 
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo
  write(cda, '(i2.2)') da
  write(chr, '(i2.2)') hr
  remainder = modulo(hr,2)
  if (remainder==0) then
     !=== Use analysis field, rather than forecast
     fhr = chr
  else if ((hr==03).or.(hr==15)) then
     write(fhr, '(i2.2)') hr-3
  else if ((hr==09).or.(hr==21)) then
     write(fhr, '(i2.2)') hr-9
  end if
  filename = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
  instfilename = trim(dir)//'/'//cyr//cmo//'/'//trim(filename)

  !time averaged filenames
  !=== Establish file name
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo
  write(cda, '(i2.2)') da
  write(chr, '(i2.2)') hr
  if ((hr==09).or.(hr==21).or.(hr==06).or.(hr==18).or.(hr==12)) then
     write(chr3,'(i2.2)') hr-3
  endif
  if (hr==0) then
     !=== roll back one day for forecast initialization time
     iyr=yr;  imo=mo;  ida=da
     ihr=hr;  imn=0;   iss=0
     ts = -24*60*60
     call LIS_tick(itime,idoy,igmt,iyr,imo,ida,ihr,imn,iss,real(ts))
     write(ciyr, '(i4.4)') iyr
     write(cimo, '(i2.2)') imo
     write(cida, '(i2.2)') ida
     file2 = 'ecmwf.'//ciyr//cimo//cida//'12.'//cmo//cda//'00.1_4'
     file1 = 'ecmwf.'//ciyr//cimo//cida//'12.'//cimo//cida//'21.1_4'
  else if ((hr==03).or.(hr==15)) then
     write(fhr, '(i2.2)') hr-3
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'none'
  else if ((hr==06).or.(hr==18)) then
     write(fhr, '(i2.2)') hr-6
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  else if ((hr==09).or.(hr==21)) then
     write(fhr, '(i2.2)') hr-9
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  else if (hr==12) then
     write(fhr, '(i2.2)') hr-12
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  end if
  if (hr==0) then
     avgfilename2 = trim(dir)//'/'//ciyr//cimo//'/'//trim(file2)
     avgfilename1 = trim(dir)//'/'//ciyr//cimo//'/'//trim(file1)
  else
     avgfilename2 = trim(dir)//'/'//cyr//cmo//'/'//trim(file2)
     avgfilename1 = trim(dir)//'/'//cyr//cmo//'/'//trim(file1)
  end if

end subroutine create_ecmwf_filename



