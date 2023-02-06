!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_mogrepsg_pcp
! \label{read_mogrepsg_pcp}
!
! !REVISION HISTORY:
! 26 Jan 2023: Yeosang Yoon, initial code
!
! !INTERFACE:
subroutine read_mogrepsg_pcp(n, m, findex, order, gribfile, rc)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use mogrepsg_forcingMod, only : mogrepsg_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none

!ARGUMENTS:
  integer,           intent(in)    :: n
  integer,           intent(in)    :: m       ! number of ensembles
  integer,           intent(in)    :: findex  ! forcing index
  integer,           intent(in)    :: order
  character(len=*),  intent(in)    :: gribfile

!DESCRIPTION:
!  For the given time, reads the precipitation forcing data from the
!  MOGREPS-G file, transforms into LIS forcing
!  parameters and interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    hourly instance, order=2, read the next hourly instance)
!  \item[gribfile]
!    name of the file to be read
!  \item[rc]
!    return error flag (0-fail, 1-success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[interp\_mogrepsg](\ref{interp_mogrepsg}) \newline
!    Performs spatial interpolation of GALEM-GE forecast data to the LIS grid
!  \end{description}

!EOP
  integer         :: ftn, igrib, ierr
  integer         :: center
  character*100   :: gtype
  integer         :: file_julhr
  integer         :: yr1, mo1, da1, hr1
  character*100   :: message     ( 20 )
  integer         :: iginfo      ( 40 )
  real            :: gridres_dlat, gridres_dlon
  integer         :: ifguess, jfguess
  integer         :: dataDate, dataTime
  integer         :: fc_hr

  real            :: prectot(LIS_rc%lnc(n),LIS_rc%lnr(n))   ! Total precipitation [kg/m^2] 

  integer, intent(out) :: rc

  ! Initialize return code to "no error".  We will change it below if
  ! necessary.
  rc = 0

#if (defined USE_GRIBAPI)

   call grib_open_file(ftn,trim(gribfile),'r',ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] Failed to open - '//trim(gribfile)
      call LIS_endrun()
   end if

   ! Read in the first grib record, unpack the header and extract
   ! section 1 and section 2 information.
   call grib_new_from_file(ftn,igrib,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] failed to read - '//trim(gribfile)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'centre',center,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grib_get: ' // &
           'centre in read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'gridType',gtype,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'gridtype in read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   if(trim(gtype).ne."regular_ll") then
      write(LIS_logunit,*)'[ERR] MOGREPS-G file not on lat/lon grid!'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'Ni',iginfo(1),ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: Ni read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'Nj',iginfo(2),ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: Nj in read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'jDirectionIncrementInDegrees',gridres_dlat,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'jDirectionIncrementInDegrees ' // &
           'in read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   ! EMK...Added dlon
   call grib_get(igrib, 'iDirectionIncrementInDegrees', gridres_dlon, ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'iDirectionIncrementInDegrees ' // &
           'in read_mogrepsg'
      call grib_release(igrib, ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'dataDate',dataDate,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'dataDate in read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'dataTime',dataTime,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'dataTime in read_mogrepsg'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   ! Here we tentatively have a file we can use. Close it for now, and
   ! prepare to pull the appropriate variables.
   call grib_release(igrib,ierr)
   call grib_close_file(ftn)

   ifguess = iginfo(1)
   jfguess = iginfo(2)

   call fldbld_read_mogrepsg_pcp(n, findex, order, gribfile, ifguess, jfguess, &
           prectot, rc)
   call assign_processed_mogrepsgforc_pcp(n, m, order, 8, prectot)
      
#endif

end subroutine read_mogrepsg_pcp

subroutine fldbld_read_mogrepsg_pcp(n, findex, order, gribfile, ifguess, jfguess, prectot, rc)                            
 
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,        intent(in)    :: n
  integer, intent(in)           :: findex  ! Forcing index
  integer,        intent(in)    :: order
  character(len=*),  intent(in) :: gribfile
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  real,           intent(out)   :: prectot(LIS_rc%lnc(n),LIS_rc%lnr(n))   !Total precipitation [kg/m^2]
  integer,        intent(out)   :: rc
!
! !DESCRIPTION:
!
!     To read MOGREPS-G data in GRIB-2 format.
!
!EOP
  character*9                   :: cstat
  character*100                 :: message     ( 20 )
  character(len=7)              :: grib_msg
  character(len=7)              :: check_mogrepsg_message
  integer                       :: count_prectot
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: igrib
  integer                       :: ftn
  integer                       :: kk, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, surface_val, level_val
  real, allocatable             :: dum1d   ( : )

  real, allocatable             :: fg_prectot ( : , : )

  logical                       :: found_inq

  rc = 1 ! Initialize as "no error"

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
  ! writing error messages to stdout/stderr, which may lead to runtime
  ! problems.
  inquire(file=trim(gribfile),exist=found_inq)
  if (.not. found_inq) then
     write(LIS_logunit,*) '[WARN] Cannot find file '//trim(gribfile)
     rc = 0
     return
  end if

#if (defined USE_GRIBAPI)

  ! If a problem occurs here, we can just return immediately since no
  ! memory has been allocated yet.
  call grib_open_file(ftn,trim(gribfile),'r',ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] Failed to open - '//trim(gribfile)
     rc = 0
     return
  end if

  allocate ( fg_prectot (ifguess, jfguess) )

  allocate ( dum1d   (ifguess*jfguess) )

  ! Initalization
  fg_prectot = LIS_rc%udef 
  dum1d = LIS_rc%udef

  prectot = LIS_rc%udef

  ! From this point, we must deallocate memory before returning.
  ! Unfortunately this means using a GOTO statement if a problem is
  ! encountered, but such is life.
  count_prectot = 0

  call grib_count_in_file(ftn,nvars,ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] in grib_count_in_file in ' // &
          'fldbld_read_mogrepsg'
     goto 100
  end if

  ! Tentatively loop through every field in GRIB file looking for the variables
  ! we want. The code below will exit the loop early if a problem is found *or*
  ! once all the required variables are found and read in.
  do kk=1,nvars

     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] failed to read - '//trim(gribfile)
        goto 100
     end if

     call grib_get(igrib,'discipline',param_disc_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grib_get: parameterNumber in ' // &
             'fldbld_read_mogrepsg'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterCategory in ' // &
             'fldbld_read_mogrepsg'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'parameterNumber',param_num_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterNumber in ' // &
             'fldbld_read_mogrepsg'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'typeOfFirstFixedSurface',surface_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_mogrepsg'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'level',level_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_mogrepsg'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! We have enough information to determine what GRIB parameter this is.
     !grib_msg = check_mogrepsg_message_pcp(param_disc_val, &
     !     param_cat_val, param_num_val, surface_val, level_val)
     if ( param_disc_val == 0 .and.      &
          param_cat_val  == 1 .and.  &
          param_num_val  == 49 .and. &
          surface_val    == 1 ) then
          grib_msg = 'prectot' ! Total precipitation [kg/m2]
     else
          grib_msg = 'none'
     endif

     ! Skip this field if GRIB parameter is not required.
     if (grib_msg == 'none') then
        call grib_release(igrib,ierr)
        if (ierr .ne. 0) then
           write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
                'fldbld_read_mogrepsg'
           goto 100
        end if
        cycle ! Not a message we are interested in.
     end if

     call grib_get(igrib,'values',dum1d,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: values for '//grib_msg// ' in ' // &
             'fldbld_read_mogrepsg_pcp'
        call grib_release(igrib,ierr)
        goto 100
     end if

     select case (grib_msg)
     case('prectot')  ! accumulated total precipitation
        fg_prectot = reshape(dum1d, (/ifguess,jfguess/))
        count_prectot = count_prectot + 1

     case default ! Internal error, we shouldn't be here
        write(LIS_logunit,*)'[WARN] Unknown grib_message ',grib_msg
        write(LIS_logunit,*)'Aborting...'
        flush(LIS_logunit)
        write(cstat,'(i9)',iostat=istat1) ierr
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  fldbld_read_mogrepsg.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(gribfile)
        if( istat1 .eq. 0 )then
           message(5) = '  Status = ' // trim(cstat)
        endif
        call LIS_abort( message)
     end select

     ! Finished with this field
     call grib_release(igrib,ierr)
     if (ierr .ne. 0) then
        write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
             'fldbld_read_mogrepsg'
        goto 100
     end if

     ! Jump out of loop early if we have everything
     if (count_prectot .eq. 1) then
        exit
     end if

  enddo ! Loop through all GRIB file fields

  ! Interpolate the fields to the LIS grid
  ! prectot
  call interp_mogrepsg_pcp(n, findex, ifguess, jfguess, .true., fg_prectot, prectot)

  ! At this point, we have everything.  Close the file and return.
  call grib_close_file(ftn)
  rc = 1
  return

  ! Jump down here to clean up memory before returning after finding a
  ! problem.
  100 continue
  call grib_close_file(ftn)

  deallocate ( dum1d      )
  deallocate ( fg_prectot )
  rc = 0
#endif

end subroutine fldbld_read_mogrepsg_pcp

!function check_mogrepsg_message_pcp(param_disc_val, &
!     param_cat_val, param_num_val, surface_val, level_val)
!! !USES:
!! none
!
!   implicit none
!! !ARGUMENTS:
!   integer, intent(in) :: param_disc_val, param_cat_val, &
!        param_num_val, surface_val, level_val
!   character(len=7)    :: check_mogrepsg_message_pcp
!!EOP
!
!   if ( param_disc_val == 0 .and.      &
!            param_cat_val  == 1 .and.  &
!            param_num_val  == 49 .and. &
!            surface_val    == 1 ) then
!      check_mogrepsg_message_pcp = 'prectot' ! Total precipitation [kg/m2]
!   else
!      check_mogrepsg_message_pcp = 'none'
!   endif
!end function check_mogrepsg_message_pcp

subroutine interp_mogrepsg_pcp(n, findex, ifguess, jfguess, pcp_flag, input, output)

! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use mogrepsg_forcingMod, only : mogrepsg_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)   :: findex
  integer, intent(in)  :: ifguess
  integer, intent(in)  :: jfguess
  logical, intent(in)  :: pcp_flag
  real,    intent(in)  :: input  ( ifguess,jfguess )
  real,    intent(out) :: output ( LIS_rc%lnc(n),LIS_rc%lnr(n) )
!
! !DESCRIPTION:
!
! This routine interpolates the MOGREPS-G data to the LIS grid.

!EOP

  integer   :: mi, mo
  integer   :: k
  integer   :: i,j
  integer   :: iret
  integer   :: midway
  character(len=50) :: method
  real, allocatable, dimension(:,:)    :: var
  logical*1, allocatable, dimension(:) :: lb
  logical*1, allocatable, dimension(:) :: lo

  allocate(var(ifguess,jfguess))
  allocate(lb(ifguess*jfguess))
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  mi = ifguess * jfguess
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  lb = .true.
  lo = .true.
  output = LIS_rc%udef

  ! translate from (0,360) to (-180,180).
  midway = ifguess/2
  do j = 1, jfguess
     do i = 1, ifguess
        if ( (i+midway) < ifguess ) then
           var(i,j) = input((i+midway),j)
        else
           var(i,j) = input((i-midway+1),j)
        endif
     enddo
  enddo

  ! Interpolate to LIS grid
  select case( LIS_rc%met_interp(findex) )

    case( "bilinear" )
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                     &
          var,lo,output,mi,mo,                                         &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          mogrepsg_struc(n)%w111,mogrepsg_struc(n)%w121,                   &
          mogrepsg_struc(n)%w211,mogrepsg_struc(n)%w221,                   &
          mogrepsg_struc(n)%n111,mogrepsg_struc(n)%n121,                   &
          mogrepsg_struc(n)%n211,mogrepsg_struc(n)%n221,LIS_rc%udef,iret)

    case( "budget-bilinear" )
     if (pcp_flag) then
        call conserv_interp(LIS_rc%gridDesc(n,:),lb,                   &
             var,lo,output,mi,mo,                                      &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
             mogrepsg_struc(n)%w112,mogrepsg_struc(n)%w122,                &
             mogrepsg_struc(n)%w212,mogrepsg_struc(n)%w222,                &
             mogrepsg_struc(n)%n112,mogrepsg_struc(n)%n122,                &
             mogrepsg_struc(n)%n212,mogrepsg_struc(n)%n222,LIS_rc%udef,iret)
     else
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                  &
             var,lo,output,mi,mo,                                      &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
             mogrepsg_struc(n)%w111,mogrepsg_struc(n)%w121,                &
             mogrepsg_struc(n)%w211,mogrepsg_struc(n)%w221,                &
             mogrepsg_struc(n)%n111,mogrepsg_struc(n)%n121,                &
             mogrepsg_struc(n)%n211,mogrepsg_struc(n)%n221,LIS_rc%udef,iret)
     endif

     case( "neighbor" )
        call neighbor_interp(LIS_rc%gridDesc(n,:),lb,                     &
             var,lo,output,mi,mo,                                         &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                        &
             mogrepsg_struc(n)%n113,LIS_rc%udef,iret)

     case DEFAULT
        write(LIS_logunit,*) 'ERR: Unexpected interpolation method'
        write(LIS_logunit,*) '     in interp_mogrepsg_first_guess'
        write(LIS_logunit,*) '     ', trim(LIS_rc%met_interp(findex))
        call LIS_endrun()  
  end select

  deallocate(var)
  deallocate(lb)
  deallocate(lo)

end subroutine interp_mogrepsg_pcp

!BOP
!
! !ROUTINE: assign_processed_mogrepsgforc
! \label{assign_processed_mogrepsgforc}
!
! !INTERFACE:
subroutine assign_processed_mogrepsgforc_pcp(n,m,order,var_index,mogrepsgforc)
! !USES:
  use LIS_coreMod
  use mogrepsg_forcingMod, only : mogrepsg_struc
!
! !DESCRIPTION:
!  This routine assigns the interpolated MOGREPS-G forcing data
!  to the module data structures to be used later for
!  time interpolation
!
!EOP
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: m 
  integer, intent(in) :: order
  integer, intent(in) :: var_index
  real,    intent(in) :: mogrepsgforc(LIS_rc%lnc(n),LIS_rc%lnr(n))

  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           if(order.eq.1) then
              mogrepsg_struc(n)%metdata1(var_index,m,&
                   LIS_domain(n)%gindex(c,r)) = &
                   mogrepsgforc(c,r)
           elseif(order.eq.2) then
              mogrepsg_struc(n)%metdata2(var_index,m,&
                   LIS_domain(n)%gindex(c,r)) = &
                   mogrepsgforc(c,r)
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_mogrepsgforc_pcp
