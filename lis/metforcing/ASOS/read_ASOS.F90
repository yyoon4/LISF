!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_ASOS.F90
! \label{read_ASOS}
!
! !REVISION HISTORY: 
!  19 Apr 2013: Jonathan Case; Corrected some initializations and comparison
!               against undefined value, incorporated here by Soni Yatheendradas
!  21 Jun 2013: Soni Yatheendradas; Based on new TRMM3B42V6 code, changes from
!               earlier code to avoid (a) alternate file skip, (b) jump to
!               previous day TRMM, and (c) absence of rain rate weighting
!  06 Sep 2017: Minwook Kim; Based on TRMM3B42V7 code, changes from
!               earlier code
!
! !INTERFACE:
subroutine read_ASOS (n, kk, fname, findex, order, ferror_ASOS)
! !USES:
  use LIS_coreMod,           only : LIS_rc, LIS_domain
  use LIS_logMod,            only : LIS_logunit, LIS_getNextUnitNumber, &
                                    LIS_releaseUnitNumber
  use LIS_metforcingMod,     only : LIS_forc
  use ASOS_forcingMod, only : ASOS_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: kk
  character(len=80)   :: fname
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_ASOS
  !integer             :: filehr ! SY

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  ASOS data and interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the 12 hour ASOS file
!  \item[ferror\_ASOS]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[interp\_ASOS](\ref{interp_ASOS}) \newline
!    spatially interpolates the ASOS data
!  \end{description}
!EOP

  integer :: index1, nASOS

  !==== Local Variables=======================

  integer :: ios
  integer :: i,j,xd,yd
  parameter(xd=33,yd=45)                            ! Dimension of original ASOS data

  real :: precip(xd,yd)                                ! Original real precipitation array
  logical*1,allocatable  :: lb(:)
  real, allocatable :: precip_regrid(:,:)                  ! Interpolated precipitation array
  integer :: ftn

  !=== End Variable Definition =======================

  !------------------------------------------------------------------------
  ! Fill necessary arrays to assure not using old ASOS data
  !------------------------------------------------------------------------
  ! J.Case (4/22/2013) -- Make consistent with Stg4/NMQ routines
  if(order.eq.1) then
     ASOS_struc(n)%metdata1 = LIS_rc%udef ! J.Case
  elseif(order.eq.2) then 
     ASOS_struc(n)%metdata2 = LIS_rc%udef ! J.Case
  endif
  allocate (precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!  precip = -1.0
!  if(order.eq.1) then
!     LIS_forc(n,findex)%metdata1 = -1.0
!  elseif(order.eq.2) then 
!     LIS_forc(n,findex)%metdata2 = -1.0
!  endif
  precip_regrid = -1.0 ! J.Case
  !------------------------------------------------------------------------
  ! Find ASOS precip data, read it in and assign to forcing precip array.
  ! Must reverse grid in latitude dimension to be consistent with LDAS grid
  !------------------------------------------------------------------------
  ftn = LIS_getNextUnitNumber()
  open(unit=ftn,file=fname, form='unformatted')

  if (ios .eq. 0) then
     read (ftn) precip
     Do j=1, xd
        Do i=1, yd
           if (precip(j, i) .LT. 0.0 ) precip(j, i) = 0.0    ! reset to 0 for weird values
        End Do
     End Do

! J.Case (4/19/2013) -- Test print out of raw precip array
! write (99,*) precip

     !------------------------------------------------------------------------
     ! Interpolating to desired LIS_domain and resolution
     ! Global precip datasets not used currently to force NLDAS
     !------------------------------------------------------------------------
     !print*, "Writing un-interpolated ASOS precipitation out "
     !open(71, file="ASOS-ungrid.1gd4r", access="direct", &
     !    recl=xd*yd*4, form="unformatted")
     ! write(71, rec=1) precip
     !close(71)

     nASOS = ASOS_struc(n)%ncold*ASOS_struc(n)%nrold
     allocate(lb(nASOS))
     lb = .true.
     call interp_ASOS(n, nASOS, precip, lb, LIS_rc%gridDesc, &
          !LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid) ! SY
          LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid, findex) ! SY
     deallocate (lb) 

     !print*, "Writing interpolated ASOS precipitation out "
     !open(73, file="ASOS-regrid.1gd4r", access="direct", &
     !    recl=LIS_rc%d%lnr*LIS_rc%d%lnc*4, form="unformatted")
     ! write(73, rec=1) precip_regrid
     !close(73)
     !print*, "Writing interpolated ASOS precipitation out finished"

! J.Case (4/19/2013) -- Test print out of the regridded precip (on LIS grid).
! write (98,*) precip_regrid

     do j = 1,LIS_rc%lnr(n)
        do i = 1,LIS_rc%lnc(n)
           if (precip_regrid(i,j) .ne. -1.0) then
              index1 = LIS_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then
                 if(order.eq.1) then 
                    ASOS_struc(n)%metdata1(kk,1,index1) = precip_regrid(i,j)   !here is mm/h
                 elseif(order.eq.2) then 
                    ASOS_struc(n)%metdata2(kk,1,index1) = precip_regrid(i,j)   !here is mm/h
                 endif
              endif
           endif
        enddo
     enddo

! J.Case (4/19/2013) -- Test print out of the suppdata precip (on LIS grid).
! write (97,*) LIS_forc(n,findex)%metdata1(1,:)

     ferror_ASOS = 1
     write(LIS_logunit,*) "Obtained ASOS precipitation data ", fname
  else
     write(LIS_logunit,*) "Missing ASOS precipitation data ", fname
     ferror_ASOS = 0
  endif
  call LIS_releaseUnitNumber(ftn)

  deallocate (precip_regrid)

end subroutine read_ASOS
