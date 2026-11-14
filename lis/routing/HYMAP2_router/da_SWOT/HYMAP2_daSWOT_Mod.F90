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
module HYMAP2_daSWOT_Mod
!BOP
!
! !MODULE: HYMAP2_daSWOT_Mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
!                          copied from HYMAP2_daWL_Mod
!
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP2_daSWOT_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: HYMAP2_daSWOT_struc
!EOP

 type, public :: daSWOT_dec

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal
     integer                :: nsites
     integer                :: useLocalUpd
     integer                :: localUpdDX
     character(len=32)      :: weightingMethod = "area-max"
     real                   :: distanceCutoffKm = 0.0
     real, allocatable      :: sites(:,:)
     real, allocatable      :: localWeight(:,:,:)

  end type daSWOT_dec

  type(daSWOT_dec), allocatable :: HYMAP2_daSWOT_struc(:)

contains
!BOP
!
! !ROUTINE: HYMAP2_daSWOT_init
! \label{HYMAP2_daSWOT_init}
!
! !INTERFACE:
  subroutine HYMAP2_daSWOT_init(k)
! !USES:
    use HYMAP2_routingMod
    use HYMAP2_initMod
    use LIS_constantsMod, only: LIS_CONST_PATH_LEN
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
!
! !DESCRIPTION:
!
!  This routine initializes the DA related data structures in HYMAP2
!  and the runtime DA configuration options
!
!EOP

    implicit none

    integer                :: k
    integer                :: n
    integer                :: nid
    integer                :: status
    integer                :: ncId, nrId, siteId
    integer                :: nc,nr,drainid,sid
    logical                :: local_upd_flag
    character(len=LIS_CONST_PATH_LEN) :: localWeightMap
    logical                :: file_exists

    if(.not.allocated(HYMAP2_daSWOT_struc)) then
       allocate(HYMAP2_daSWOT_struc(LIS_rc%nnest))

       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP2 use localization update in DA:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP2_daSWOT_struc(n)%useLocalUpd,rc=status)
          call LIS_verify(status,&
               "HYMAP2 use localization update in DA: not defined")
       enddo

       ! Optional labels: old configurations retain the exact area-map path.
       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP2 localization weighting method:",rc=status)
       if(status.eq.0) then
          do n=1,LIS_rc%nnest
             call ESMF_ConfigGetAttribute(LIS_config,&
                  HYMAP2_daSWOT_struc(n)%weightingMethod,rc=status)
             call LIS_verify(status,"Invalid HYMAP2 localization weighting method")
          enddo
       endif
       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP2 localization distance cutoff (km):",rc=status)
       if(status.eq.0) then
          do n=1,LIS_rc%nnest
             call ESMF_ConfigGetAttribute(LIS_config,&
                  HYMAP2_daSWOT_struc(n)%distanceCutoffKm,rc=status)
             call LIS_verify(status,"Invalid HYMAP2 localization distance cutoff")
          enddo
       endif
       do n=1,LIS_rc%nnest
          select case(trim(HYMAP2_daSWOT_struc(n)%weightingMethod))
          case("area-max")
          case("area-distance-max")
             if(HYMAP2_daSWOT_struc(n)%useLocalUpd.eq.1) then
                if(.not.ieee_is_finite(HYMAP2_daSWOT_struc(n)%distanceCutoffKm)) then
                   call HYMAP2_localization_error("Distance cutoff must be finite and positive")
                endif
                if(HYMAP2_daSWOT_struc(n)%distanceCutoffKm.le.0.0) then
                   call HYMAP2_localization_error("Distance cutoff must be finite and positive")
                endif
#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
                call HYMAP2_localization_error("Distance localization requires NetCDF support")
#endif
             endif
          case default
             call HYMAP2_localization_error("Unknown localization weighting method")
          end select
          if(HYMAP2_daSWOT_struc(n)%useLocalUpd.eq.1) then
             write(LIS_logunit,*) "[INFO] HYMAP2 localization weighting: nest=",n,&
                  " method=",trim(HYMAP2_daSWOT_struc(n)%weightingMethod),&
                  " cutoff(km)=",HYMAP2_daSWOT_struc(n)%distanceCutoffKm
          endif
       enddo

       local_upd_flag = .false.
       do n=1,LIS_rc%nnest
          if(HYMAP2_daSWOT_struc(n)%useLocalUpd.eq.1) then
             local_upd_flag = .true.
          endif
       enddo

       if(local_upd_flag) then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP2 localization update window size:",rc=status)
          do n=1, LIS_rc%nnest
             call ESMF_ConfigGetAttribute(LIS_config,&
                  HYMAP2_daSWOT_struc(n)%localupdDX,rc=status)
             call LIS_verify(status,&
                  "HYMAP2 localization update window size: not defined")
          enddo

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP2 localization weight map:",rc=status)
          do n=1, LIS_rc%nnest
             call ESMF_ConfigGetAttribute(LIS_config,&
                  localWeightMap,rc=status)
             call LIS_verify(status,&
                  "HYMAP2 localization weight map: not defined")
             ! Consume each nest's map token, but disabled nests need no map.
             if(HYMAP2_daSWOT_struc(n)%useLocalUpd.ne.1) cycle
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=localWeightMap, exist=file_exists)
             if(file_exists) then

                write(LIS_logunit,*)'[INFO] Reading localization map from '&
                     //trim(localWeightMap)
                status = nf90_open(path=localWeightMap,&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(status, &
                     'Error in nf90_open in HYMAP2_daSWOT_Mod')

                status = nf90_inq_dimid(nid,"east_west",ncId)
                call LIS_verify(status, &
                     'Error in nf90_inq_dimid in HYMAP2_daSWOT_Mod')

                status = nf90_inq_dimid(nid,"north_south",nrId)
                call LIS_verify(status, &
                     'Error in nf90_inq_dimid in HYMAP2_daSWOT_Mod')

                status = nf90_inq_dimid(nid,"count",siteId)
                call LIS_verify(status, &
                     'Error in nf90_inq_dimid in HYMAP2_daSWOT_Mod')

                status = nf90_inquire_dimension(nid,ncId, len=nc)
                call LIS_verify(status, &
                     'Error in nf90_inquire_dimension in HYMAP2_daSWOT_Mod')

                status = nf90_inquire_dimension(nid,nrId, len=nr)
                call LIS_verify(status, &
                     'Error in nf90_inquire_dimension in HYMAP2_daSWOT_Mod')

                status = nf90_inquire_dimension(nid,siteId, &
                     len=HYMAP2_daSWOT_struc(n)%nsites)
                call LIS_verify(status, &
                     'Error in nf90_inquire_dimension in HYMAP2_daSWOT_Mod')

                status = nf90_inq_varid(nid,'weight',drainid)
                call LIS_verify(status, &
                     'distance field not found in the localization weight file')

                status = nf90_inq_varid(nid,'sites',sid)
                call LIS_verify(status, &
                     'sites field not found in the localization weight file')
                allocate(HYMAP2_daSWOT_struc(n)%localWeight(&
                     LIS_rc%gnc(n), LIS_rc%gnr(n), &
                     HYMAP2_daSWOT_struc(n)%nsites))

                allocate(HYMAP2_daSWOT_struc(n)%sites( &
                     LIS_rc%lnc(n), LIS_rc%lnr(n)))

                status = nf90_get_var(nid,drainid,&
                     HYMAP2_daSWOT_struc(n)%localWeight)
                call LIS_verify(status, &
                     'Error in nf90_get_var in HYMAP2_daSWOT_Mod')

                status = nf90_get_var(nid,sid,&
                     HYMAP2_daSWOT_struc(n)%sites, &
                     start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
                     LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                     count = (/LIS_ewe_halo_ind(n,LIS_localPet+1) - &
                     LIS_ews_halo_ind(n,LIS_localPet+1)+1, &
                     LIS_nse_halo_ind(n,LIS_localPet+1) - &
                     LIS_nss_halo_ind(n,LIS_localPet+1)+1/))

                call LIS_verify(status, &
                     'Error in nf90_get_var in HYMAP2_daSWOT_Mod')

                if(trim(HYMAP2_daSWOT_struc(n)%weightingMethod).eq."area-distance-max") then
                   call HYMAP2_taper_localization_weights(n,nid,drainid,ncId,nrId,siteId)
                endif

                status = nf90_close(nid)
                call LIS_verify(status, &
                     'Error in nf90_close in HYMAP2_daSWOT_Mod')

             else
                write(LIS_logunit,*) '[ERR] localization map: ', &
                     trim(localWeightMap), ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
          enddo
       endif
    endif
  end subroutine HYMAP2_daSWOT_init

  subroutine HYMAP2_localization_error(message)
    use LIS_logMod, only: LIS_logunit,LIS_endrun
    implicit none
    character(len=*),intent(in) :: message
    write(LIS_logunit,*) "[ERR] HYMAP2 localization: ",trim(message)
    call LIS_endrun
  end subroutine HYMAP2_localization_error

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  subroutine HYMAP2_taper_localization_weights(n,nid,weightId,ncId,nrId,siteDimId)
    use netcdf
    use LIS_coreMod
    use LIS_logMod, only: LIS_verify
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: n,nid,weightId,ncId,nrId,siteDimId
    integer :: distanceId,status,ndims,weightDims(3),distanceDims(3)
    integer :: c,r,s,nx,ny,gx,gy,nc,nr
    character(len=32) :: units
    real(kind=8),allocatable :: distance(:,:)
    real(kind=8) :: missingDistance,cutoff,pi,taper,d
    real :: w

    ! Validate only the opt-in path; do not add NetCDF lookups in area-max.
    status=nf90_inquire_variable(nid,weightId,ndims=ndims)
    call LIS_verify(status,"Unable to inspect localization weight dimensions")
    if(ndims.ne.3) call HYMAP2_localization_error("Weight must have three dimensions")
    status=nf90_inquire_variable(nid,weightId,dimids=weightDims)
    call LIS_verify(status,"Unable to inspect localization weight axes")
    if(any(weightDims.ne.(/ncId,nrId,siteDimId/))) then
       call HYMAP2_localization_error("Weight axes must be east_west,north_south,count in Fortran")
    endif
    status=nf90_inquire_dimension(nid,ncId,len=nc)
    call LIS_verify(status,"Unable to inspect localization grid width")
    status=nf90_inquire_dimension(nid,nrId,len=nr)
    call LIS_verify(status,"Unable to inspect localization grid height")
    if(nc.ne.LIS_rc%gnc(n).or.nr.ne.LIS_rc%gnr(n)) then
       call HYMAP2_localization_error("Localization map must match the global model grid")
    endif
    status=nf90_inq_varid(nid,"network_distance_km",distanceId)
    call LIS_verify(status,"area-distance-max requires network_distance_km")
    status=nf90_inquire_variable(nid,distanceId,ndims=ndims)
    call LIS_verify(status,"Unable to inspect network_distance_km dimensions")
    if(ndims.ne.3) call HYMAP2_localization_error("network_distance_km must have three dimensions")
    status=nf90_inquire_variable(nid,distanceId,dimids=distanceDims)
    call LIS_verify(status,"Unable to inspect network_distance_km axes")
    if(any(distanceDims.ne.weightDims)) then
       call HYMAP2_localization_error("network_distance_km axes must match weight")
    endif
    units=""
    status=nf90_get_att(nid,distanceId,"units",units)
    call LIS_verify(status,"network_distance_km requires units=km")
    if(trim(adjustl(units)).ne."km") call HYMAP2_localization_error("network_distance_km units must be km")
    status=nf90_get_att(nid,distanceId,"_FillValue",missingDistance)
    call LIS_verify(status,"network_distance_km requires _FillValue=-9999")
    if(.not.ieee_is_finite(missingDistance)) call HYMAP2_localization_error("Invalid distance fill value")
    if(missingDistance.ne.-9999.0_8) call HYMAP2_localization_error("Distance fill value must be -9999")

    nx=size(HYMAP2_daSWOT_struc(n)%localWeight,1)
    ny=size(HYMAP2_daSWOT_struc(n)%localWeight,2)
    gx=1
    gy=1
    if(gx.lt.1.or.gy.lt.1.or.gx+nx-1.gt.nc.or.gy+ny-1.gt.nr) then
       call HYMAP2_localization_error("Distance hyperslab lies outside the map")
    endif
    allocate(distance(nx,ny))
    cutoff=dble(HYMAP2_daSWOT_struc(n)%distanceCutoffKm)
    pi=acos(-1.0_8)
    do s=1,HYMAP2_daSWOT_struc(n)%nsites
       ! Reuse one 2-D distance slice; never allocate another full 3-D map.
       status=nf90_get_var(nid,distanceId,distance,&
            start=(/gx,gy,s/),count=(/nx,ny,1/))
       call LIS_verify(status,"Unable to read network_distance_km slice")
       do r=1,ny
          do c=1,nx
             w=HYMAP2_daSWOT_struc(n)%localWeight(c,r,s)
             if(w.eq.-9999.0) cycle
             if(.not.ieee_is_finite(w)) call HYMAP2_localization_error("Nonfinite base localization weight")
             if(w.lt.0.0.or.w.gt.1.0) call HYMAP2_localization_error("Base localization weight outside [0,1]")
             d=distance(c,r)
             if(.not.ieee_is_finite(d)) then
                HYMAP2_daSWOT_struc(n)%localWeight(c,r,s)=-9999.0
             elseif(d.lt.0.0_8.or.d.ge.cutoff) then
                HYMAP2_daSWOT_struc(n)%localWeight(c,r,s)=-9999.0
             else
                taper=0.5_8*(1.0_8+cos(pi*d/cutoff))
                HYMAP2_daSWOT_struc(n)%localWeight(c,r,s)=real(dble(w)*taper)
             endif
          enddo
       enddo
    enddo
    deallocate(distance)
  end subroutine HYMAP2_taper_localization_weights
#endif
end module HYMAP2_daSWOT_Mod
