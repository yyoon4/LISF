! Minimal framework doubles for the standalone SNIP interface checks.
! NoahMP types, parameters and snow-layer physics are compiled from source.
module ESMF
  implicit none
  type ESMF_Field
     real, pointer :: values(:) => null()
     real :: minimum = 0.0, maximum = 10000.0
  end type
  type ESMF_State
     type(ESMF_Field) :: swe, snowdepth, observation
  end type
contains
  subroutine ESMF_StateGet(state, name, field, itemCount, itemNameList, rc)
    type(ESMF_State), intent(in) :: state
    character(*), optional, intent(in) :: name
    type(ESMF_Field), optional, intent(out) :: field
    integer, optional, intent(out) :: itemCount
    character(*), optional, intent(out) :: itemNameList(:)
    integer, intent(out) :: rc
    rc = 0
    if (present(itemCount)) itemCount = 1
    if (present(itemNameList)) itemNameList = 'Observation01'
    if (present(field)) then
       select case (name)
       case ('SWE')
          field = state%swe
       case ('Snowdepth')
          field = state%snowdepth
       case ('Observation01')
          field = state%observation
       case default
          error stop 'Unknown test field'
       end select
    end if
  end subroutine
  subroutine ESMF_FieldGet(field, localDE, farrayPtr, rc)
    type(ESMF_Field), intent(in) :: field
    integer, intent(in) :: localDE
    real, pointer :: farrayPtr(:)
    integer, intent(out) :: rc
    farrayPtr => field%values
    rc = 0
  end subroutine
  subroutine ESMF_AttributeGet(field, name, value, rc)
    type(ESMF_Field), intent(in) :: field
    character(*), intent(in) :: name
    real, intent(out) :: value
    integer, intent(out) :: rc
    value = field%minimum
    if (name == 'Max Value') value = field%maximum
    rc = 0
  end subroutine
end module

module LIS_constantsMod
  implicit none
  integer, parameter :: LIS_CONST_PATH_LEN = 1024
  real, parameter :: LIS_CONST_TKFRZ = 273.15
end module

module LIS_coreMod
  implicit none
  type runtime_type
     integer :: nnest = 1, lsm_index = 1
     integer :: npatch(1,1) = 5, ngrid(1) = 1, nensem(1) = 5
     integer :: obs_ngrid(1) = 1, snowsrc(1) = 0
  end type
  type tile_type
     integer :: col = 1, row = 1, index = 1, tile_id = 1, vegt = 1
  end type
  type surface_type
     type(tile_type) :: tile(5)
  end type
  type domain_type
     integer :: gindex(1,1) = 1
  end type
  type(runtime_type) :: LIS_rc
  type(surface_type) :: LIS_surface(1,1)
  type(domain_type) :: LIS_domain(1)
end module

module LIS_logMod
  implicit none
contains
  subroutine LIS_verify(status, message)
    integer, intent(in) :: status
    character(*), optional, intent(in) :: message
    if (status /= 0) error stop 'Framework error'
  end subroutine
  subroutine LIS_endrun()
    error stop 'Unexpected LIS_endrun'
  end subroutine
end module

module LIS_snowMod
  implicit none
  type snow_type
     real :: snowdepth(1) = 0.0, sneqv(1) = 0.0
  end type
  type(snow_type) :: LIS_snow_struc(1)
end module

module LIS_DAobservationsMod
  implicit none
contains
  subroutine LIS_convertPatchSpaceToObsEnsSpace(n, k, lsm, values, result)
    integer, intent(in) :: n, k, lsm
    real, intent(in) :: values(:)
    real, intent(out) :: result(:,:)
    result(1,:) = values
  end subroutine
  subroutine LIS_convertPatchSpaceToObsSpace(n, k, lsm, values, result)
    integer, intent(in) :: n, k, lsm
    real, intent(in) :: values(:)
    real, intent(out) :: result(:)
    result = sum(values) / size(values)
  end subroutine
end module

module LIS_lsmMod
  implicit none
contains
  subroutine LIS_lsm_DAmapTileSpaceToObsSpace(n, k, t, st_id, en_id)
    integer, intent(in) :: n, k, t
    integer, intent(out) :: st_id, en_id
    st_id = 1
    en_id = 1
  end subroutine
end module
