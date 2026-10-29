! Minimal data containers for compiling the actual SNIP updater in isolation.
! No physics is stubbed: use the unchanged NoahMP50_snowphys_updateMod.F90.
module LisNoahmpParamType
  implicit none
  type LisNoahmpParam_type
     real :: BEXP(4)=4.0, PSISAT(4)=0.1, SMCMAX(4)=0.4
  end type
end module

module NoahMP50_lsmMod
  use LisNoahmpParamType
  implicit none
  type tile_type
     integer :: isnow=0, soiltype(4)=1
     real :: snowh=0.0, sneqv=0.0, sfctmp=260.0
     real :: snowice(3)=0.0, snowliq(3)=0.0, tsno(3)=0.0, zss(7)=0.0
     real :: tslb(4)=270.0, sh2o(4)=0.2, smc(4)=0.3
     type(LisNoahmpParam_type) :: param
  end type
  type nest_type
     integer :: nsnow=3, nsoil=4
     real :: ts=1800.0
     real :: sldpth(4)=(/0.1,0.3,0.6,1.0/)
     type(tile_type) :: noahmp50(1)
  end type
  type(nest_type) :: Noahmp50_struc(1)
end module
