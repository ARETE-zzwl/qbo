! Instrumented CAM interfaces only; not the CAM tropopause algorithms.
module shr_kind_mod
  use, intrinsic :: iso_fortran_env, only: real64
  integer, parameter :: shr_kind_r8=real64
end module shr_kind_mod

module ppgrid
  integer, parameter :: pcols=4, pver=6
end module ppgrid

module physics_types
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  implicit none
  type physics_state
    integer :: ncol=3, psetcols=pcols, lchnk=1
    real(r8) :: pint(pcols,pver+1), trop_pa(pcols,3)
    integer :: trop_level(pcols,3)
  end type
  type physics_ptend
    real(r8) :: u(pcols,pver)
  end type
contains
  subroutine physics_ptend_init(ptend, psetcols, name, lu)
    type(physics_ptend), intent(out) :: ptend
    integer, intent(in) :: psetcols
    character(len=*), intent(in) :: name
    logical, intent(in) :: lu
    if (psetcols /= pcols .or. name /= 'qbo' .or. .not. lu) error stop 'Unexpected tendency initialization'
    ptend%u = 0._r8
  end subroutine
end module physics_types

module physics_buffer
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  implicit none
  type physics_buffer_desc
    integer :: unused=0
  end type
  real(r8), target :: zonal_wind(pcols,pver)=0._r8
contains
  subroutine pbuf_get_field(pbuf, index, field)
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer, intent(in) :: index
    real(r8), pointer :: field(:,:)
    if (.not. associated(pbuf) .or. index /= 1) error stop 'Unexpected buffer access'
    field => zonal_wind
  end subroutine
end module physics_buffer

module phys_grid
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  implicit none
  real(r8) :: latitudes(pcols) = [0._r8, .1_r8, .6_r8, 0._r8]
contains
  function get_rlat_p(lchnk, i) result(latitude)
    integer, intent(in) :: lchnk, i
    real(r8) :: latitude
    if (lchnk /= 1) error stop 'Unexpected chunk'
    latitude = latitudes(i)
  end function
end module phys_grid

module tropopause
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use physics_types, only: physics_state
  implicit none
  integer, parameter :: TROP_ALG_NONE=1, TROP_ALG_CLIMATE=3, TROP_ALG_TWMO=5, TROP_ALG_CPP=8
  integer :: trop_calls=0
contains
  subroutine tropopause_find(pstate, tropLev, tropP, primary, backup)
    type(physics_state), intent(in) :: pstate
    integer, intent(out) :: tropLev(pcols)
    real(r8), intent(out) :: tropP(pcols)
    integer, intent(in) :: primary, backup
    integer :: definition
    if (primary == TROP_ALG_TWMO .and. backup == TROP_ALG_CLIMATE) then
      definition=1
    else if (primary == TROP_ALG_TWMO .and. backup == TROP_ALG_NONE) then
      definition=2
    else if (primary == TROP_ALG_CPP .and. backup == TROP_ALG_NONE) then
      definition=3
    else
      error stop 'Unexpected tropopause algorithm or history-writing path'
    end if
    trop_calls = trop_calls + 1
    tropLev = pstate%trop_level(:,definition)
    tropP = pstate%trop_pa(:,definition)
  end subroutine
end module tropopause

module cam_history
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  implicit none
  real(r8) :: saved_tend(pcols,pver), saved_u0(pcols,pver)
  real(r8) :: saved_mask(pcols,pver), saved_rate(pcols,pver)
  real(r8) :: saved_pressure(pcols,3), saved_found(pcols,3)
  integer :: history_calls=0
  interface outfld
    module procedure outfld_2d, outfld_1d
  end interface
contains
  subroutine outfld_2d(name, field, ldim, lchnk)
    character(len=*), intent(in) :: name
    integer, intent(in) :: ldim, lchnk
    real(r8), intent(in) :: field(ldim,pver)
    if (ldim /= pcols .or. lchnk /= 1) error stop 'Invalid history shape/chunk'
    history_calls = history_calls + 1
    select case(name)
    case('QBOTEND')
      saved_tend = field
    case('QBO_U0')
      saved_u0 = field
    case('QBOMASK')
      saved_mask = field
    case('QBO_RATE')
      saved_rate = field
    case default
      error stop 'Unexpected 2D history write'
    end select
  end subroutine
  subroutine outfld_1d(name, field, ldim, lchnk)
    character(len=*), intent(in) :: name
    integer, intent(in) :: ldim, lchnk
    real(r8), intent(in) :: field(ldim)
    if (ldim /= pcols .or. lchnk /= 1) error stop 'Invalid history shape/chunk'
    history_calls = history_calls + 1
    select case(name)
    case('QBO_TROP')
      saved_pressure(:,1) = field
    case('QBO_TRPP')
      saved_pressure(:,2) = field
    case('QBO_TRPF')
      saved_pressure(:,3) = field
    case('QBO_FD')
      saved_found(:,1) = field
    case('QBO_FDP')
      saved_found(:,2) = field
    case('QBO_FDF')
      saved_found(:,3) = field
    case default
      error stop 'Unexpected 1D history write'
    end select
  end subroutine
end module cam_history
