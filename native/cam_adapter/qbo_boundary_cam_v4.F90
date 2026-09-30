module qbo_boundary_cam_v4
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver
  use physics_types, only: physics_state
  use tropopause, only: tropopause_find, TROP_ALG_TWMO, TROP_ALG_CLIMATE, &
       TROP_ALG_CPP, TROP_ALG_NONE
  use qbo_boundary_v4, only: qbo_column_weight
  implicit none
  private
  public :: qbo_boundary_weights
contains
  subroutine qbo_boundary_weights(state, weight, pressure_pa, found)
    type(physics_state), intent(in) :: state
    real(r8), intent(out) :: weight(pcols,pver), pressure_pa(pcols,3), found(pcols,3)
    real(r8) :: column_pressure_hpa(3), column_found(3)
    integer :: levels(pcols), i, k, ncol

    ncol = state%ncol
    weight = 0._r8
    found = 0._r8
    call tropopause_find(state, levels, tropP=pressure_pa(:,1), &
         primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)
    found(:ncol,1) = merge(1._r8, 0._r8, levels(:ncol) > 0 .and. levels(:ncol) <= pver)
    call tropopause_find(state, levels, tropP=pressure_pa(:,2), &
         primary=TROP_ALG_TWMO, backup=TROP_ALG_NONE)
    found(:ncol,2) = merge(1._r8, 0._r8, levels(:ncol) > 0 .and. levels(:ncol) <= pver)
    call tropopause_find(state, levels, tropP=pressure_pa(:,3), &
         primary=TROP_ALG_CPP, backup=TROP_ALG_NONE)
    found(:ncol,3) = merge(1._r8, 0._r8, levels(:ncol) > 0 .and. levels(:ncol) <= pver)
    pressure_pa(ncol+1:,:) = 0._r8

    do i = 1, ncol
      column_pressure_hpa = pressure_pa(i,:) / 100._r8
      column_found = found(i,:)
      do k = 1, pver
        weight(i,k) = qbo_column_weight(state%pint(i,k+1) / 100._r8, &
             column_pressure_hpa, column_found)
      end do
    end do
  end subroutine qbo_boundary_weights
end module qbo_boundary_cam_v4
