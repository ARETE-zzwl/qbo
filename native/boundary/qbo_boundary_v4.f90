module qbo_boundary_v4
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: qbo_column_weight
  real(real64), parameter :: margin_hpa = 0.001_real64
  real(real64), parameter :: log_width = 0.2_real64
contains
  pure function qbo_column_weight(bottom_hpa, trop_hpa, found) result(weight)
    real(real64), intent(in) :: bottom_hpa, trop_hpa(3), found(3)
    real(real64) :: weight, safe, x

    weight = 0.0_real64
    ! Separate finite checks avoid invalid comparisons with NaN under FPE traps.
    if (.not. ieee_is_finite(bottom_hpa)) return
    if (.not. all(ieee_is_finite(trop_hpa))) return
    if (.not. all(ieee_is_finite(found))) return
    if (bottom_hpa <= 0.0_real64) return
    if (any(trop_hpa <= 0.0_real64)) return
    if (any(found <= 0.5_real64)) return
    safe = minval(trop_hpa)
    if (bottom_hpa + margin_hpa >= safe) return
    x = (log(safe) - log(bottom_hpa + margin_hpa)) / log_width
    x = min(1.0_real64, max(0.0_real64, x))
    weight = x * x * (3.0_real64 - 2.0_real64 * x)
  end function qbo_column_weight
end module qbo_boundary_v4
