program probe_cam_adapter_v4
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use physics_types, only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc, zonal_wind
  use tropopause, only: trop_calls
  use cam_history
  use qbo_test, only: qbo_relax, qbo_use_forcing, tauz, u_tstep, ktop, kbot
  implicit none
  type(physics_state) :: state, original
  type(physics_ptend) :: tendency
  type(physics_buffer_desc), target :: buffer_storage(1)
  type(physics_buffer_desc), pointer :: buffer(:)
  real(r8), parameter :: base_rate=1._r8/864000._r8
  real(r8) :: interfaces(pver+1), expected(pcols,pver)
  integer :: i, k, assertions=0

  buffer => buffer_storage
  interfaces = [100._r8, 1000._r8, 5000._r8, &
       (100._r8*exp(-.1_r8)-.001_r8)*100._r8, 10000._r8, 12000._r8, 20000._r8]
  do i=1,pcols
    state%pint(i,:) = interfaces
  end do
  state%trop_pa = 10000._r8
  state%trop_level = 4
  state%trop_level(2,2) = -1
  tauz = 1._r8
  tauz(1) = 2._r8
  tauz(pver) = 2._r8
  u_tstep = [-10._r8, 20._r8, 30._r8, 40._r8, 50._r8, 60._r8]
  zonal_wind = 5._r8
  original = state
  call qbo_relax(state, buffer, tendency)
  call check(trop_calls == 3, 'synchronous_three_definitions')
  call check(history_calls == 10, 'diagnostic_write_count')
  call check(abs(saved_mask(1,3)-.5_r8) <= 5e-14_r8, 'pa_hpa_half_transition')
  call check(all(saved_mask(1,4:) == 0._r8), 'lower_interface_and_margin')
  call check(all(saved_mask(2,:) == 0._r8), 'missing_primary_blocks_fallback')
  call check(all(saved_mask(3,:) == 0._r8), 'original_latitude_support')
  call check(all(saved_mask(4,:) == 0._r8) .and. all(saved_found(4,:) == 0._r8), 'padding_zero')
  expected = 0._r8
  expected(1,1:3) = [.5_r8, 1._r8, .5_r8] * base_rate
  call check(maxval(abs(saved_rate-expected)) < 1e-18_r8, 'applied_rate')
  do k=1,pver
    expected(:,k) = saved_rate(:,k) * (u_tstep(k)-zonal_wind(:,k))
  end do
  call check(maxval(abs(tendency%u-expected)) < 1e-18_r8, 'tendency_uses_same_rate')
  do k=1,pver
    expected(:,k) = saved_rate(:,k) * u_tstep(k) / base_rate
  end do
  call check(maxval(abs(saved_u0-expected)) < 1e-12_r8, 'u0_uses_same_weight')
  call check(all(saved_tend == tendency%u), 'saved_tendency_matches_applied')
  call check(saved_found(2,2) == 0._r8 .and. saved_pressure(1,1) == 10000._r8, 'pressure_and_found_units')
  call check(all(state%pint == original%pint) .and. all(state%trop_pa == original%trop_pa), 'state_unchanged')

  state%trop_pa = 30000._r8
  state%trop_level = 4
  call qbo_relax(state, buffer, tendency)
  call check(trop_calls == 6 .and. saved_mask(1,3) == 1._r8, 'fresh_state_not_cached')
  call check(all(saved_mask(1,5:6) == 0._r8), 'target_skip_at_or_above_50')

  u_tstep = 10._r8
  call qbo_relax(state, buffer, tendency)
  call check(abs(saved_rate(1,1)-base_rate/2) < 1e-18_r8 .and. &
       abs(saved_rate(1,6)-base_rate/2) < 1e-18_r8, 'both_half_strength_buffers')
  call check(abs(saved_rate(2,2)-base_rate*exp(-.1_r8**2/(2*.174532925_r8**2))) < 1e-18_r8, &
       'original_gaussian_latitude_rate')
  call check(all(ieee_is_finite(tendency%u)) .and. maxval(saved_rate) <= base_rate, 'finite_no_amplification')

  qbo_use_forcing = .false.
  trop_calls = 0
  history_calls = 0
  call qbo_relax(state, buffer, tendency)
  call check(trop_calls == 0 .and. history_calls == 0 .and. all(tendency%u == 0._r8), 'disabled_forcing')
  qbo_use_forcing = .true.
  state%trop_level(1,3) = -1
  call qbo_relax(state, buffer, tendency)
  call check(all(tendency%u(1,:) == 0._r8), 'missing_coldpoint')
  state%trop_level(1,3) = 4
  state%pint(1,3) = ieee_value(0._r8, ieee_quiet_nan)
  state%trop_pa(2,1) = 0._r8
  call qbo_relax(state, buffer, tendency)
  call check(tendency%u(1,2) == 0._r8 .and. saved_mask(1,2) == 0._r8, 'invalid_interface_fails_closed')
  call check(all(tendency%u(2,:) == 0._r8), 'invalid_tropopause_fails_closed')
  state%pint(1,:) = interfaces
  state%trop_pa = 30000._r8
  ktop = 3
  kbot = 4
  tauz = 1._r8
  tauz(2) = 2._r8
  tauz(5) = 2._r8
  call qbo_relax(state, buffer, tendency)
  call check(all(saved_mask(:,1) == 0._r8) .and. all(saved_mask(:,6) == 0._r8), 'outside_pressure_support')
  write(*,'(A,I0)') 'CAM_ADAPTER_INSTRUMENTED_PASS assertions=', assertions
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(A,A)') 'ASSERT_FAIL ', label
      error stop 1
    end if
    assertions = assertions + 1
    write(*,'(A,A)') 'ASSERT_PASS ', label
  end subroutine check
end program probe_cam_adapter_v4
