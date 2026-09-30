program read_qbo_native
  use iso_fortran_env, only: real64
  use netcdf
  implicit none
  integer :: ncid, varid, n, nt, k, t
  integer, allocatable :: dates(:), secs(:)
  real(real64) :: p0
  real(real64), allocatable :: u(:,:), p(:), a(:), b(:), ai(:), bi(:)
  character(len=16) :: mode
  character(len=1024) :: filename

  call get_command_argument(1, mode)
  call get_command_argument(2, filename)
  call check(nf90_open(trim(filename), nf90_nowrite, ncid))
  select case (trim(mode))
  case ("qbo")
    n = dim_size("level")
    nt = dim_size("time")
    allocate(u(n,nt), p(n), dates(nt), secs(nt))
    call check(nf90_inq_varid(ncid, "qbo", varid))
    ! This is the same allocation and untransposed call as pinned WACCM.
    call check(nf90_get_var(ncid, varid, u))
    call read_vector("level", p)
    call check(nf90_inq_varid(ncid, "date", varid))
    call check(nf90_get_var(ncid, varid, dates))
    call check(nf90_inq_varid(ncid, "secs", varid))
    call check(nf90_get_var(ncid, varid, secs))
    print '(A)', "# time_index level_index date secs pressure_hpa wind_mps"
    do t = 1, nt
      do k = 1, n
        print '(4(I8,1X),2(ES25.17E3,1X))', t, k, dates(t), secs(t), p(k), u(k,t)
      end do
    end do
  case ("grid")
    n = dim_size("lev")
    allocate(a(n), b(n), ai(n+1), bi(n+1))
    call read_vector("hyam", a)
    call read_vector("hybm", b)
    call read_vector("hyai", ai)
    call read_vector("hybi", bi)
    call check(nf90_inq_varid(ncid, "P0", varid))
    call check(nf90_get_var(ncid, varid, p0))
    print '(A,ES25.17E3)', "# P0_Pa = ", p0
    print '(A)', "# index midpoint_hpa top_interface_hpa bottom_interface_hpa hybm"
    do k = 1, n
      print '(I8,1X,4(ES25.17E3,1X))', k, (a(k)+b(k))*p0/100, &
        (ai(k)+bi(k))*p0/100, (ai(k+1)+bi(k+1))*p0/100, b(k)
    end do
  case default
    error stop "mode must be qbo or grid"
  end select
  call check(nf90_close(ncid))

contains
  subroutine check(code)
    integer, intent(in) :: code
    if (code /= nf90_noerr) then
      print *, "NETCDF_STATUS", code, trim(nf90_strerror(code))
      stop 2
    end if
  end subroutine

  integer function dim_size(name) result(length)
    character(len=*), intent(in) :: name
    integer :: dimid
    call check(nf90_inq_dimid(ncid, name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, len=length))
  end function

  subroutine read_vector(name, values)
    character(len=*), intent(in) :: name
    real(real64), intent(out) :: values(:)
    call check(nf90_inq_varid(ncid, name, varid))
    call check(nf90_get_var(ncid, varid, values))
  end subroutine
end program
