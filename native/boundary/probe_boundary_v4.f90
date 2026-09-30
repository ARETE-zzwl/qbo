program probe_boundary_v4
  use, intrinsic :: iso_fortran_env, only: real64
  use qbo_boundary_v4, only: qbo_column_weight
  implicit none
  integer :: input, n, i, status
  real(real64) :: bottom, trop(3), found(3), weight
  character(len=1024) :: filename

  call get_command_argument(1, filename)
  open(newunit=input, file=trim(filename), status='old', action='read')
  read(input, *) n
  if (n < 1 .or. n > 20000) error stop 'Unexpected vector count'
  do i = 1, n
    read(input, *, iostat=status) bottom, trop, found
    if (status /= 0) error stop 'Invalid or incomplete input vector'
    weight = qbo_column_weight(bottom, trop, found)
    write(*, '(I8,1X,ES26.17E3)') i, weight
  end do
  close(input)
end program probe_boundary_v4
