program probe_zonal_v5_mpi
  use, intrinsic :: iso_fortran_env, only: real64, output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi
  implicit none
  integer :: ierr, rank, ranks, unit, ios, nrings, nlon, i, j, owner
  integer :: local_count, global_count
  real(real64), allocatable :: weights(:,:), local_minimum(:), global_minimum(:)
  character(len=1024) :: input_path
  character(len=32) :: partition

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, ranks, ierr)
  if (command_argument_count() /= 2) call MPI_Abort(MPI_COMM_WORLD, 10, ierr)
  call get_command_argument(1, input_path)
  call get_command_argument(2, partition)
  open(newunit=unit, file=trim(input_path), status='old', action='read', iostat=ios)
  if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
  read(unit, *, iostat=ios) nrings, nlon
  if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 12, ierr)
  if (nrings <= 0 .or. nlon <= 0) call MPI_Abort(MPI_COMM_WORLD, 13, ierr)
  if (real(nrings,real64)*real(nlon,real64) > 2000000._real64) then
    call MPI_Abort(MPI_COMM_WORLD, 14, ierr)
  end if
  allocate(weights(nlon,nrings), local_minimum(nrings), global_minimum(nrings))
  do i = 1, nrings
    read(unit, *, iostat=ios) weights(:,i)
    if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 15, ierr)
  end do
  close(unit)
  if (.not. all(ieee_is_finite(weights))) call MPI_Abort(MPI_COMM_WORLD, 16, ierr)
  if (any(weights < 0._real64) .or. any(weights > 1._real64)) then
    call MPI_Abort(MPI_COMM_WORLD, 16, ierr)
  end if

  ! Empty owners contribute one, the neutral minimum for valid [0,1] weights.
  local_minimum = 1._real64
  local_count = 0
  do j = 1, nlon
    select case (trim(partition))
    case ('cyclic', 'incomplete')
      owner = mod(j-1, ranks)
      if (trim(partition) == 'incomplete' .and. j == nlon) cycle
    case ('block')
      owner = min(ranks-1, ((j-1)*ranks)/nlon)
    case ('single')
      owner = 0
    case default
      call MPI_Abort(MPI_COMM_WORLD, 18, ierr)
      owner = -1
    end select
    if (owner /= rank) cycle
    local_count = local_count + 1
    local_minimum = min(local_minimum, weights(j,:))
  end do
  call MPI_Allreduce(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (global_count /= nlon) then
    if (rank == 0) then
      write(output_unit,'(a,i0,a,i0)') 'PARTITION_COVERAGE_REJECTED got=', global_count, ' expected=', nlon
      flush(output_unit)
    end if
    call MPI_Finalize(ierr)
    if (rank == 0) stop 17
    stop
  end if
  call MPI_Allreduce(local_minimum, global_minimum, nrings, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  if (rank == 0) then
    do i = 1, nrings
      write(output_unit,'(i8,1x,es25.17e3)') i, global_minimum(i)
    end do
  end if
  deallocate(weights, local_minimum, global_minimum)
  call MPI_Finalize(ierr)
end program probe_zonal_v5_mpi
