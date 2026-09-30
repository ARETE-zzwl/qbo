module qbo_contract
  use mpi
  implicit none
  private
  public :: contract_context, contract_begin, contract_observe, contract_collect
  integer, parameter :: capacity=8
  type :: contract_context
    private
    integer :: epoch=-1, consumed=-1, mode=0, levels=0, ng=0, error=0
    logical :: active=.false.
    integer :: ids(capacity)=0, lat(capacity)=0, lon(capacity)=0, counts(capacity)=0
  end type
contains
  subroutine contract_begin(ctx,epoch,mode,levels,ids,lat,lon)
    type(contract_context), intent(inout) :: ctx
    integer, intent(in) :: epoch,mode,levels,ids(:),lat(:),lon(:)
    integer :: i,j
    if(ctx%active .or. epoch<=ctx%consumed) then
      ctx%error=10
      return
    end if
    ctx%epoch=epoch;ctx%mode=mode;ctx%levels=levels;ctx%ng=size(ids)
    ctx%counts=0;ctx%ids=0;ctx%lat=0;ctx%lon=0;ctx%error=0;ctx%active=.true.
    if(epoch<0 .or. mode<0 .or. mode>1 .or. levels<1 .or. levels>3) ctx%error=12
    if(ctx%ng<1 .or. ctx%ng>capacity .or. size(lat)/=ctx%ng .or. size(lon)/=ctx%ng) then
      ctx%error=12
      return
    end if
    ctx%ids(:ctx%ng)=ids;ctx%lat(:ctx%ng)=lat;ctx%lon(:ctx%ng)=lon
    if(any(ids<=0) .or. any(lat<=0) .or. any(lon<=0)) ctx%error=12
    do i=1,ctx%ng
      do j=1,i-1
        if(ids(i)==ids(j)) ctx%error=12
        if(lat(i)==lat(j) .and. lon(i)==lon(j)) ctx%error=12
      end do
    end do
  end subroutine

  subroutine contract_observe(ctx,ncol,gcol,lat,lon)
    type(contract_context), intent(inout) :: ctx
    integer, intent(in) :: ncol,gcol(:),lat(:),lon(:)
    integer :: i,j,slot
    if(.not.ctx%active) ctx%error=10
    if(ctx%error/=0) return
    if(ncol<0 .or. ncol>size(gcol) .or. ncol>size(lat) .or. ncol>size(lon)) then
      ctx%error=12
      return
    end if
    do i=1,ncol
      slot=0
      do j=1,ctx%ng
        if(gcol(i)==ctx%ids(j)) slot=j
      end do
      if(slot==0) then
        ctx%error=12
        return
      end if
      if(lat(i)/=ctx%lat(slot) .or. lon(i)/=ctx%lon(slot)) then
        ctx%error=14
        return
      end if
      if(ctx%counts(slot)>=capacity) then
        ctx%error=12
        return
      end if
      ctx%counts(slot)=ctx%counts(slot)+1
    end do
  end subroutine

  subroutine contract_collect(ctx,comm,status,global_counts)
    type(contract_context), intent(inout) :: ctx
    integer, intent(in) :: comm
    integer, intent(out) :: status,global_counts(capacity)
    integer :: header(4),lo(4),hi(4),table(3*capacity),tlo(3*capacity),thi(3*capacity),ierr
    global_counts=-1
    if(.not.ctx%active) ctx%error=10
    header=[ctx%epoch,ctx%mode,ctx%ng,ctx%levels]
    ! Every caller executes the same fixed-size prefix, even on local errors.
    call MPI_Allreduce(ctx%error,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(header,lo,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(header,hi,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    ctx%active=.false.;ctx%consumed=max(ctx%consumed,ctx%epoch)
    if(status/=0) return
    if(any(lo/=hi)) then
      status=11
      return
    end if
    ! Compare the independently supplied expectation on all ranks as well.
    table=[ctx%ids,ctx%lat,ctx%lon]
    call MPI_Allreduce(table,tlo,3*capacity,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(table,thi,3*capacity,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(any(tlo/=thi)) then
      status=14
      return
    end if
    global_counts=0
    call MPI_Allreduce(ctx%counts,global_counts,ctx%ng,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(global_counts(:ctx%ng)/=1)) status=13
  end subroutine
end module qbo_contract
