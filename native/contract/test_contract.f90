program test_contract
  use mpi
  use qbo_contract
  implicit none
  type(contract_context) :: ctx
  integer, parameter :: ids(8)=[41,7,103,2,88,19,54,66]
  integer, parameter :: lat(8)=[2,1,2,1,2,1,2,1]
  integer, parameter :: lon(8)=[3,4,1,2,4,1,2,3]
  integer, parameter :: expected(12)=[0,0,0,13,13,12,14,11,11,11,10,10]
  integer :: rank,np,ierr,k,i,j,n,ep,mode,levels,code,lo,hi,fail,allfail
  integer :: g(10),la(10),ln(10),order(8),counts(8),padding_ok,allpadding
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  fail=0
  if (np/=2) then
    call MPI_Finalize(ierr)
    stop 2
  end if
  do k=1,12
    ep=100*k; mode=1; levels=3
    if(k==8 .and. rank==1) ep=ep+1
    if(k==9 .and. rank==1) mode=0
    if(k==10 .and. rank==1) levels=2
    call contract_begin(ctx,ep,mode,levels,ids,lat,lon)
    n=4
    do i=1,4
      order(i)=2*i-1+rank
    end do
    if(k==2) order(1:4)=order(4:1:-1)
    if(k==3) then
      n=0
      if(rank==0) then
        n=8
        order=[1,2,3,4,5,6,7,8]
      end if
    end if
    if(k==4 .and. rank==1) order(4)=1
    if(k==5 .and. rank==1) n=3
    g=-999999; la=-999999; ln=-999999
    do i=1,n
      j=order(i);g(i)=ids(j);la(i)=lat(j);ln(i)=lon(j)
    end do
    if(k==6 .and. rank==1) g(1)=999999
    if(k==7 .and. rank==1) la(1)=3
    ! Two observations exercise chunk accumulation; padding must never be read.
    j=n/2
    call contract_observe(ctx,j,g(:j),la(:j),ln(:j))
    call contract_observe(ctx,n-j,g(j+1:),la(j+1:),ln(j+1:))
    call contract_collect(ctx,MPI_COMM_WORLD,code,counts)
    if(k<=3 .or. k>=11) then
      if(code/=0 .or. any(counts/=1)) fail=fail+1
    end if
    if(k==11) call contract_collect(ctx,MPI_COMM_WORLD,code,counts)
    if(k==12) then
      call contract_begin(ctx,ep,mode,levels,ids,lat,lon)
      call contract_collect(ctx,MPI_COMM_WORLD,code,counts)
    end if
    padding_ok=0
    if(all(g(n+1:)==-999999) .and. all(la(n+1:)==-999999) .and. all(ln(n+1:)==-999999)) padding_ok=1
    call MPI_Allreduce(code,lo,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(code,hi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(padding_ok,allpadding,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    if(lo/=expected(k) .or. hi/=expected(k) .or. allpadding/=1) fail=fail+1
    if(rank==0) write(*,'(a,6(1x,i0))') 'CASE',k,lo,hi,expected(k),allpadding,sum(counts)
  end do
  call MPI_Allreduce(fail,allfail,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(rank==0) then
    write(*,'(a,1x,i0)') 'CONTRACT_FAILURES',allfail
    if(allfail==0) write(*,'(a)') 'CONTRACT_SYNTHETIC_PASS'
  end if
  call MPI_Finalize(ierr)
  if(allfail/=0) stop 1
end program test_contract
