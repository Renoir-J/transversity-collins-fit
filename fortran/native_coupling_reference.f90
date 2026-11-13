! Standalone comparator using unchanged pinned HOPPET running-coupling code.
! A fresh process owns every batch. No grid or splitting function is initialized.
program native_coupling_reference
  use iso_fortran_env, only: input_unit, output_unit, error_unit
  use, intrinsic :: ieee_arithmetic
  use types, only: dp
  use qcd_coupling, only: running_coupling, InitRunningCoupling, Value, Delete, &
       SetDefaultCouplingDt, DefaultCouplingDt
  implicit none
  character(len=256) :: magic, extra
  integer :: n,nf,nloop,i,ios
  real(dp) :: alpha0,q02,dt,saved_dt
  real(dp), allocatable :: q2(:),alpha(:)
  type(running_coupling) :: coupling

  read(input_unit,'(A)',iostat=ios) magic
  call require(ios==0 .and. trim(magic)=='CF-COUPLING-1','protocol header')
  read(input_unit,*,iostat=ios) n,nf,nloop
  call require(ios==0,'integer request header')
  call require(n>=1 .and. n<=20000,'query count')
  call require(nf>=3 .and. nf<=6,'fixed nf')
  call require(nloop==1 .or. nloop==2,'loop order')
  read(input_unit,*,iostat=ios) alpha0,q02,dt
  call require(ios==0,'coupling request header')
  call require(all(ieee_is_finite([alpha0,q02,dt])),'finite settings')
  call require(alpha0>0.0_dp .and. alpha0<=0.4_dp,'alpha0 domain')
  call require(q02>=2.4_dp .and. q02<=1.0e8_dp,'Q02 domain')
  call require(dt>=0.0001_dp .and. dt<=0.2_dp,'dt domain')
  allocate(q2(n),alpha(n))
  do i=1,n
    read(input_unit,*,iostat=ios) q2(i)
    call require(ios==0,'query row')
  end do
  call require(all(ieee_is_finite(q2)),'finite query scales')
  call require(all(q2>=q02) .and. all(q2<=1.0e8_dp),'forward query domain')
  do
    read(input_unit,'(A)',iostat=ios) extra
    if (ios<0) exit
    call require(ios==0,'request trailer')
    call require(len_trim(extra)==0,'unexpected request trailer')
  end do

  saved_dt=DefaultCouplingDt()
  call SetDefaultCouplingDt(dt)
  call InitRunningCoupling(coupling,alfas=alpha0,Q=sqrt(q02),nloop=nloop,fixnf=nf, &
       Qmax=max(100.0_dp,sqrt(maxval(q2)))*1.01_dp)
  do i=1,n
    alpha(i)=Value(coupling,sqrt(q2(i)))
  end do
  call Delete(coupling)
  call SetDefaultCouplingDt(saved_dt)
  call require(all(ieee_is_finite(alpha)) .and. all(alpha>0.0_dp),'native finite positive values')
  ! Emit only after the entire batch has completed and its object was deleted.
  write(output_unit,'(A)') 'CF-COUPLING-1'
  write(output_unit,'(I0,1X,I0,1X,I0)') n,nf,nloop
  do i=1,n
    write(output_unit,'(ES26.17E3)') alpha(i)
  end do
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(error_unit,'(A)') 'REFERENCE_INPUT_OR_RESULT_REJECTED: '//label
      stop 2
    end if
  end subroutine require
end program native_coupling_reference
