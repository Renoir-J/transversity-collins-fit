! Fixed Float64 ABI for already validated operators. This does not build a
! transversity HOPPET kernel and is never claimed to be KPSY's missing code.
module collins_operator_bridge
  use iso_c_binding
  implicit none
contains
  subroutine cf_apply_operator(nout, nin, nbatch, op, initial, result, status) bind(C)
    integer(c_int), value :: nout, nin, nbatch
    real(c_double), intent(in) :: op(nout,nin), initial(nin,nbatch)
    real(c_double), intent(out) :: result(nout,nbatch)
    integer(c_int), intent(out) :: status
    integer :: i, j, k
    status = 1_c_int
    if (nout <= 0 .or. nin <= 0 .or. nbatch <= 0) return
    result = 0.0_c_double
    do j = 1, nbatch
      do k = 1, nin
        do i = 1, nout
          result(i,j) = result(i,j) + op(i,k) * initial(k,j)
        end do
      end do
    end do
    status = 0_c_int
  end subroutine
end module

