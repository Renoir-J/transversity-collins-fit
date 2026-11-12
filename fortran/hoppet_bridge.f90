! Project-owned adapter to unchanged, pinned HOPPET. Reconstructed LO only.
! The numerical job is a single isolated process. No callbacks into Python/Julia.
module collins_hoppet_bridge
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use hoppet
  use sort, only: indexx
  use convolution_communicator, only: cc_piece, cc_REAL, cc_VIRT, cc_REALVIRT, cc_DELTA
  implicit none
  private
  type(grid_def), save :: grid
  type(grid_conv), allocatable, save :: kernels(:)
  real(dp), allocatable, save :: xs(:), unique_xs(:)
  integer, allocatable, save :: node_map(:), representative(:)
  real(dp), allocatable, save :: times(:)
  integer(c_int64_t), save :: generation=0_c_int64_t
  logical, save :: active=.false.
  public :: cf_hoppet_abi, cf_hoppet_open, cf_hoppet_open_batched, cf_hoppet_close, cf_hoppet_grid
  public :: cf_hoppet_forward, cf_hoppet_adjoint, cf_hoppet_generator, cf_hoppet_evaluate
  public :: cf_hoppet_info, cf_hoppet_rows, cf_hoppet_collins, cf_hoppet_collins_rows
contains
  integer(c_int) function cf_hoppet_abi() bind(C)
    cf_hoppet_abi=1_c_int
  end function

  ! Regular matching R(v)=2*CF*v*ln(v)/(1-v), with the native extra v.
  ! ln(1)=0 removes the plus subtraction and its delta contribution.
  function collins_regular_kernel(y) result(value)
    real(dp), intent(in) :: y
    real(dp) :: value,z,ratio
    value=0.0_dp
    if (cc_piece/=cc_REAL .and. cc_piece/=cc_REALVIRT) return
    z=exp(-y)
    if (y==0.0_dp) then
      ratio=1.0_dp
    else
      ratio=y/(2.0_dp*exp(-0.5_dp*y)*sinh(0.5_dp*y))
    end if
    value=-(8.0_dp/3.0_dp)*z*z*ratio
  end function

  ! Plus acts on 1/(1-z). HOPPET's non-delta pieces carry one extra z.
  function transversity_kernel(y) result(value)
    real(dp), intent(in) :: y
    real(dp) :: value, z, omz
    real(dp), parameter :: cft=4.0_dp/3.0_dp
    z=exp(-y)
    select case(cc_piece)
    case(cc_REAL)
      omz=2.0_dp*exp(-0.5_dp*y)*sinh(0.5_dp*y)
      value=2.0_dp*cft*z*z/omz
    case(cc_VIRT)
      omz=2.0_dp*exp(-0.5_dp*y)*sinh(0.5_dp*y)
      value=-2.0_dp*cft*z/omz
    case(cc_REALVIRT)
      value=-2.0_dp*cft*z
    case(cc_DELTA)
      value=1.5_dp*cft
    case default
      value=0.0_dp
    end select
  end function

  subroutine cf_hoppet_open(ymax,lower_x,dy,order,layout,eps,du,nf,ntau,tau,n,handle,status) bind(C)
    real(c_double), value :: ymax,lower_x,dy,eps,du
    integer(c_int), value :: order,layout,nf,ntau
    real(c_double), intent(in) :: tau(ntau)
    integer(c_int), intent(out) :: n,status
    integer(c_int64_t), intent(out) :: handle
    call open_context(ymax,lower_x,dy,order,layout,eps,du,nf,ntau,tau,n,handle,status,.false.)
  end subroutine

  subroutine cf_hoppet_open_batched(ymax,lower_x,dy,order,layout,eps,du,nf,ntau,tau,n,handle,status) bind(C)
    real(c_double), value :: ymax,lower_x,dy,eps,du
    integer(c_int), value :: order,layout,nf,ntau
    real(c_double), intent(in) :: tau(ntau)
    integer(c_int), intent(out) :: n,status
    integer(c_int64_t), intent(out) :: handle
    call open_context(ymax,lower_x,dy,order,layout,eps,du,nf,ntau,tau,n,handle,status,.true.)
  end subroutine

  subroutine open_context(ymax,lower_x,dy,order,layout,eps,du,nf,ntau,tau,n,handle,status,batched)
    real(c_double), intent(in) :: ymax,lower_x,dy,eps,du
    integer(c_int), intent(in) :: order,layout,nf,ntau
    real(c_double), intent(in) :: tau(ntau)
    integer(c_int), intent(out) :: n,status
    integer(c_int64_t), intent(out) :: handle
    logical, intent(in) :: batched
    type(dglap_holder), target :: holder
    type(running_coupling), target :: auxiliary
    type(evln_operator) :: evop
    type(grid_conv) :: generator
    type(grid_def) :: subgrids(7)
    real(dp), parameter :: limits(6)=[0.01_dp,0.03_dp,0.08_dp,0.2_dp,0.5_dp,2.0_dp]
    integer :: k,ny,i,u,ii
    integer, allocatable :: sorted(:), targets(:)
    real(dp), allocatable :: ordered(:)
    real(dp) :: saved_du
    status=1_c_int; n=0_c_int; handle=0_c_int64_t
    if (active) then
      status=2_c_int
      return
    end if
    if (.not. all(ieee_is_finite([ymax,lower_x,dy,eps,du]))) return
    if (lower_x<=0.0_dp .or. lower_x>=1.0_dp) return
    if (abs(log(lower_x)+ymax)>8*epsilon(ymax)*max(1.0_dp,abs(ymax))) return
    if (ymax<=0 .or. ymax>30 .or. dy<=0 .or. eps<=0 .or. eps>1e-4_dp .or. du<=0 .or. du>0.2_dp) return
    if (order > -2 .or. order < -10 .or. nf<3 .or. nf>6 .or. ntau<1 .or. ntau>4096) return
    if (layout/=0 .and. layout/=1) return
    if (layout==1 .and. (ymax<2.5_dp .or. dy>0.2_dp)) return
    if (.not. all(ieee_is_finite(tau)) .or. any(abs(tau)>1.0_dp)) return
    if (ymax/dy>20000.0_dp) return
    ny=nint(ymax/dy)
    if (ny<abs(order)+1) return
    if (layout==0) then
      call InitGridDef(grid,dy=dy,ymax=ymax,order=order,eps=eps)
    else
      do k=1,6
        call InitGridDef(subgrids(k),dy=dy/(3.0_dp**(7-k)),ymax=limits(k),order=order,eps=eps)
      end do
      call InitGridDef(subgrids(7),dy=dy,ymax=ymax,order=order,eps=eps)
      call InitGridDef(grid,subgrids,locked=.true.)
      call Delete(subgrids)
    end if
    allocate(xs(0:grid%ny),times(ntau),kernels(ntau))
    xs=xValues(grid); times=tau
    ! The declared fraction boundary is canonical. log/exp round trips need
    ! not reproduce e.g. 0.05 bit-for-bit; this is not a source extrapolation.
    xs(grid%ny)=lower_x
    ! Preserve exact x values, including distinct floating point knots. Shared
    ! nodes are sampled once; choose their finest native representation.
    allocate(sorted(size(xs)),ordered(size(xs)),targets(size(xs)),node_map(0:grid%ny))
    call indexx(xs,sorted)
    u=0
    do i=1,size(xs)
      ii=sorted(i)-1
      if (i==1) then
        u=1; ordered(u)=xs(ii); targets(u)=ii
      else if (xs(ii)/=ordered(u)) then
        u=u+1; ordered(u)=xs(ii); targets(u)=ii
      else
        targets(u)=min(targets(u),ii)
      end if
      node_map(ii)=u
    end do
    allocate(unique_xs(u),representative(u))
    unique_xs=ordered(:u); representative=targets(:u)
    call qcd_SetNf(nf)
    call InitDglapHolder(grid,holder,nloop=1,nflo=nf,nfhi=nf)
    call Delete(holder%P_LO)
    call InitSplitMat(grid,holder%P_LO,nf)
    call InitGridConv(grid,generator,transversity_kernel)
    ! Every quark/antiquark evolves independently, with no gluon component.
    call InitGridConv(holder%P_LO%qq,generator)
    call InitGridConv(holder%P_LO%NS_plus,generator)
    call InitGridConv(holder%P_LO%NS_minus,generator)
    call InitGridConv(holder%P_LO%NS_V,generator)
    call InitRunningCoupling(auxiliary,alfas=1.0_dp,Q=100.0_dp,nloop=0,fixnf=nf)
    ! In this pinned revision ev_evolve accepts du but omits it when invoking
    ! ev_SetModuleConsts. Use the public default setter, then restore it.
    saved_du=DefaultEvolutionDu()
    call SetDefaultEvolutionDu(du)
    if (batched) then
      call prepare_batched_kernels(holder,auxiliary,nf,du)
    else
      do k=1,ntau
        if (tau(k)==0.0_dp) cycle
        call InitEvlnOperator(holder,evop,auxiliary,100.0_dp,100.0_dp*exp(pi*tau(k)),nloop=1,du=du)
        call InitGridConv(kernels(k),evop%P%NS_plus)
        call Delete(evop)
      end do
    end if
    call SetDefaultEvolutionDu(saved_du)
    call Delete(generator)
    call Delete(holder)
    call Delete(auxiliary)
    generation=generation+1_c_int64_t
    active=.true.; n=u; handle=generation; status=0_c_int
  end subroutine

  ! Advance retained derived probes, not an evop that EvolveGeneric would reset.
  subroutine prepare_batched_kernels(holder,auxiliary,nf,du)
    type(dglap_holder), intent(in), target :: holder
    type(running_coupling), intent(in), target :: auxiliary
    integer, intent(in) :: nf
    real(dp), intent(in) :: du
    real(dp), pointer :: working(:,:,:)=>null(), exported(:,:,:)=>null()
    type(split_mat) :: snapshot
    integer, allocatable :: permutation(:)
    integer :: direction,it,k,previous,probe,position
    real(dp) :: Q_previous,Q_current
    logical :: duplicate
    allocate(permutation(size(times)))
    call indexx(times,permutation)
    do direction=1,-1,-2
      previous=0; Q_previous=100.0_dp
      do it=1,size(times)
        position=it
        if (direction==-1) position=size(times)+1-it
        k=permutation(position)
        if (direction==1 .and. times(k)<=0.0_dp) cycle
        if (direction==-1 .and. times(k)>=0.0_dp) cycle
        if (.not. associated(working)) call GetDerivedSplitMatProbes(grid,nf,working)
        duplicate=.false.
        if (previous>0) duplicate=times(k)==times(previous)
        if (duplicate) then
          call InitGridConv(kernels(k),kernels(previous))
        else
          Q_current=100.0_dp*exp(pi*times(k))
          call DisableGridLocking()
          do probe=1,size(working,3)
            call EvolveGeneric(holder,auxiliary,Q_previous,Q_current,pdf=working(:,:,probe),nloop=1,du=du)
          end do
          allocate(exported(lbound(working,1):ubound(working,1),lbound(working,2):ubound(working,2), &
                            lbound(working,3):ubound(working,3)))
          exported=working
          call AllocSplitMat(grid,snapshot,nf)
          ! This call consumes exported and restores the grid-lock override.
          call SetDerivedSplitMat(snapshot,exported)
          nullify(exported)
          call InitGridConv(kernels(k),snapshot%NS_plus)
          call Delete(snapshot)
          Q_previous=Q_current
        end if
        previous=k
      end do
      if (associated(working)) then
        deallocate(working)
        nullify(working)
        call RestoreGridLocking()
      end if
    end do
  end subroutine

  logical function valid(handle,n)
    integer(c_int64_t), intent(in) :: handle
    integer(c_int), intent(in) :: n
    valid=.false.
    if (.not. active) return
    valid=handle==generation .and. n==size(unique_xs)
  end function

  subroutine cf_hoppet_close(handle,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), intent(out) :: status
    integer :: k
    status=3_c_int
    if (.not. active .or. handle/=generation) return
    do k=1,size(kernels)
      if (times(k)/=0.0_dp) call Delete(kernels(k))
    end do
    deallocate(kernels,xs,unique_xs,times,node_map,representative)
    call Delete(grid)
    active=.false.; status=0_c_int
  end subroutine

  subroutine cf_hoppet_grid(handle,n,x,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n
    real(c_double), intent(out) :: x(n)
    integer(c_int), intent(out) :: status
    status=3_c_int
    if (.not. valid(handle,n)) return
    x=unique_xs; status=0_c_int
  end subroutine

  subroutine cf_hoppet_forward(handle,n,nb,ntau,h,result,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,nb,ntau
    real(c_double), intent(in) :: h(n,nb)
    real(c_double), intent(out) :: result(n,nb,ntau)
    integer(c_int), intent(out) :: status
    real(dp), allocatable :: weighted(:), evolved(:)
    integer :: j,k
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (nb<1 .or. ntau/=size(times)) return
    if (.not. all(ieee_is_finite(h)) .or. any(h(n,:)/=0.0_dp)) return
    allocate(weighted(0:grid%ny),evolved(0:grid%ny))
    do j=1,nb
      weighted=xs*h(node_map,j)
      do k=1,ntau
        if (times(k)==0.0_dp) then
          result(:,j,k)=h(:,j)
        else
          evolved=kernels(k)*weighted
          result(:,j,k)=evolved(representative)/unique_xs
          result(n,j,k)=0.0_dp
        end if
      end do
    end do
    if (.not. all(ieee_is_finite(result))) return
    status=0_c_int
  end subroutine

  ! Transpose locked-grid decanting in reverse order, then reuse upstream
  ! Toeplitz convolutions on each native subgrid. No forward solver rewrite.
  recursive function convolution_transpose(kernel,g) result(r)
    type(grid_conv), intent(in) :: kernel
    real(dp), intent(in) :: g(0:)
    real(dp) :: r(0:size(g)-1), work(0:size(g)-1)
    integer :: sub,step,iy,low,high,coarse,fine
    if (kernel%grid%nsub==0) then
      work=kernel*g(size(g)-1:0:-1)
      r=work(size(g)-1:0:-1)
      return
    end if
    work=g
    do sub=kernel%grid%nsub,2,-1
      step=nint(kernel%grid%subgd(sub)%dy/kernel%grid%subgd(sub-1)%dy)
      do iy=0,kernel%grid%subgd(sub-1)%ny/step
        coarse=kernel%grid%subiy(sub)+iy
        fine=kernel%grid%subiy(sub-1)+iy*step
        work(fine)=work(fine)+work(coarse)
        work(coarse)=0.0_dp
      end do
    end do
    do sub=1,kernel%grid%nsub
      low=kernel%grid%subiy(sub); high=kernel%grid%subiy(sub+1)-1
      r(low:high)=convolution_transpose(kernel%subgc(sub),work(low:high))
    end do
  end function

  subroutine cf_hoppet_adjoint(handle,n,nb,ntau,g,result,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,nb,ntau
    real(c_double), intent(in) :: g(n,nb,ntau)
    real(c_double), intent(out) :: result(n,nb)
    integer(c_int), intent(out) :: status
    real(dp), allocatable :: native(:), evolved(:)
    integer :: j,k,i
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (nb<1 .or. ntau/=size(times)) return
    if (.not. all(ieee_is_finite(g))) return
    allocate(native(0:grid%ny),evolved(0:grid%ny))
    result=0.0_dp
    do j=1,nb
      do k=1,ntau
        if (times(k)==0.0_dp) then
          result(:,j)=result(:,j)+g(:,j,k)
        else
          native=0.0_dp
          native(representative)=g(:,j,k)/unique_xs
          evolved=convolution_transpose(kernels(k),native)
          do i=0,grid%ny
            result(node_map(i),j)=result(node_map(i),j)+xs(i)*evolved(i)
          end do
        end if
      end do
      result(n,j)=0.0_dp
    end do
    if (.not. all(ieee_is_finite(result))) return
    status=0_c_int
  end subroutine

  subroutine cf_hoppet_generator(handle,n,nb,h,result,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,nb
    real(c_double), intent(in) :: h(n,nb)
    real(c_double), intent(out) :: result(n,nb)
    integer(c_int), intent(out) :: status
    type(grid_conv) :: generator
    real(dp), allocatable :: weighted(:), evolved(:)
    integer :: j
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (nb<1) return
    if (.not. all(ieee_is_finite(h)) .or. any(h(n,:)/=0.0_dp)) return
    allocate(weighted(0:grid%ny),evolved(0:grid%ny))
    call InitGridConv(grid,generator,transversity_kernel)
    do j=1,nb
      weighted=xs*h(node_map,j)
      evolved=generator*weighted
      result(:,j)=evolved(representative)/unique_xs
      result(n,j)=0.0_dp
    end do
    call Delete(generator)
    if (.not. all(ieee_is_finite(result))) return
    status=0_c_int
  end subroutine

  subroutine cf_hoppet_evaluate(handle,n,nb,nq,h,xq,result,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,nb,nq
    real(c_double), intent(in) :: h(n,nb),xq(nq)
    real(c_double), intent(out) :: result(nq,nb)
    integer(c_int), intent(out) :: status
    real(dp), allocatable :: weighted(:)
    integer :: j,k
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (nb<1 .or. nq<1) return
    if (.not. all(ieee_is_finite(h)) .or. .not. all(ieee_is_finite(xq))) return
    if (any(h(n,:)/=0.0_dp) .or. any(xq<unique_xs(1)) .or. any(xq>1.0_dp)) return
    allocate(weighted(0:grid%ny))
    do j=1,nb
      weighted=xs*h(node_map,j)
      do k=1,nq
        if (xq(k)==1.0_dp) then
          result(k,j)=0.0_dp
        else
          result(k,j)=EvalGridQuant(grid,weighted,-log(xq(k)))/xq(k)
        end if
      end do
    end do
    if (.not. all(ieee_is_finite(result))) return
    status=0_c_int
  end subroutine

  recursive integer(c_int64_t) function coefficient_bytes(kernel) result(total)
    type(grid_conv), intent(in) :: kernel
    integer :: sub
    total=0_c_int64_t
    if (associated(kernel%conv)) total=int(size(kernel%conv),c_int64_t)*8_c_int64_t
    if (associated(kernel%subgc)) then
      do sub=1,size(kernel%subgc)
        total=total+coefficient_bytes(kernel%subgc(sub))
      end do
    end if
  end function

  subroutine cf_hoppet_info(handle,n,native_nodes,payload,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n
    integer(c_int), intent(out) :: native_nodes,status
    integer(c_int64_t), intent(out) :: payload
    integer :: k
    status=3_c_int; payload=0_c_int64_t; native_nodes=0_c_int
    if (.not. valid(handle,n)) return
    do k=1,size(kernels)
      if (times(k)/=0.0_dp) payload=payload+coefficient_bytes(kernels(k))
    end do
    native_nodes=grid%ny+1; status=0_c_int
  end subroutine
  ! Export one selected output row per scale with linear workspace. This is
  ! the same constrained adjoint, without a node-by-scale-by-scale cotangent.
  subroutine cf_hoppet_rows(handle,n,ntau,node,rows,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,ntau,node
    real(c_double), intent(out) :: rows(ntau,n)
    integer(c_int), intent(out) :: status
    real(dp), allocatable :: native(:), evolved(:)
    integer :: i,k
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (ntau/=size(times) .or. node<1 .or. node>n) return
    allocate(native(0:grid%ny),evolved(0:grid%ny))
    rows=0.0_dp
    do k=1,ntau
      if (times(k)==0.0_dp) then
        rows(k,node)=1.0_dp
      else
        native=0.0_dp
        native(representative(node))=1.0_dp/unique_xs(node)
        evolved=convolution_transpose(kernels(k),native)
        do i=0,grid%ny
          rows(k,node_map(i))=rows(k,node_map(i))+xs(i)*evolved(i)
        end do
      end if
      rows(k,n)=0.0_dp
    end do
    if (.not. all(ieee_is_finite(rows))) return
    status=0_c_int
  end subroutine

  ! Unweighted regular Collins matching on a single retained grid. This full
  ! native path is used to validate row preparation, not inside Julia fitting.
  subroutine cf_hoppet_collins(handle,n,nb,h,result,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,nb
    real(c_double), intent(in) :: h(n,nb)
    real(c_double), intent(out) :: result(n,nb)
    integer(c_int), intent(out) :: status
    type(grid_conv) :: regular
    real(dp), allocatable :: weighted(:),matched(:)
    integer :: j
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (nb<1 .or. grid%nsub/=0) return
    if (.not. all(ieee_is_finite(h)) .or. any(h(n,:)/=0.0_dp)) return
    allocate(weighted(0:grid%ny),matched(0:grid%ny))
    call InitGridConv(grid,regular,collins_regular_kernel)
    do j=1,nb
      weighted=xs*h(node_map,j)
      matched=regular*weighted
      result(:,j)=matched(representative)/unique_xs
      result(n,j)=0.0_dp
    end do
    call Delete(regular)
    if (.not. all(ieee_is_finite(result))) return
    status=0_c_int
  end subroutine

  ! Selected rows of R*U: apply R^T once, then the retained U(tau)^T.
  ! No composed dense matrix or node-by-scale-by-scale tensor is allocated.
  subroutine cf_hoppet_collins_rows(handle,n,ntau,node,rows,status) bind(C)
    integer(c_int64_t), value :: handle
    integer(c_int), value :: n,ntau,node
    real(c_double), intent(out) :: rows(ntau,n)
    integer(c_int), intent(out) :: status
    type(grid_conv) :: regular
    real(dp), allocatable :: native(:),cotangent(:),evolved(:)
    integer :: i,k
    status=3_c_int
    if (.not. valid(handle,n)) return
    status=1_c_int
    if (ntau/=size(times) .or. node<1 .or. node>n .or. grid%nsub/=0) return
    allocate(native(0:grid%ny),cotangent(0:grid%ny),evolved(0:grid%ny))
    native=0.0_dp
    native(representative(node))=1.0_dp/unique_xs(node)
    call InitGridConv(grid,regular,collins_regular_kernel)
    cotangent=convolution_transpose(regular,native)
    call Delete(regular)
    rows=0.0_dp
    do k=1,ntau
      if (times(k)==0.0_dp) then
        evolved=cotangent
      else
        evolved=convolution_transpose(kernels(k),cotangent)
      end if
      do i=0,grid%ny
        rows(k,node_map(i))=rows(k,node_map(i))+xs(i)*evolved(i)
      end do
      rows(k,n)=0.0_dp
    end do
    if (.not. all(ieee_is_finite(rows))) return
    status=0_c_int
  end subroutine
end module
