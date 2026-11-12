! Process-isolated driver for the unchanged, locked legacy source.
program check_dssv08
  implicit none
  integer :: count, i
  real(8) :: x, q2, uv, dv, ub, db, str, glu
  call DSSVINI(0)
  read(*,*) count
  do i=1,count
    read(*,*) x,q2
    call DSSVFIT(x,q2,uv,dv,ub,db,str,glu)
    write(*,'(a,1p,6(es25.16,1x))') 'RESULT ',uv,dv,ub,db,str,glu
  end do
end program
