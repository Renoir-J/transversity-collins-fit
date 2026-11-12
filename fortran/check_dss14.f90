! Process-isolated driver; every run starts a fresh COMMON/SAVE state.
program check_dss14
  implicit none
  integer :: count, i, charge, fini
  real(8) :: z,q2,u,ub,d,db,s,sb,c,b,glu
  common /fragini/ fini
  fini=0
  read(*,*) count
  do i=1,count
    read(*,*) charge,z,q2
    call fDSSH(0,0,charge,1,z,q2,u,ub,d,db,s,sb,c,b,glu)
    write(*,'(a,1p,9(es25.16,1x))') 'RESULT ',u,ub,d,db,s,sb,c,b,glu
  end do
end program
