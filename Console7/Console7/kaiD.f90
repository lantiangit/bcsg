!c*****************************************************************************
!c kai * d
!c*****************************************************************************
	   module kaiD_mod
        contains
    
       subroutine kaiD(xn)
	   use model_mod
       implicit none
       
       complex, intent(out) :: xn(mx,mz)
       complex :: xks(mx,mz)
       
       integer :: i,j,k
!c
!c the x,y,z component
!c
       
       do i=1,mx
		   do k=1,mz
			  xn(i,k)=xn(i,k)*xks(i,k)
		   end do
       end do

       
       end subroutine kaiD
	   
	   end module kaiD_mod
