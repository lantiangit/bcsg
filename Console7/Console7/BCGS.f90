!c***********************************************************************
!c  Solve the operator equation by the CG, BCG, BCGS method
!c  and obtain the induced  electric current inside the scatterer
!c
!c Input parameter:
!c     cEi(3, mx+1, my+1, mz+1): the incident electric field
!c      mode: TE or TM mode, to generate different output file names
!c     ceps(mx+1, my+1, mz+1): the dielectric constant matrix for the
!c     computational domain
!c     ckai(mx, mz): the contrast matrix of the scatterer
!c     to the backgroud
!c     of layer nobj
!c     dx, dz: the increment in the x, y, z direction
!c
!c Ouput parameter:
!c     xn: the electric flux density inside the computational domain
!c**********************************************************************
	module bcgs_mod
	contains
      subroutine bcgs(dx,dz,xn)
      ! 
	  use model_mod
	  !use tds1_mod
	  !use contrast_mod
	  use para_mod
	  use ffts_mod
      use general_constants
      use geo_background
      use iter_control_mod
	  use field_mod
      implicit none


      real, intent(in) :: dx, dz	  
      complex, intent(out) :: xn(mx,mz)	  
      character(len=2) :: mode
      integer :: iexx,i,j,k
      integer :: istep, memm, itabnormal



	  real :: dxz
      real :: rerror0,rnormb,rerror,rnorms



      complex :: pn(mx,mz)
      complex :: rn(mx,mz)
      complex :: vn(mx,mz)

      complex :: tn(mx,mz)

      complex :: rn_hat(mx,mz)
      complex :: dabu,c1,c2,beta,rho1,rho

      character(len=20) :: filename


!c initialization for iteration
!c some constants for the operator equation.
!c
       xk02=xomega**2*xmu0*cer(nobj)*xepsilon0
  
       dxz=dx*dz

       do i=1,mx
		   do k=1,mz
			  rn(i,k)=cei(i,k)
			  pn(i,k)=0.0
			  rn_hat(i,k)=conjg(rn(i,k))
			  vn(i,k)=0.0
			  xn(i,k)=0.0
		   enddo
       enddo


       rho=(1.,0)
       beta=(1.,0)
       dabu=(1.,0)



!cc  Coefs used to control Iterations!

       itabnormal=1000
       memm=0
       rError0=1.0e9
       rNormB=0
       do k=1,mz
		   do i=1,mx
			rNormB=rNormB+abs(rn(i,k))**2
		   enddo
       enddo


!c
!c     starting iteration
!c

        do istep=1,istepmax

			rho1=0.0
			do k=1,mz
				do i=1,mx
					rho1=rho1+rn_hat(i,k)*rn(i,k)
				enddo
			enddo


			beta=rho1/rho*beta/dabu

			do k=1,mz
				do i=1,mx
					pn(i,k)=rn(i,k)+beta*(pn(i,k)-dabu*vn(i,k))
				enddo
			enddo


			call Loperator(dxz,pn,vn)


			beta=0.0
			do k=1,mz
				do i=1,mx
					beta=beta+rn_hat(i,k)*vn(i,k)
				enddo
			enddo


			beta=rho1/beta
			rNormS=0.00


			do k=1,mz
				do i=1,mx
					rn(i,k)=rn(i,k)-beta*vn(i,k)
					rNormS=rNormS+abs(rn(i,k))**2
				enddo
			enddo


			call Loperator(dxz,rn,tn)

			c1=0.0
			c2=0.0
			do k=1,mz
				do i=1,mx
					c1=c1+conjg(tn(i,k))*rn(i,k)
					c2=c2+conjg(tn(i,k))*tn(i,k)
				enddo
			enddo


			dabu=c1/c2
			do k=1,mz
				do i=1,mx
					xn(i,k)=xn(i,k)+beta*pn(i,k)+dabu*rn(i,k)
					rn(i,k)=rn(i,k)-dabu*tn(i,k)
				enddo
			enddo


			rError=sqrt(rNormS/rNormB)

			if(rError < rError0) then
				memm=memm+1
				rError0=rError
			endif



			if (rError <= ctrmax .or. memm >= itabnormal) exit

			rho=rho1

		enddo



      end subroutine bcgs
	  
	  end module bcgs_mod
