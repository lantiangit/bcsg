    module FFTW3
        use, intrinsic :: iso_c_binding
        include 'fftw3.f03'
    end module
        
    
!************************************************************************
! This code computes the 2D fft of the Green's function in spatial domain.
! Input parameters:
!     The 8 components of Green's function in the common block define in 
!     layer.h
! Output parameters:
!     After finish the 2-D fft, write them back to the common block defined
!     in layer.h.
!************************************************************************


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   2D INVERSE FOURIER TRANSFORM      CC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	module ffts_mod
	contains
      subroutine in2dfft(xopen)
	  !use fftc_mod
	  use model_mod
      use FFTW3
	  implicit none 

      complex, intent(inout) :: xopen(nx,nz)
     ! complex :: xx1(nx),xx3(nz)
      integer :: i, k 
	  real :: acz, acx
      type(C_PTR) :: plan
      complex(C_DOUBLE_COMPLEX),allocatable, dimension(:) :: xx3_in, xx3_out
      complex(C_DOUBLE_COMPLEX),allocatable, dimension(:) :: xx1_in, xx1_out
      allocate(xx3_in(1:nz))
      allocate(xx3_out(1:nz))
      allocate(xx1_in(1:nx))
      allocate(xx1_out(1:nx))
	  acz = 1.0/real(nz)
      acx = 1.0/real(nz)

         do i=1,nx
            do k=1,nz
				xx3_in(k)=xopen(i,k)
            enddo
            !call cfftb(nz,xx3,zwsave)
            plan = fftw_plan_dft_1d(nz, xx3_in, xx3_out, FFTW_BACKWARD ,FFTW_ESTIMATE) !FFTW_BACKWARD


            call fftw_execute_dft(plan, xx3_in, xx3_out)
            call fftw_destroy_plan(plan) 
            do k=1,nz
				xopen(i,k)=xx3_out(k)*acz
            enddo
         enddo


      do k=1,nz
   
         do i=1,nx
			xx1_in(i)=xopen(i,k)
          enddo
  
          !call cfftb(nx,xx1,xwsave)
           plan = fftw_plan_dft_1d(nz, xx1_in, xx1_out, FFTW_BACKWARD ,FFTW_ESTIMATE) !FFTW_BACKWARD


            call fftw_execute_dft(plan, xx1_in, xx1_out)
            call fftw_destroy_plan(plan)
           do i=1,nx
			xopen(i,k)=xx1_out(i)*acx
           enddo
       enddo
      end subroutine in2dfft

!....................
!  2D FORWARD FFT
!...................

      subroutine fw2dfft(xopen)
      use model_mod
	  !use fftc_mod
      use FFTW3
	  implicit none

      complex, intent(inout) :: xopen(nx,nz)
      !complex :: xx1(nx),xx3(nz)
	  !complex :: xx1(nx)
      integer :: i, k
      type(C_PTR) :: plan
      complex(C_DOUBLE_COMPLEX),allocatable, dimension(:) :: xx3_in, xx3_out
      complex(C_DOUBLE_COMPLEX),allocatable, dimension(:) :: xx1_in, xx1_out
      allocate(xx3_in(1:nz))
      allocate(xx3_out(1:nz))
      allocate(xx1_in(1:nx))
      allocate(xx1_out(1:nx))

		 do i=1,nx
			do k=1,nz
				xx3_in(k)=xopen(i,k)
			enddo
			 !call cfftf(nz,xx3,zwsave)
            
             plan = fftw_plan_dft_1d(nz, xx3_in, xx3_out, FFTW_FORWARD,FFTW_ESTIMATE) !FFTW_BACKWARD


            call fftw_execute_dft(plan, xx3_in, xx3_out)
            call fftw_destroy_plan(plan) 
             
             
			do k=1,nz
				xopen(i,k)=xx3_out(k)
			enddo
		 enddo


      do k=1,nz
		do i=1,nx
			xx1_in(i)=xopen(i,k)
        enddo
        
             plan = fftw_plan_dft_1d(nx, xx1_in,xx1_out, FFTW_FORWARD,FFTW_ESTIMATE) !FFTW_BACKWARD


            call fftw_execute_dft(plan, xx1_in, xx1_out)
            call fftw_destroy_plan(plan) 
		 !call cfftf(nx,xx1,xwsave)
		do i=1,nx
			xopen(i,k)=xx1_out(i)
		enddo
      enddo

   
      end subroutine fw2dfft


!c************************************************************************
!c This program finish the function of three-dimensional fftshift
!c It is used in computing the Green's function
!c Input parameters: 
!c      Nx, Ny, Nz--the dimension of the max
!c      cA--the matrix to be fftshifted
!c Output parameters:
!c      cA--the matrix shifted
!c************************************************************************c
      subroutine fftshift(ca)
      use model_mod
      implicit none
 

      complex, intent(inout) :: ca(nx, nz)

      integer :: i, j, k, inc
      complex :: xx3(nz)


	 do i=1,nx
		do k=1,nz
		   inc=k+nz/2
		   if(inc > nz) then
		   inc=inc-nz
		   endif
		   xx3(inc)=ca(i,k)
		enddo
		do k=1,nz
		   ca(i,k)=xx3(k)
		enddo
	 enddo


      call fftshiftxy(ca)
      
      end subroutine fftshift
!c************************************************************************
!c This program finish the function of two-dimensional fftshift in x and y
!c direction.
!c It is used in computing the Green's function
!c Input parameters: 
!c      Nx, Ny, Nz--the dimension of the max
!c      cA--the matrix to be fftshifted
!c Output parameters:
!c      cA--the matrix shifted
!c************************************************************************c

      subroutine fftshiftxy(ca)
      use model_mod
      implicit none
	  
      complex, intent(inout) :: ca(nx,nz)

      integer i, j, k, inc

      complex xx1(nx)

        do k=1,nz


         do i=1,nx
			 inc=i+nx/2
			 if(inc > nx) inc=inc-nx
				xx1(inc)=ca(i,k)
         enddo

          do i=1,nx
			ca(i,k)=xx1(i)
          enddo

        enddo
  
      end subroutine fftshiftxy


!C*******************************************************************
!C     2D FORWARD FFT for the Green's function of the correlation term
!c     The difference from the fw2dfft1 is it do adjustment along 
!c     z direction.
!C*******************************************************************

      subroutine fw2dfft1(xopen)
      !use fftc_mod
	  use model_mod
      use FFTW3
	  implicit none 


      complex, intent(inout) :: xopen(nx,nz)
      !complex :: xx1(nx), xx3(nz)
      integer :: i, k
      type(C_PTR) :: plan
      complex(C_DOUBLE_COMPLEX),allocatable, dimension(:) :: xx3_in, xx3_out
      complex(C_DOUBLE_COMPLEX),allocatable, dimension(:) :: xx1_in, xx1_out
      allocate(xx3_in(1:nz))
      allocate(xx3_out(1:nz))
      allocate(xx1_in(1:nx))
      allocate(xx1_out(1:nx))

         do i=1,nx
            do k=1,nz
				xx3_in(k)=xopen(i,k)
            enddo
            !call cfftf(nz,xx3,zwsave)
             plan = fftw_plan_dft_1d(nz, xx3_in,xx3_out, FFTW_FORWARD,FFTW_ESTIMATE) !FFTW_BACKWARD


            call fftw_execute_dft(plan, xx3_in, xx3_out)
            call fftw_destroy_plan(plan) 
            xopen(i,1)=xx3_out(1)
            do k=2,nz
               xopen(i,k)=xx3_out(nz+2-k)  
            enddo						  
         enddo


      do k=1,nz
		do i=1,nx
			xx1_in(i)=xopen(i,k)
		enddo
		!call cfftf(nx,xx1,xwsave)
        plan = fftw_plan_dft_1d(nz, xx1_in,xx1_out, FFTW_FORWARD,FFTW_ESTIMATE) !FFTW_BACKWARD


            call fftw_execute_dft(plan, xx1_in, xx1_out)
            call fftw_destroy_plan(plan) 
		do i=1,nx
			xopen(i,k)=xx1_out(i)
		enddo
      enddo      
      end subroutine fw2dfft1
	  
	  
	  
	  
	  
!	  c***********************************************************************
!c  Operation: L[aE]---->bE
!c  
!c Input parameter: aE(mx,mz)
!c     
!c Ouput parameter: bE(mx,mz)
!c**********************************************************************

      subroutine Loperator(dxz,aE,bE)
		use para_mod
        use model_mod
        use green_global_parameters
		use geo_background

		implicit none


		integer :: i,k

		!c---the Variables for BCGSTAB

		real, intent(in) :: dxz
		complex, intent(in) :: aE(mx,mz)
		complex, intent(out) :: bE(mx,mz)
		
		real :: sc
		complex :: cJ(nx,nz),temp
		complex :: fopenm(nx,nz), fopenp(nx,nz)
		!type(C_PTR) :: plan
		!complex(C_DOUBLE_COMPLEX), dimension(nx,nz) :: fopenm, fopenp

        sc=real(nx*nz)
!---------------------------------------------
		do k=1,nz
			do i=1,nx
				fopenm(i,k)=cmplx(0.0,0.0)
				fopenp(i,k)=cmplx(0.0,0.0)
				cJ(i,k)=cmplx(0.0,0.0)
			enddo
		enddo

		do k=1,mz
			do i=1,mx
				cJ(i,k)=aE(i,k)*ckai(i,k)
				fopenm(i,k)=cJ(i,k)
			enddo
		enddo
         
  

!---------convolution----------minus part------  
		!plan = fftw_plan_dft_2d(nz,nx, fopenm, fopenm, FFTW_FORWARD,FFTW_ESTIMATE)
		!call fftw_execute_dft(plan, fopenm, fopenm)
		!call fftw_destroy_plan(plan)

		call fw2dfft(fopenm)
    
     
		do k=1,nz
			do i=1,nx
				fopenm(i,k)=fopenm(i,k)*fgreenkm(i,k)
			enddo
		enddo
         
         
		call in2dfft(fopenm)
		!plan = fftw_plan_dft_2d(nz,nx, fopenm, fopenm, FFTW_BACKWARD,FFTW_ESTIMATE)
		!call fftw_execute_dft(plan, fopenm, fopenm)
		!call fftw_destroy_plan(plan)
		


!---------correlation----plus part---------------  


		do k=1,mz+2
			do i=1,mx+2
				fopenp(i,k)=cJ(i,k)
			enddo
		enddo



         
		call fw2dfft1(fopenp)


		do k=1,nz
			do i=1,nx
				fopenp(i,k)=fopenp(i,k)*fgreenkp(i,k)
			enddo
		enddo

		call in2dfft(fopenp)

!-----operator

		do k=1,mz
			do i=1,mx
				temp=(fopenm(i,k)+fopenp(i,k))*dxz  
				bE(i,k)=aE(i,k)-xk02*temp
			enddo
		enddo
      
      end subroutine Loperator  
	  
	end module ffts_mod
