!c**************************************************************************
!c  This subroutine is to calculate the scattered electric field
!c  in Layer 1 from the inhomogeneous scattering object.
!c     Input parameters:
!c        dx,dz: cell size in x and z direction
!c        xx0,xz0: the location of the dipole source
!c        x_cell, z_cell: coordinates of the scattering object in x,y,z direction
!c        ckaiD: the electric flux density inside the target multiplied by 
!c               the contrast constant kai.
!c        file: Rloc.inp, the receiver locations in 3D rectangular coordinates
!c     Output parameters:
!c        file Escat.dat: the scattered electric field Ey
!c                 at the receiver locations 
!c  kinv: the control parameter for the output
!c  c**************************************************************************

	module EsField_mod
	contains
		subroutine EsField(ckaiD,cEs,ir,is,kinv,ifr,nfield)
		use model_mod
		use field_mod
		use geo_scatterer_mod
		
		implicit none
      
      
      
		complex, intent(in) :: ckaiD(mx,mz)	
		complex, intent(out) :: cEs
		complex :: cEs_matrix 
		integer, intent(in) :: ir, is, kinv, ifr, nfield
		integer :: i, k, irow, jcoL  
		complex :: Eyinc_pre    



		cEs=0.0

!c      using the new scattering approximation

		irow=ir+(is-1)*nrecey
		do k=1,mz
			do i=1,mx
				jcoL=i+(k-1)*mx     
				cEs=cEs+EyRinc(i,k,ir,ifr)*ckaiD(i,k)
			enddo
		enddo

		cEs=cEs*econst1


      end subroutine EsField
	end module EsField_mod
	
		module EiDipole_mod
	contains
		subroutine EiDipole(xxf,xzf,finc,is)

		use model_mod
		use weakform_mod
		use geo_background
		use green_mod
        use geo_background
        use general_constants
		use bessel_mod
		use green_mod
		use precision_mod, wp => dp
		implicit none


		integer :: i, k, is, Nrct
		logical :: Test_Green
		real :: x_cell(mx), z_cell(mz), xxs, xzs
		real :: xzb, xzf, xxf, xrho

		complex :: finc 


		!    variables for temp green's functions
		!complex :: Gxx_minusF
!--------------------------------------------------

		complex :: ca(nsegmax), cb(nsegmax), cG(2)

		!-------Variables for 2D weakform
		!complex*16 cjd0,chd0,kr
		complex (wp):: cjd0, chd0, kr    
		real ::  xr

		!c-----------------------------------------------




		Test_Green = .true.
		xzb = 0.0

		!c--------Note--------------------

		Gxx_minusF = cmplx(0.0, 0.0)

		xxs = Tloc(is,1)
		xzs = Tloc(is,2)

          

		!c----in the case of nobj=nobs, we need process the direct term separately
		!c-----------2D weak form--for direct term-----------

		xr = sqrt((xxf-xxs)**2+(xzf-xzs)**2)
		kr = xr*ck(nobj)

		call bes(kr,0,cjd0,chd0,0)
	
		Gxx_minusF=-(0.0,0.25)*(2.0*cjd0-chd0)

		xrho=abs(xxf-xxs)

		call InsegA(ca,cb,xrho,Nrct)
		call GreenOnePointA(xrho,xzf,xzs,xzb,ca,cb,cG,Nrct)

		finc=-(Gxx_minusF+cG(1)+cG(2))*ci*xomega*xmu0                      

          
      end subroutine EiDipole

	end module EiDipole_mod
	

	!  This part to calculate the incident field,i.e. scaled green's function

	module EiDipoleYBC_mod
	contains
		subroutine EiDipoleYBC(xxs,xzs,is,ifr,nfield)
		use field_mod
		use Ddomain_mod 
		use weakform_mod
		use model_mod
        use geo_background
        use general_constants
		use bessel_mod
		use green_mod
		use precision_mod, wp => dp
		implicit none

		integer :: i,k,Nrct
		logical :: Test_Green
		real, intent(in) :: xxs, xzs
		integer, intent(in) :: is, ifr, nfield
		real :: xzb, xzf, xxf, xrho



		!    variables for temp green's functions
		!complex :: Gxx_minusF
		!--------------------------------------------------

		complex :: ca(nsegmax),cb(nsegmax),cG(2)

		!-------Variables for 2D weakform
		complex(wp) :: cjd0, chd0, kr
	    
		real :: xr

		!-----------------------------------------------


		complex :: temp(mx,mz)

		Test_Green = .true.
		xzb = 0.0

		!--------Note--------------------




		Gxx_minusF=(0.0,0.0)

		do k=1,mz
			xzf=z_cell(k)
			do i=1,mx
				xxf=x_cell(i)

				!----in the case of nobj=nobs, we need process the direct term separately
				!-----------2D weak form--for direct term-----------

				xr=sqrt((xxf-xxs)**2+(xzf-xzs)**2)
				kr=xr*ck(nobj)

				if(nobj.eq.nobs) then
					if(xr <= 0.5*a0) then
						Gxx_minusF=xgm
					else
						call bes(kr,0,cjd0,chd0,0)
						Gxx_minusF=(2.0*cjd0-chd0)
					endif 
					
				else
					call bes(kr,0,cjd0,chd0,0)
					Gxx_minusF=(2.0*cjd0-chd0)
				endif

				Gxx_minusF=-(0.0,0.25)*Gxx_minusF

				xrho=abs(xxf-xxs)

				call InsegA(ca,cb,xrho,Nrct)

				call GreenOnePointA(xrho,xzs,xzf,xzb,ca,cb,cG,Nrct)

				EyRinc(i,k,is,ifr)=-(Gxx_minusF+cG(1)+cG(2))*ci*xomega*xmu0   
				Egrinc(i,k,ifr)=Gxx_minusF+cG(1)+cG(2)                   
			enddo
		enddo


		end subroutine EiDipoleYBC
	end module EiDipoleYBC_mod
	
	
	!-------------just for Hz component ------------------------------------
	module EsFieldhz_mod
	contains
		subroutine EsFieldhz(ckaiD,cEs,ir,is,kinv,ifr,nfield)

		use model_mod
		use field_mod
		use geo_scatterer_mod
		implicit none


		complex, intent(in) :: ckaiD(mx,mz)
		complex, intent(out) :: cEs
		integer, intent(in) :: ir, is, kinv, ifr, nfield
		complex :: cEs_matrix
		integer :: i,k,irow,jcoL
		complex :: Eyinc_pre


		cEs=0.0

		!c     this constant is put in the inputpara.f      
		!c      econst1=dx*dz*ci*xomega*cer(nobj)*xepsilon0
		!c      using the new scattering approximation



		do k=1,mz
			do i=1,mx
				cEs=cEs+EyRinchz(i,k,ir,ifr)*ckaiD(i,k)
			enddo
		enddo

		cEs=cEs*econst1

		end subroutine EsFieldhz

	end module EsFieldhz_mod
	
	
		!   This part to calculate the incident field,i.e. scaled green's function

		module EiDipoleYBChz_mod
		contains
			subroutine EiDipoleYBChz(xxs,xzs,is,ifr,nfield)
			use field_mod
			use Ddomain_mod 
			use model_mod
			use weakform_mod
            use general_constants
            use geo_background
			use bessel_mod
			use green_mod
			use precision_mod, wp => dp
			implicit none
			
			integer, intent(in) :: is, ifr, nfield
			integer :: i, k, Nrct
			logical :: Test_Green
			real, intent(in) :: xxs,xzs
			real :: xzb,xzf,xxf,xrho


			!    variables for temp green's functions
			!complex :: Gxx_minusF
			!--------------------------------------------------

			complex :: ca(nsegmax),cb(nsegmax),cG(2)

			!-------Variables for 2D weakform
			complex(wp) :: cjd0, chd0, kr      
			real :: xr

			!-----------------------------------------------


			complex temp(mx,mz)

			Test_Green=.true.
			xzb=0.0

			!--------Note--------------------




			Gxx_minusF=(0.0,0.0)

			do k=1,mz
				xzf=z_cell(k)
				do i=1,mx
					xxf=x_cell(i)

					!----in the case of nobj=nobs, we need process the direct term separately
					!-----------2D weak form--for direct term-----------

					xr=sqrt((xxf-xxs)**2+(xzf-xzs)**2)
					kr=xr*ck(nobj)

					if(nobj == nobs)then

						if(xr <= 0.5*a0)then
							Gxx_minusF=xgm
						else
							call bes(kr,0,cjd0,chd0,0)
							Gxx_minusF=(2.0*cjd0-chd0)
						endif 

					else
						call bes(kr,0,cjd0,chd0,0)
						Gxx_minusF=(2.0*cjd0-chd0)
					endif

					Gxx_minusF=-(0.0,0.25)*Gxx_minusF

					xrho=abs(xxf-xxs)

					call InsegA(ca,cb,xrho,Nrct)

					call GreenOnePointA(xrho,xzs,xzf,xzb,ca,cb,cG,Nrct)

					EyRinchz(i,k,is,ifr)=-(Gxx_minusF+cG(1)+cG(2))*ci*xomega*xmu0   
					Egrinc(i,k,ifr)=Gxx_minusF+cG(1)+cG(2)                   
				enddo
			enddo



			end subroutine EiDipoleYBChz
		end module EiDipoleYBChz_mod
	
	
	module Escat_mod
	contains
      subroutine Escat(ckaiD,kinv,is,ifr,nfield)
		use model_mod
		!use receivers_mod
		use field_mod
		use Ddomain_mod  
		 
		use field_mod
		!use receiver_mod
        use general_constants
		use EsField_mod
		use EiDipole_mod
		use EiDipoleYBC_mod
		use EiDipoleYBChz_mod
		use EsFieldhz_mod
		implicit none
      

		complex, intent (in) :: ckaiD(mx,mz)
		integer, intent (in) :: kinv,is,ifr,nfield
		complex :: cEs

		real :: rx,rz  !Rloc,
		complex :: cHx, cHz
        real :: rx1, rz1 !dx, dz,
		complex :: ces1, ces2, ces3, ces4
		complex :: finc, finc1, finc2, finc3, finc4
		complex :: cHcx, cHcy, cHsx, cHsy	  
		complex :: einc1(mx,mz), einc2(mx,mz)	  
		complex :: einc3(mx,mz), einc4(mx,mz)        
		integer :: i, k, kqr, m
      


!c      complex xsct
!c      common/Caldata/xsct(nrec,nsou)
!c*************************************************************************
!c calculate the coordinates of the cell center
!c*************************************************************************

!           
!c***********************************************************************
!c calculating the scattered Efield and write it to file 'Escat.dat'
!c***********************************************************************
!c for Ey component

		
		kqr=0        
		do m=1,nrey  
			kqr=kqr+1
			rx=Rlocey(m,1)
			rz=Rlocey(m,2)
			if(kinv == 0)then
				call EsField(ckaiD,cEs,m,is,kinv,ifr,nfield)
				call EiDipole(rx,rz,finc,is)
				write(333,'(2e30.14)')cEs
				write(334,*)real(finc),aimag(finc),real(ces),aimag(ces)
			endif
		enddo 


		if(nfield == 1)then
		! for Hx component       
			kqr=0
			do m=1,nrhx

			kqr=kqr+1         
			rx=Rlochx(m,1)         
			rz=Rlochx(m,2)       
			rz1=rz-dz

			call EiDipoleYBC(rx,rz1,m,ifr,nfield)

			do i=1,mx
				do k=1,mz
				einc2(i,k)=Egrinc(i,k,ifr)
				enddo
			enddo

			if(kinv == 0)then
				call EsField(ckaiD,ces2,m,is,kinv,ifr,nfield)
			endif
	   
			rz1=rz+dz

			call EiDipoleYBC(rx,rz1,m,ifr,nfield)

			do i=1,mx
				do k=1,mz
					einc4(i,k)=Egrinc(i,k,ifr)
				end do
			end do

			if(kinv == 0)then
				call EsField(ckaiD,ces4,m,is,kinv,ifr,nfield)
			endif


			if(kinv == 0)then

			cHsx=(0.0,1.0)/(xomega*xmu0)*(ces4-ces2)/dz/2.0


!  for output to inversion
			write(133,'(2e30.14)')cHsx
			endif


			do i=1,mx
				do k=1,mz
					Reinc_x(m,i,k,ifr)=(einc4(i,k)-einc2(i,k))/dz/2.0
					write(255,'(1x,2e30.14)') Reinc_x(m,i,k,ifr)
				end do
			end do

		enddo  
       

	 
       kqr=0
        do m=1,nrhz
			kqr=kqr+1
			rx=Rlochz(m,1)
			rz=Rlochz(m,2)


			rx1=rx-dx


			call EiDipoleYBChz(rx1,rz,m,ifr,nfield)

			do i=1,mx
				do k=1,mz
					einc1(i,k)=Egrinc(i,k,ifr)
				end do
			end do

			if(kinv == 0)then
				call EsFieldhz(ckaiD,ces1,m,is,kinv,ifr,nfield)
			endif


			rx1=rx+dx


			call EiDipoleYBChz(rx1,rz,m,ifr,nfield)

			do i=1,mx
				do k=1,mz
					einc3(i,k)=Egrinc(i,k,ifr)
				end do
			end do

			if(kinv == 0) then
				call EsFieldhz(ckaiD,ces3,m,is,kinv,ifr,nfield)
			endif


			if(kinv == 0) then
				cHsy=(0.0,1.0)/(xomega*xmu0)*(ces3-ces1)/dx/2.0
			!   for output to inversion
				write(144,'(2e30.14)')cHsy
			endif


			do i=1,mx
				do k=1,mz
					Reinc_z(m,i,k,ifr)=(einc3(i,k)-einc1(i,k))/dx/2.0
					write(266,'(1x,2e30.14)') Reinc_z(m,i,k,ifr)
				enddo
			enddo

		enddo
 
       
      endif


      if(kinv == 0)then
		write(*,*)'No. of Rx points: ',kqr
      end if 
   
      end subroutine Escat
	end module Escat_mod

	
