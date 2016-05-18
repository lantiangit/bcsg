	module input_mod
    contains
		subroutine input(x_cell,z_cell,dx,dz,model_str,fre)
		use model_mod
        use general_constants
        use geo_background
        use iter_control_mod
        use geo_scatterer_mod
        use green_mod
		implicit none

		real, intent(in) ::  fre
        real, intent(out) :: x_cell(mx),z_cell(mz),dx,dz

		character(len=20), intent(in) :: model_str
		integer :: i, start, last
		real :: ermax,xto

		!   variables for the test domain

		real :: ztopD, zbotD
		integer :: model_unit
		

		
		!ci=(0.0d0,1.0d0)
		xfreq=fre
		print*,'The operating frequency is:', xfreq
		print*, 'Input model file is ', model_str
		!open(unit=model_unit,file=model_str,status='unknown')
		open(newunit= model_unit, file= model_str, status='old', action='read')

		read(model_unit,*) 
		!nlayer, nobj, number of layers, the layer for the test domain, measurement layer
		read(model_unit,*) nlayer, nobj, nobs
		if (nlayer > nlayermax) then
			print*, 'ERROR: parameter nlayer larger than memory allocated'
			print*, 'Check dimen.h and change nlayermax'
		endif
		print*,'The number of layer is:', nlayer

		read(model_unit,*)
		print*,'The (er,sig,z_i) for each layer are:'
		do i=1,nlayer-1
			read(model_unit,*) xer(i), xsig(i), xloc(i)
			cer(i) = cmplx(xer(i), 0)			
			print *, i, xer(i), xsig(i), xloc(i)
        enddo
		
        read(model_unit,*) xer(nlayer), xsig(nlayer)
        cer(nlayer) = cmplx(xer(nlayer), 0)
        print*, nlayer, cer(nlayer), xsig(nlayer)
		

		read(model_unit,*)
		read(model_unit,*) pec_nlayer
		print*,'The nth layer boudary is pec? : ',pec_nlayer
		read(model_unit,*)
		read(model_unit,*) pec_1layer
		print*,'The first layer boudary is pec? : ',pec_1layer



		!c------Find the k_max based on the background------
		!c    to form a complex permittivity

		xomega=2.0*xpi*xfreq
		xto=xomega*xepsilon0

		do i=1,nlayer
			cer(i)=cer(i)-ci*xsig(i)/xto
		enddo

		!cks is squared wave vector
		do i=1,nlayer
			cks(i)=xomega**2*(cer(i)*xepsilon0*xmu0)
			ck(i)=xomega*sqrt(cer(i)*xepsilon0*xmu0)
			call process_ckz(ck(i))
		enddo
		!-------------------------------------------
		start=1
		last=nlayer
		if (pec_nlayer == .true.) last=nlayer-1
		if (pec_1layer == .true.) start=2
		ermax=abs(cer(start))
		do i=start+1,last
			if (abs(cer(i)) > ermax)	ermax=abs(cer(i))
		enddo
		k_max=xomega*(sqrt(xepsilon0*xmu0))*sqrt(ermax)



		!c     the parameter for BCGSTAB solver
		read(model_unit, *)
		read(model_unit, *)istepmax,ctrmax
		print*,'The maximum number of iterations is:', istepmax
		print*,'The maximum residual is:',ctrmax




		!c   read size of the D domain 
		read(model_unit,*)
		read(model_unit, *) xa,za

		!c************************************************
		!c read the parameters for a Austria Profile
		!c************************************************

		
		
		read(model_unit,*) 
		read(model_unit,*) scat_rad, scat_eps, scat_sig
		scat_ceps = cmplx(scat_eps, 0)
		print *, 'the radius of the scatter is :::', scat_rad
		print *, 'the relative epsilon of the scatter is :', scat_eps
		print *, 'the sigma of the scatter is :', scat_sig


		dx=xa/real(mx)
		dz=za/real(mz)


		read(model_unit,*)
		read(model_unit,*) xc, zc
		print *, 'center of the test domain:', xc,zc
		!c----the error tolerance for the Sommerfeld Integration
		read(model_unit,*)
		read(model_unit,*) tolerans, tolerans2
		!c----------------modified part by LPS------------


		!C--------------This is the center coordinates of each cell
		do i=1,mx
			x_cell(i)=real(i)*dx-0.5*dx-0.5*xa+xc
		enddo

		do i=1,mz
			z_cell(i)=real(i)*dz-0.5*dz-0.5*za+zc
		enddo

		!c---------top side of the D domain
		ztopD=real(mz)*dz-0.5*za+zc
		!c---------bottom side of the D domain
		zbotD=-real(mz)*dz+0.5*za+zc
		print *, 'topside of the test domain', ztopD
		print *, 'bottomside of the test domain', zbotD
		if(ztopD > xloc(nobj-1))then
			print *,'top side of inverse doamin' 
			print *,'is across the interface'
		endif

		print *, 'The dist.top side of D to the interface'
		print *, xloc(nobj-1)-ztopD
		econst1=dx*dz*ci*xomega*cer(nobj)*xepsilon0
        
        print *,'The layer for measured domain: ',nobs
		print *,'The layer for object domain: ',nobj

		print *,'Wavelength in measurement layer',2.0*xpi/abs(ck(nobs))
		print *,'Wavelength in test layer',2.0*xpi/abs(ck(nobj))

      end subroutine input	  
    end module input_mod
