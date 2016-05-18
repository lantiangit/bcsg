!c***************************************************************************
!c The 2D forward TM program in layered media
!c The green functions' part is from Ergun's part
!c 
!c subroutines:
!c
!c   input: read the input file and discretize the computational domain
!c
!c   kaieps: calculate the contrast matric and 1/epsilon matrix
!c
!c   GreenDdomain: Set computational domain and the source location
!c             for the Green's function then call subroutine to 
!c             calculate the Green's function in spatial domain
!c
!c   GreensforwardFFT: Do the FFT for the Green's function in Spatial domain
!c
!c   EiDipoleY: Calculate the incident electric field inside 
!c              the computational domain from a y-directed dipole
!c              in the top layer.
!c
!c    
!c   bcgs: Using Bi-CGSTAB to solve the linear operator eqation
!c
!c   kaiD: Calculate the contrast matrix multiply the electric flux density
!c         matrix
!c   Escat: Calculate the incident and scattered field at specified 
!c          receiver locations.
!c   kinv = 0 forward operation at the ground truth model，  the control parameter for the output，
!    kinv = 1 forward operation at the current model or updated model:
!c   produce incident field, Green's functions in the background, the data
!c	
!c   fre是频率数组，ifr是第几个频率,fre貌似没有什么用
!nfield  !表示只计算0=only E field,1=both E and H field,2=only H field
!inc_str 入射场的文件名
!c 
!c**************************************************************************
	module fwd2D_mod
		implicit none
		contains
		
		subroutine fwd2D(kinv,fre,ifr,nfield,inc_str)
			use model_mod
			use geo_background	
			use field_mod
            use Ddomain_mod
			use model_mod
			use green_mod
			use EiDipoleY_mod
			use read_mod
            use bcgs_mod
            use kaiD_mod
            use Escat_mod
            
            
			implicit none
			
			integer, intent(in) :: kinv, ifr, nfield 
			real, intent(in) :: fre    
			character(len = 20), intent(in) :: inc_str
			
			integer :: i,k, mx2,irow,jcoL
			real :: xx0, xz0 !source location		

			complex :: xn(mx,mz)
			integer :: i_ite,m,j,ns      !i_ite循环编号
			real :: px1, px2
			
			integer :: einc_unit, erecinc_unit, escat_unit, re_escat_unit, hxcat_unit, hzcat_unit
			open(newunit= einc_unit, file='einc.dat', action='write')
			open(newunit= erecinc_unit, file='erecinc.dat', action='write')
			open(newunit= escat_unit, file='escat.dat', action='write')
			open(newunit= re_escat_unit, file='re_escat.dat', action='write')
			open(newunit= hxcat_unit, file='hxcat.dat', action='write')
			open(newunit= hzcat_unit, file='hzcat.dat', action='write')
			

	  
!c**********************************************************************
!c Compute the Green's function in spatial domain                     *
!c**********************************************************************

!  initizlazing Green's function--------------

			print*, 'Compute the Greens function in D domain'
!    set up background Green's function
			call GreenDdomain(dx,dz)


!***********************************************************************
!  solve the operator equation by CG, BCG, BCGS and obtain the induced 
!  electric current inside the scatterer
!***********************************************************************

!compute the incident field and other constants                       * 



		mx2=mx+2

!----------begin to calculate the incident field and scattered field 
        
!       first step calculate the incident field (Green's function) for all sources
		if(kinv == 0) then            
			do  i_ite=1,NT
				xx0=Tloc(i_ite,1)
					xz0=Tloc(i_ite,2)

					print*,'Tx point: ', i_ite
					write(*,*)'xx0,xz0, nobs',xx0,xz0,nobs

					call EiDipoleY(x_cell,z_cell,xx0,xz0,i_ite,ifr)


					do i=1,mx
						do k=1,mz
							write(einc_unit,'(1x,2e30.14)') Eyinc(i,k,i_ite,ifr)
						enddo
					enddo
			enddo
		endif
          
		if(kinv == 0) then            
			do i_ite=1, nrey
				xx0=Rlocey(i_ite,1)
				xz0=Rlocey(i_ite,2)

				print*,'Rx point: ', i_ite
				write(*,*)'xx0,xz0, nobs',xx0,xz0,nobs
				call EiDipoleYB(x_cell,z_cell,xx0,xz0,i_ite,ifr,nfield)

				do i=1,mx
					do k=1,mz
						write(erecinc_unit,'(1x,2e30.14)') EyRinc(i,k,i_ite,ifr)
					enddo
				enddo
			enddo
		endif
      
          
		if(kinv == 2) then
			call einc_read(inc_str,nfield)
		endif

          
!    To solve the total field in the D domain for each source, given incident field.

		do  i_ite=1,NT
			xx0=Tloc(i_ite,1)
			xz0=Tloc(i_ite,2)
			if(kinv == 0)	write(escat_unit,*) i_ite, ifr
			if(kinv == 0)	write(re_escat_unit,*) i_ite, ifr

			if(nfield == 1)	write(hxcat_unit,*)	i_ite, ifr
			if(nfield == 1)	write(hzcat_unit,*)	i_ite, ifr
		


			if(kinv == 0)then
				do i=1,mx
					do k=1,mz
						cei(i,k)= Eyinc(i,k,i_ite,ifr)
					end do
				end do

				call bcgs(dx,dz,xn)
			endif


			call kaiD(xn) ! cd is current (contrast source)

	!c Now using reciprocity to calculate scattered field in the receiver

	!c       print*, 'Computing the scattered E-Field in Layer 1'
	!c        call Escat(xx0,xz0,dx,dz,x_cell,z_cell,cD,kinv,i_ite)

			call Escat(xn,kinv,i_ite,ifr,nfield)
		enddo
		
		close(einc_unit)
		close(erecinc_unit)
		close(escat_unit)
		close(re_escat_unit)
		close(hxcat_unit)
		close(hzcat_unit)
		
    end subroutine fwd2D
	end module fwd2D_mod






