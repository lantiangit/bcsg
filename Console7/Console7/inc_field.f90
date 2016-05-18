!   This part to calculate the incident field,i.e. scaled green's function


	module EiDipoleY_mod
		use model_mod
		use field_mod
		use weakform_mod
		use green_mod
        use geo_background 
        use general_constants
        !use field_mod
        use bessel_mod
		use precision_mod, wp => dp
       
	contains
      subroutine EiDipoleY(x_cell,z_cell,xxs,xzs,is,ifr)

		implicit none

		integer :: i,k,Nrct
		logical :: Test_Green
		integer, intent(in) :: is, ifr
		real, intent(in) :: x_cell(mx), z_cell(mz),xxs,xzs
		real :: xzb, xzf, xxf, xrho

   

		!    variables for temp green's functions
		complex :: Gxx_minusF
		!--------------------------------------------------

		complex :: ca(nsegmax),cb(nsegmax),cG(2)

		!-------Variables for 2D weakform
		complex(wp) :: cjd0, chd0, kr   
		real :: xr
		!-----------------------------------------------


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
				!           kr=xr*ck01
				call bes(kr,0,cjd0,chd0,0)
				Gxx_minusF=-(0.0,0.25)*(2.0*cjd0-chd0)


				xrho=abs(xxf-xxs)

				call InsegA(ca,cb,xrho,Nrct)
				call GreenOnePointA(xrho,xzf,xzs,xzb,ca,cb,cG,Nrct)

				Eyinc(i,k,is,ifr)=-(Gxx_minusF+cG(1)+cG(2))*ci*xomega*xmu0                      
			enddo
		enddo

    end subroutine EiDipoleY



    
    
    
    
    !   This part to calculate the incident field,i.e. scaled green's function

      subroutine EiDipoleYB(x_cell,z_cell,xxs,xzs,is,ifr,nfield)
	
        
		implicit none


		integer i,k,is,Nrct,ifr,nfield
		logical Test_Green
		real x_cell(mx),z_cell(mz),xxs,xzs
		real xzb,xzf,xxf,xrho

      
      

		!    variables for temp green's functions
		complex :: Gxx_minusF
		!--------------------------------------------------

		complex ca(nsegmax),cb(nsegmax),cG(2)

		!-------Variables for 2D weakform
		complex(wp) :: cjd0,chd0,kr    
		real :: xr

		!-----------------------------------------------
   
       

		Test_Green=.true.
		xzb=0.0

		!--------Note--------------------

          
          


		Gxx_minusF=(0.0,0.0)

		do k=1,mz
			xzf=z_cell(k)
			do i=1,mx
				xxf=x_cell(i)

				!c----in the case of nobj=nobs, we need process the direct term separately
				!c-----------2D weak form--for direct term-----------

				xr=sqrt((xxf-xxs)**2+(xzf-xzs)**2)
				kr=xr*ck(nobj)
				!c           kr=xr*ck01
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


				EyRinc(i,k,is,ifr)=-(Gxx_minusF+cG(1)+cG(2))*ci*xomega*xmu0    
					  
						 
			enddo
		enddo
          

      end subroutine EiDipoleYB

	  
	  end module EiDipoleY_mod

