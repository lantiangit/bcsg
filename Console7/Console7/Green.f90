



module green_mod
	  use model_mod
	  use green_global_parameters
      use general_constants
      use geo_background
      use weakform_mod
      use bessel_mod
      use xu_mod
      use precision_mod, wp => dp
	  implicit none
	  
	  
	  contains
!c************************************************************************
!c This program changes the sign of k_z,i.
!c The sign of k_z,i satisfy: real(k_z,i)<0 imag(k_z,i)>0
!c************************************************************************

	  subroutine process_ckz(ckz)
      implicit none
	  complex, intent(inout) ::ckz
      real :: xrkz, xikz
      complex :: ci

      ci = cmplx(0,1)
      xrkz=real(ckz)
      xikz=imag(ckz)
      if (xrkz < 0) then
         xrkz=-xrkz
      endif
      if (xikz > 0) then
         xikz=-xikz
      endif
      ckz=xrkz+ci*xikz
 
      end  subroutine process_ckz
	  
	  
	   complex function refl(cZ, cZ0) 
		  implicit none
		  complex, intent(in) :: cZ, cZ0
		  refl=(cZ-cZ0)/(cZ+cZ0)	  
	  end function refl  
	  
	  
	  subroutine InputImpedance(ckz,n,cZTE,cRTE,cpsi)
		  implicit none
		
		  complex, intent(in) :: ckz
		  integer, intent(in) :: n
		  complex, intent(inout) :: cZTE
		  complex, intent(out) :: cRTE, cpsi
		  complex :: cphi,ctant,cZTM,cZ0TE,cZ0TM

		  complex :: cRTM


		  cZ0TE=1.0/ckz
		  cRTE=refl(cZTE,cZ0TE)

		  cpsi=0
		  if (n > 1 .and. n < nlayer) then
			 cpsi=exp(-ci*ckz*(xloc(n-1)-xloc(n)))
			 cphi=cpsi**2
			 ctant=(1-cphi)/(1+cphi)
			 cZTE=cZ0TE*(cZTE+cZ0TE*ctant)/(cZ0TE+cZTE*ctant)

		  endif   
      end subroutine InputImpedance
	  
        !c*************************************************************************
!c The subroutine calculating the reflection coefficient at the left
!c interface of layer i
!c input parameters
!c    mode--TE or TM mode
!c    ckrho2--k_rho^2 =k^2-k_z^2
!c    ckzs--wave vector for the layer i
!c Output parameters      
!c    cR-leftTE--reflection coefficient at the left interface of layer i
!c               for TE mode
!c    cR_leftTM--reflection coefficient at the left interface of layer i 
!c               for TM mode
!c    cLTTEV---trasmission coefficient from target layer to obsevation layer
!c              for TE mode: voltage response to a current source.
!c    cLTTMV---trasmission coefficient from target layer to obsevation layer
!c              for TM mode: voltage response to a current source. 
!c***************************************************************************	  
	  
	  subroutine R_left(ckrho2,z_cell,ni,nm,cR_leftTE,cLTTEV)
	  !use geo_background
	
      !use general_constants
      implicit none
      
	  complex, intent(in) :: ckrho2
	  real, intent(in) :: z_cell
	  integer, intent(in) :: ni, nm
	  complex, intent(out) :: cR_leftTE, cLTTEV
      integer :: n, nUpperLimit
      
      complex :: cR_leftTM

      complex :: cZTE, cZTM, ctant, cphi

      complex :: ckz(1:nlayermax), cRTE(1:nlayermax-1), cRTM(1:nlayermax-1)

      complex :: cpsi(1:nlayermax-1)
      complex :: cLTTMI,cLTTMV

      complex :: fconst1,fconst2

!c*****************************************************************
!c The bottom layer is PEC. For TE mode, cZ=0; for TE mode, cZ=inf
!c The input impedance are calculated in different way according 
!c to the target is located at nlayer-1 or not.  
!c*****************************************************************   
      if (pec_nlayer) then
         if (ni /= nlayer-1) then
            ckz(nlayer-1)=csqrt(ck(nlayer-1)**2-ckrho2)
            call process_ckz(ckz(nlayer-1))
            cpsi(nlayer-1)=exp(-ci*ckz(nlayer-1)*(xloc(nlayer-2)-xloc(nlayer-1)))
            cphi=cpsi(nlayer-1)**2
            ctant=(1-cphi)/(1+cphi)
            cZTE=ctant/ckz(nlayer-1)
            cRTE(nlayer-1)=-1
!c            cZTM=ckz(nlayer-1)/cer(nlayer-1)*ctant
!c            cRTM(nlayer-1)=-1
            nUpperlimit=nlayer-2
            do n=nUpperLimit,ni,-1
               ckz(n)=csqrt(ck(n)**2-ckrho2)
               call process_ckz(ckz(n))

               call InputImpedance(ckz(n),n,cZTE,CRTE(n),cpsi(n))
            enddo
			
			
            cR_leftTE=cRTE(ni)
         else
            CR_leftTE=-1
         endif
      else
         ckz(nlayer)=csqrt(ck(nlayer)**2-ckrho2)
         call process_ckz(ckz(nlayer))
         cZTE=1.0/ckz(nlayer)
!c         cZTM=ckz(nlayer)/(cer(nlayer))
         nUpperLimit=nlayer-1

!c******************************************************************
!c Calculate up towards the left boundary of the target layer
!c for each dielectric layer.
!c******************************************************************
         do n=nUpperLimit,ni,-1
            ckz(n)=csqrt(ck(n)**2-ckrho2)
            call process_ckz(ckz(n))

            call InputImpedance(ckz(n),n,cZTE,CRTE(n),cpsi(n))
         enddo
         cR_leftTE=cRTE(ni)
!c         cR_leftTM=cRTM(ni)
      endif

      if (nm > ni) then
         fconst1=exp(-ci*ckz(nm)*(xloc(nm-1)-z_cell))
         if (nm == nlayer) then
            cLTTEV=fconst1
         else
            fconst2=exp(-ci*2*ckz(nm)*(z_cell-xloc(nm)))
            cLTTEV=fconst1/(1+cRTE(nm)*cpsi(nm)**2)*(1+cRTE(nm)*fconst2)

         endif

         do n=ni+1, nm-1

            cLTTEV=cLTTEV*(1+cRTE(n))*cpsi(n)/(1+cRTE(n)*cpsi(n)**2)

         enddo
      endif


      end subroutine R_left
	      !c************************************************************************
!c The subroutine calculating the reflection coefficient at the right
!c interface of layer i
!c input parameters
!c    mode-- TE or TM mode
!c    ckrho2--k_rho^2 =k^2-k_z^2
!c    ncal--the layer of which the reflection coef. (R_right) is calculated.
!c Output parameters      
!c    cR-rightTE--reflection coefficient at the right interface of layer ncal
!c       for TE mode
!c    cR-rightTM--reflection coefficient at the right interface of layer ncal
!c       for TM mode
!c***************************************************************************
	  
	  subroutine R_right(ckrho2,z_cell,ni,nm,cR_rightTE,cRTTEV)

	  !use model_mod
	  !use geo_background
      !use general_constants
      implicit none
      
		
	  complex, intent(in) :: ckrho2
	  real, intent(in) :: z_cell
	  integer, intent(in) :: ni, nm
	  complex, intent(out) :: cR_rightTE, cRTTEV
      integer :: n, nUpperLimit      
      complex ::  cR_rightTM
      complex :: cZTE,cZTM,ctant,cphi
      complex :: ckz(nlayermax), cRTE(1:nlayermax-1),cRTM(1:nlayermax-1)
      complex :: cpsi(2:nlayermax-1), cRTTMI, cRTTMV
      complex :: fconst1,fconst2

      if (pec_1layer) then
         if (ni /= 2) then
            ckz(2)=csqrt(ck(2)**2-ckrho2)
            call process_ckz(ckz(2))
            cpsi(2)=exp(-ci*ckz(2)*(xloc(1)-xloc(2)))
            cphi=cpsi(2)**2
            ctant=(1.0e0-cphi)/(1.0e0+cphi)
            cZTE=ctant/ckz(2)
            cZTM=ckz(2)/cer(2)*ctant
            cRTE(2)=-1.0e0
            cRTM(2)=-1.0e0
            nUpperlimit=3
            do n=nUpperLimit,ni
               ckz(n)=csqrt(ck(n)**2-ckrho2)
               call process_ckz(ckz(n))

            call InputImpedance(ckz(n),n,cZTE,cRTE(n),cpsi(n))
            enddo 
            cR_rightTE=cRTE(ni)
!c            cR_rightTM=cRTM(ni)
         else
            CR_rightTE=-1.0e0
!c            cR_rightTM=-1.0e0
         endif
      else

!c**********************************************************************
!c compute characteristic impedance for layer 1
!c**********************************************************************
      ckz(1)=csqrt(ck(1)**2-ckrho2)
      call process_ckz(ckz(1))
      cZTE=1.0/ckz(1)
!c      cZTM=ckz(1)/cer(1)

!c***************************************************************************
!c compute the input impednace at the right interface of layer i recursively
!c***************************************************************************
      do n=2,ni
         ckz(n)=csqrt(ck(n)**2-ckrho2)
         call process_ckz(ckz(n))
         call InputImpedance(ckz(n),n,cZTE,cRTE(n-1),cpsi(n))


      enddo 

      cR_rightTE=cRTE(ni-1)
       endif

      if (nm < ni) then 
         fconst1=exp(-ci*ckz(nm)*(z_cell-xloc(nm)))
         if (nm == 1) then
            cRTTEV=fconst1
         else
            fconst2=exp(-ci*2*ckz(nm)*(xloc(nm-1)-z_cell))
            cRTTEV= fconst1/(1+cRTE(nm-1)*cpsi(nm)**2)*(1+cRTE(nm-1)*fconst2)

         endif
         
         do n=nm,ni-2
            cRTTEV=cRTTEV*(1+cRTE(n))*cpsi(n+1)/(1+cRTE(n)*cpsi(n+1)**2)
         enddo
      endif   
      end subroutine R_right
	  
	  
	  !c*********************************************
!c the direct term
!c*********************************************
      complex function fun1(ckzi,xzcall,xzs)

      implicit none
      real, intent(in) :: xzcall,xzs
      complex, intent(in) :: ckzi
	  complex :: ci

      ci=cmplx(0,1)
      fun1=exp(-ci*ckzi*abs(xzcall-xzs))

      end function fun1

!c*******************************************************************
!c function the reflection from the left interface--the second term
!c*******************************************************************
      complex function fun2(ckzi,xzcall,xzs,xzb,ni)

	  !use geo_background
	  !use general_constants
		implicit none
      
      

      integer, intent(in) :: ni
      real, intent(in) :: xzcall,xzs,xzb
      complex :: ckzi

      fun2=exp(-ci*ckzi*((xzcall+xzs+xzb)-2.0*xloc(ni)))

      end function fun2

!c******************************************************************
!c function the reflection from the right interface--the third term
!c******************************************************************
      complex function fun3(ckzi,xzcall,xzs,xzb,ni)
	  !use geo_background
	  !use general_constants
      implicit none
      

      integer, intent(in) :: ni
      real, intent(in) :: xzcall, xzs,xzb
      complex, intent(in) :: ckzi

      fun3=exp(-ci*ckzi*(2.0*xloc(ni-1)-(xzb+xzcall+xzs))) 

      end function fun3

!c******************************************************************
!c function the reflection from the right interface--the third term
!c******************************************************************
      complex function fun4(ckzi,xzcall,xzs,ni)

	  !use geo_background
	  !use general_constants
      implicit none

      integer, intent(in) :: ni
      real, intent(in) :: xzcall, xzs
      complex, intent(in) ::  ckzi
	  complex :: cphi

      cphi=ckzi*(xloc(ni-1)-xloc(ni))
      fun4=exp(-2*ci*cphi+ci*ckzi*(xzcall-xzs))

      end function fun4
	  
!c******************************************************************
!c function the reflection from the right interface--the third term
!c******************************************************************
      complex function fun5(ckzi,xzcall,xzs,ni)
	  !use geo_background
	  !use general_constants

      implicit none
      integer, intent(in) :: ni
      real, intent(in) :: xzcall, xzs
      complex, intent(in) :: ckzi
	  complex :: cphi

      cphi=ckzi*(xloc(ni-1)-xloc(ni))
      fun5=exp(-2*ci*cphi-ci*ckzi*(xzcall-xzs))
 
      end function fun5
	  
	  
	  
	  
	     
    !c****************************************************************      
!c compute the Green's function in spectral domain for multi-layer 
!c background
!c****************************************************************

	  subroutine  multilayer(ckrho2,ckzi,z_cell,xzs,xzb,ni,nm,GTEVI_minus, GTEVI_plus)
     ! use general_constants
      !use geo_background

      implicit none
      
      complex, intent(in) :: ckrho2, ckzi
      real, intent(in) :: z_cell, xzs, xzb
      integer, intent(in) :: ni,nm
      complex, intent(out) :: GTEVI_minus,GTEVI_plus	
	  
	  real :: xzcall

      complex :: cR_leftTE,cR_leftTM,cR_rightTE,cR_rightTM
      complex :: cLTTEV,cLTTMV,cLTTMI,cRTTEV,cRTTMV,cRTTMI


      complex :: cRRTE,cRRTM,cD0TE,cD0TM,cfphi0

      complex :: ft1,ft2,ft3,ft4,ft5
      complex :: ckzm, GTMVI_plus,GTMVI_minus,direct_temp

      if (nm < ni) then 
         xzcall=xloc(ni-1)
      else if (nm > ni) then
         xzcall=xloc(ni)
      else
         xzcall=z_cell
      endif

       cLTTEV=1.0
       cRTTEV=1.0

       cLTTMI=1.0
       cLTTMV=1.0
       cRTTMI=1.0
       cRTTMV=1.0

!c=================================================
!c the scatterer lies in layer 1
!c cR_right=0, the fourth and fifth terms are zeros
!c=================================================
      if (ni == 1) then
         call R_left(ckrho2,z_cell,ni,nm,cR_leftTE,cLTTEV)

         ft2=fun2(ckzi,xzcall,xzs,xzb,ni)
         GTEVI_plus=cR_leftTE*ft2
         GTEVI_minus=0

!c------------------------------------

!c======================================================
!c the scatterer lies in layer N
!c cR_left=0, the fourth term is zero         
!c======================================================
      else if (ni == nlayer) then
         call R_right(ckrho2,z_cell,ni,nm,cR_rightTE,cRTTEV)

         ft3=fun3(ckzi,xzcall,xzs,xzb,ni)
         GTEVI_plus=cR_rightTE*ft3
         GTEVI_minus=0
      else
!c*********************************************************************
!c     the scatterer lies in the middle
!c*********************************************************************
         call R_right(ckrho2,z_cell,ni,nm,cR_rightTE,cRTTEV)


         call R_left(ckrho2,z_cell,ni,nm,cR_leftTE,cLTTEV)


         cfphi0=exp(-ci*2*ckzi*(xloc(ni-1)-xloc(ni)))

         cRRTE=cR_rightTE*cR_leftTE
         cD0TE=1-cRRTE*cfphi0

         ft2=fun2(ckzi,xzcall,xzs,xzb,ni)
         ft3=fun3(ckzi,xzcall,xzs,xzb,ni)
         ft4=fun4(ckzi,xzcall,xzs,ni)
         ft5=fun5(ckzi,xzcall,xzs,ni)

         GTEVI_minus=cRRTE*(ft4+ft5)/cD0TE

         GTEVI_plus=(cR_leftTE*ft2+cR_rightTE*ft3)/cD0TE

      endif

!c*************************************************************
!c the source and observation points are not in the same layer
!c Transmission coeffiecients: cLTTEV,cLTTMV,cLTTMI,cRTTEV,
!c cRTTMV,cRTTMI are used
!c*************************************************************
      if (ni /= nm) then
         ft1=fun1(ckzi,xzcall,xzs)
         direct_temp=exp(-ci*ckzi*abs(xzs-z_cell))
        GTEVI_minus=(ft1+GTEVI_minus)*cRTTEV*CLTTEV-direct_temp

         GTEVI_plus=GTEVI_plus*cRTTEV*CLTTEV

         ckzm=csqrt(ck(nm)**2-ckrho2)
         call process_ckz(ckzm)

      endif      
      end subroutine multilayer

  
    
    !c**********************************************************************
!c Calculate the Green's funciotn in spectral domain.
!c Input parameters
!c    ckrho2-k_rho^2
!c    z_cell: z coordinate for observation point 
!c    xzs:z coordinater for source point
!c    xzb: a shift in z for correlation calculation
!c    ni, the source layer
!c    nm, the observation layer      
!c output parameters
!c    GTEVI: Gxx in spectral domain
!c    GTMVI: Gzz in spectral domain
!c    QII: part of Gzx Gzy term in spectral domain
!c Note: each term is split intp plus term (z+z') and minus term (z-z')
!c       for future computation convenience

!c**********************************************************************

      subroutine SpectralA(ckrho2,z_cell,xzs,xzb,ni,nm, GTEVI_minus,GTEVI_plus)
      !use general_constants
      !use geo_background
      implicit none
      	  
	  complex, intent(in) :: ckrho2
	  real, intent(in) :: z_cell, xzs, xzb
      integer, intent(in) :: ni,nm
      complex, intent(out) :: GTEVI_minus, GTEVI_plus

      complex :: ckzi, cc1

!c**********************************************************
!c k_z,i=sqrt(k_i^2-k_rho^2), cc1=2jk_z,i is conatant for a krho
!c***********************************************************
      ckzi=csqrt(ck(ni)**2-ckrho2)
      call process_ckz(ckzi)
      cc1=2.0*ci*ckzi
!C     check the sign of cc1

      call multilayer(ckrho2,ckzi,z_cell,xzs,xzb,ni,nm,GTEVI_minus,GTEVI_plus)


       GTEVI_minus=GTEVI_minus/cc1
       GTEVI_plus=GTEVI_plus/cc1
	   
      end subroutine SpectralA
	  
	  


!--------------------------------------------------------

      subroutine InsegA(ca,cb,wrho,Nrecta)
	  !use model_mod
	  !use geo_background
	  !use general_constants
      implicit none
     
		
	  integer, intent(inout) :: Nrecta      
      complex, intent(inout) :: ca(nsegmax),cb(nsegmax)
	  real, intent(inout) :: wrho
	  
	  integer :: N,i,N1
      real :: xL1,dk,dk1,kpmax,dk2
      
      if (wrho < 1.0e4/xfreq) wrho=1.0e4/xfreq

      kpmax = 1.2e0*real(k_max)
      dk=0.1e0*xpi/wrho
      dk2=xpi/wrho
      if (dk > kpmax) then
           dk1=kpmax/2e0
           N1=1
      else
          N1=(1+int(kpmax/dk))
          dk1=0.5e0*kpmax/N1
      endif
      xL1=0.1e0*xpi/wrho           
      N=10001

      ca(1)=0.0e0
      if (N1 == 1) then
         cb(1)=dk1+ci*xL1
         cb(2)=kpmax
         ca(2)=cb(1)
         do i=2,N
           ca(i+1)=cb(i)
           cb(i+1)=cb(i)+dk
         end do
      else
         do i=1,N1
            cb(i)=i*(dk1+ci*xl1/N1)
            ca(i+1)=cb(i)
            cb(N1+i)=kpmax/2+ci*xL1+i*(dk1-ci*xl1/N1)
			ca(N1+i+1)=cb(N1+i)
         end do
         cb(2*N1+1)=real(cb(2*N1)+dk)
         do i=2*N1,N
           ca(i+1)=real(cb(i))
           cb(i+1)=real(cb(i)+dk2)
         end do
      end if
      Nrecta=N1+3
     
      end subroutine InsegA
	  

    
    !*************************************************
!     32 POINT GAUSS QUADRATURE
!*************************************************

      subroutine GQ32A (FCT,Y)
	  implicit none
	  complex, intent(in) :: FCT(32)
	  complex, intent(inout) :: Y
      Y=.35093050047350483E-2*(FCT(32)+FCT(1))
      Y=Y+.8137197365452835E-2*(FCT(31)+FCT(2))
      Y=Y+.12696032654631030E-1*(FCT(30)+FCT(3))
      Y=Y+.17136931456510717E-1*(FCT(29)+FCT(4))
      Y=Y+.21417949011113340E-1*(FCT(28)+FCT(5))
      Y=Y+.25499029631188088E-1*(FCT(27)+FCT(6))
      Y=Y+.29342046739267774E-1*(FCT(26)+FCT(7))
      Y=Y+.32911111388180923E-1*(FCT(25)+FCT(8))
      Y=Y+.36172897054424253E-1*(FCT(24)+FCT(9))
      Y=Y+.39096947893535153E-1*(FCT(23)+FCT(10))
      Y=Y+.41655962113473378E-1*(FCT(22)+FCT(11))
      Y=Y+.43826046502201906E-1*(FCT(21)+FCT(12))
      Y=Y+.45586939347881942E-1*(FCT(20)+FCT(13))
      Y=Y+.46922199540402283E-1*(FCT(19)+FCT(14))
      Y=Y+.47819360039637430E-1*(FCT(18)+FCT(15))
      Y=Y+.48270044257363900E-1*(FCT(17)+FCT(16))
      
      end subroutine GQ32A
	  

    
    !c***********************************************************************
!c integrate over a straight line from ca to cb (parallel to imaginary axis
!c or on real axis) using 32-point Gaussian quadrature. 
!c Input parameters:
!c     inx, iny, inz: the indexes for the point (x_cell, xy, z_cell) in the matrix
!c                    for Green's function
!c Ouput parameters:
!c     The computed components of the Green's function for point (x_cell,yy,zz)
!c     are written to the matrixes at loacation (inx, iny, inz)
!c***********************************************************************

      subroutine lineinteA(xtrho,z_cell,xzs,xzb,caK,cbK,cr1)
	  
	 ! use model_mod
   
      !use general_constants
      !use geo_background  
      implicit none
      
      real, intent(in) :: xtrho, z_cell, xzs, xzb
	  complex, intent(in) :: caK, cbK
	  complex, intent(inout) :: cr1(2)
	  
      integer :: Nf  
      
      integer :: m,j,kk,Nseg
      real :: xerr,xerr2,xerr4,xinvpi
      complex :: cr0(2),cf(2,32), ctf(32)
      complex :: GTEVI_minus,GTEVI_plus,ctx
      complex :: ckrho,ckrho2,ctr,clen,cnewa
      xinvpi=1.0/xpi
      clen=cbK-caK
!c******************************************************************************
!c Integrate from ca to cb to evaluate the Sommerfeld integrals using 32-point 
!c Graussian quadrature
!c******************************************************************************
      Nf=2
      do j=1,32
         ckrho=clen*XU(j)+caK
         ckrho2=ckrho**2
            call SpectralA(ckrho2,z_cell,xzs,xzb,nobj,nobs,GTEVI_minus,GTEVI_plus)
         ctx = cos(ckrho*xtrho)*xinvpi
         cf(2,j) = GTEVI_plus*ctx
         cf(1,j) = GTEVI_minus*ctx
      enddo

      do m=1,Nf
         do j=1,32
            ctf(j)=cf(m,j)
         enddo
         call GQ32A(ctf,ctr)
         cr0(m)=ctr*clen
      enddo
 
   
      Nseg=1
	  
      do
       Nseg=Nseg*2      
       clen=(cbK-caK)/real(Nseg)
          do m=1,Nf
             cr1(m)=0
          end do
          do kk=1,Nseg
             cnewa=caK+real(kk-1)*clen
             do j=1,32
                ckrho=clen*XU(j)+cnewa
                ckrho2=ckrho**2
                call SpectralA(ckrho2,z_cell,xzs,xzb,nobj,nobs,GTEVI_minus,GTEVI_plus)

             ctx = cos(ckrho*xtrho)*xinvpi
             cf(1,j) = GTEVI_minus*ctx
             cf(2,j) = GTEVI_plus*ctx
             enddo

             do m=1,Nf
                do j=1,32
                   ctf(j)=cf(m,j)
                enddo
			
                call GQ32A(ctf,ctr)
                cr1(m)=cr1(m)+ctr*clen
             enddo
		 
            enddo

          xerr2=0.0
          xerr4=0.0
          if (abs(cr1(1)) > (1.0+tolerans)*abs(cr0(1)) .or. abs(cr1(1)) < (1.0-tolerans)*abs(cr0(1))) xerr2=1.0
          if (abs(cr1(2)) > (1.0+tolerans)*abs(cr0(2)) .or. abs(cr1(2)) < (1.0-tolerans)*abs(cr0(2))) xerr4=1.0
          if (abs(cr1(1))+abs(cr1(2)) < 1.0e-10) then
              xerr2=0.0
              xerr4=0.0
          endif
          xerr=max(xerr2,xerr4)
          if (xerr == 1.0 .and. Nseg <= 2049) then
    !         print*, Nseg, abs(cr1(1))+abs(cr1(2))
             do m=1,Nf
                cr0(m)=cr1(m)
             enddo
             cycle
          else
              exit
          end if
      
      enddo
   
      end subroutine lineinteA
      	  

    
    !cc******************************************************************
!cc calculate Green's function  for a point (x_cell, xy, z_cell),
!c The source point is located at (xxs, xys, xzs)
!c*******************************************************************

      subroutine GreenOnePointA(xrho,z_cell,xzs,xzb,ca,cb,cG,Nrct)

	  !use model_mod
      !use geo_background
      implicit none
      
		
	  real, intent(in) :: xrho, z_cell, xzs, xzb
	  complex, intent(in) :: ca(nsegmax), cb(nsegmax)
	  complex, intent(inout) :: cG(2)
	  integer, intent(in) :: Nrct
      integer :: i, n, Nf, qN
      real :: xerr
      complex :: cr(2)


      qN=10000
      Nf=2
      do  i=1,Nf
         cG(i)=0.0e0
      enddo
	  

      do i=1,Nrct
         call lineinteA(xrho,z_cell,xzs,xzb,ca(i),cb(i),cr)
         do  n=1,Nf
            cG(n)=cG(n)+cr(n)
         enddo		 
	  enddo
	  
      xerr=1.0
      i=Nrct+1
      do while ((xerr == 1.0) .and. (i < qN)) 
        call lineinteA(xrho,z_cell,xzs,xzb,ca(i),cb(i),cr)
         xerr=0.0
         i=i+1
         if (abs(cr(1)) > tolerans2*abs(cG(1))) then 
			xerr=1.0		
		 endif
		 
         if (abs(cr(2)) > tolerans2*abs(cG(2))) then 
			xerr=1.0
		 endif
		 
         if (abs(cG(1))+abs(cr(2)) < 1.0d-16) then 
			xerr=0.0
		 endif

         cG(1)=cG(1)+cr(1)
         cG(2)=cG(2)+cr(2)
        
      enddo  
	  
     end subroutine GreenOnePointA	  
	  

    
    
    !C****************************************************************************
!c Compute the Green's function inside the scattering domain which is 
!c   embeded in a layered medium background. 
!c*****************************************************************************

      subroutine GreensFunction_CompA(Test_Green,x_cell,z_cell,xxs,xzs,xzb)
	   
	  !use model_mod
	  use weakform_mod
      use bessel_mod
      use precision_mod, wp => dp
      implicit none
     
      integer Nrect
      logical, intent(in) :: Test_Green
	  real, intent(in) ::  x_cell,z_cell
      real, intent(in) ::  xxs, xzs, xzb


!     local variables:
      integer :: t
      real :: xdxmin, xr, xconst1, xrho
      complex :: ca(nsegmax),cb(nsegmax),cG(2)
      complex :: cep, cxm, ck02



!c-------Variables for 2D weakform

      complex (wp) :: kr, cjd0, chd0  
!c-----------------------------------------------

       Gxx_minusF=0.0
       Gxx_plusF=0.0

!c-----------2D weak form--for direct term-----------

           xconst1=(x_cell-xxs)**2
           xr=sqrt(xconst1+(z_cell-xzs)**2)
          if(xr <= (0.5*a0) )then
			Gxx_minusF=xgm
           else
           kr=xr*ck01
           call bes(kr,0,cjd0,chd0,0)  
           
		   !cjd0=bessel_j0(kr)
		   !chd0=bessel_y0(kr)
           Gxx_minusF=2.0*cjd0-chd0
           endif

           Gxx_minusF=-(0.0,0.25)*Gxx_minusF
!c       end if
!c----------------------------------------------
         xrho=abs(x_cell-xxs)
         call InsegA(ca,cb,xrho,Nrect)
         call GreenOnePointA(xrho,z_cell,xzs,xzb,ca,cb,cG,Nrect)
!c----------------below for H_0^2
             Gxx_minusF=Gxx_minusF+cG(1)
              Gxx_plusF=Gxx_plusF+cG(2)			 
         
	   end subroutine GreensFunction_CompA


    
    !c******************************************************************
!c copy cFull(NFx,NFz) to cG(Ngx,Ngz)
!c*******************************************************************

      subroutine copyG(CFull,Nfx,Nfz,CG,Ngx,Ngz)

      implicit none
      integer, intent(in) :: Nfx, Nfz, Ngx, Ngz
      complex, intent(in) :: CFull(Nfx,Nfz)
	  complex, intent(out) :: cG(Ngx,Ngz)
	  integer :: i,k

      do k=1,Ngz
            do i=1,Ngx
               cG(i,k)=CFull(i,k)
            enddo
       enddo   
   
      end  subroutine copyG  
	  
	  
!c************************************************************************
!c Using symmetry property in x, y and z directions to fill out the full 
!c dimension Green's function.
!c************************************************************************
      subroutine FillG0(cG,cGFull)
	  !use model_mod
      implicit none      

      integer :: i, ii, k, kk
      complex, intent(in) :: cG(mx+2, mz+2)
	  complex, intent(out) :: cGFull(Nx,Nz)

       do k=1,Nz
		   if (k <= mz+2) then
			kk=k
		   else
			kk=nz+2-k
		   endif
		   
		   do i=1,Nx
			   if(i<= mx+2) then
				ii=i
			   else
				ii=nx+2-i
			   endif
		   CGFull(i,k) = cG(ii,kk)
		   enddo
       enddo       
        
      end subroutine FillG0
      

!c************************************************************************
!c Using symmetry property in x and y direction to fill out the full 
!c dimension Green's function.
!c************************************************************************
      subroutine FillG1(cG,N1,N3,cGFull,Nxx,Nzz)
	  !use model_mod
      implicit none

		integer, intent(in) :: N1, N3, Nxx, Nzz
        complex, intent(in) :: cG(N1,N3)
		complex, intent(out) :: cGFull(Nxx, Nzz)
		integer :: i, ii, k

        do k=1,Nz
			do i=1,Nx
				if(i < mx+2) then
					ii=i
				else
					ii=nx+2-i
				endif
			CGFull(i,k)=cG(ii,k)
			end do
        end do
       
       end subroutine FillG1
		

    
    !c*****************************************************************
!c Using symmetric property to fill the Green's function in spatial
!c domain  and  do the forward fft
!c***************************************************************
	

      subroutine GreenForwardFFT
		!use model_mod
        use ffts_mod
		implicit none
     

		integer :: i,j, k
		real :: xr, xt1, xt2
		real :: x_cell(Nx)
		complex :: czxt,czyt

		complex :: cGtmp(mx+2,Nz), copen(Nx,Nz)
	
		call FillG0(Gxx_minus,Copen)

        !c-----this is the old part
		call fftshift(Copen)
		call fw2dfft(Copen)
        !c-------this is the old part

      do k=1,nz
		  do i=1,nx
			fgreenkm(i,k)=Copen(i,k)
		  enddo
      enddo

      call copyG(Copen,Nx,Nz,Gxx_minus,mx+2,mz+2)


      call FillG1(Gxx_plus,mx+2,Nz,Copen,Nx,Nz)
        !c-----------this is the old part
      call fftshiftxy(Copen)
      call fw2dfft(Copen)
        !c------------this is the old part
        !c      call fdfft(Copen)

      do k=1,nz
		  do i=1,nx
			fgreenkp(i,k)=Copen(i,k)
		  enddo
      enddo

      call copyG(Copen,Nx,Nz,Gxx_plus,mx+2,Nz)   
      end subroutine GreenForwardFFT

    !c*****************************************************************
    !c Using symmetric property to fill the Green's function in spatial
    !c domain  and  do the forward fft
    !c***************************************************************
      subroutine GreenForwardFFTb
		!use model_mod
        use ffts_mod
		implicit none
		
		integer :: i, j, k
		real :: xr, xt1, xt2
		real :: x_cell(Nx)
		complex :: czxt,czyt

		complex :: cGtmp(mx+2,Nz), copen(Nx,Nz)
		

    !c Initialize FFT codes, already done in the main program
    !c
    !c      acx=1.0e0/float(nx)
    !c      acz=1.0e0/float(nz)
    !c      call cffti(nx,xwsave)
    !c      call cffti(nz,zwsave)

    !c
    !c Gxx- 
    !c
		call FillG0(Gxx_minuss,Copen)

    !c-----this is the old part
		call fftshift(Copen)
		call fw2dfft(Copen)
!c-------this is the old part

!c      call fdfft(Copen)

		do k=1,nz
			do i=1,nx
				fgreenkmb(i,k)=Copen(i,k)
			enddo
		enddo

      call copyG(Copen,Nx,Nz,Gxx_minuss,mx+2,mz+2)

!c
!cGxx + 
!c
      call FillG1(Gxx_pluss,mx+2,Nz,Copen,Nx,Nz)
!c-----------this is the old part
      call fftshiftxy(Copen)
      call fw2dfft(Copen)
!c------------this is the old part
!c      call fdfft(Copen)

      do k=1,nz
		  do i=1,nx
			fgreenkpb(i,k)=Copen(i,k)
		  end do
      end do

      call copyG(Copen,Nx,Nz,Gxx_pluss,mx+2,Nz)

      end subroutine GreenForwardFFTb
	  
	  

    
    
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The subroutine to obtain the complex conjugate of Green
!c   functions in spatial domain
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		subroutine conjuGreen
		!use model_mod
		implicit none


		integer :: i,j,k

		 do k=1,Nz
			 do i=1,mx+2
				Gxx_pluss(i,k)=conjg(Gxx_plus(i,k))
			 enddo
		 enddo

		 do k=1,mz+2
			 do i=1,mx+2
				Gxx_minuss(i,k)=conjg(Gxx_minus(i,k))
			 enddo
		 enddo

		end subroutine conjuGreen
	
    
!*****************************************************************************
! This subrooutine sets up the domain for Green's function
!  subrouine: GreensFunction_Comp--calculating the Green's function in 
!             spatial domain.
! 
! Input parameters: 
!     dx, dz: increment in x and z direction 
!      
!*****************************************************************************

		subroutine GreenDdomain(dx,dz)		  
		implicit none

		real, intent(in) :: dx,dz		
		character(len=32) :: ich1
		integer :: i, j, k, mx2
		logical :: Test_Green		 
		real :: xto
		real :: xxs, xzs, xf(Nx), zf(Nz)
		real :: lamda, xr, xzb
		!c-------Variables for 2D weakform
		complex (wp) :: cjd0, chd0, cjd1, chd1, ka
		complex :: ck02     
		integer :: greenk_unit, greenkb_unit
		open(newunit= greenk_unit, file='greenk.dat', action='write')
		open(newunit= greenkb_unit, file='greenkb.dat', action='write')

		  
	!c-----------------------------------------------
		do i=1,Nx
			xf(i)=real(i-1-Nx/2)*dx
		enddo

		do i=1,Nz
			zf(i)=real(i-1-Nz/2)*dz
		enddo

		xxs=0.0
		xzs=0.0


		Test_Green=.true.

		mx2=mx+2
		xzb=2.0*(zc+0.5*dz)


		!c-----weak form for small distance
		ck02=ck(nobj)**2
		ck01=ck(nobj)
		a0=sqrt(dx*dz/xpi)
		ka=ck(nobj)*a0

		call bes(ka,0,cjd0,chd0,0)
		call bes(ka,1,cjd1,chd1,0)		   
		   
		xgm=2.*xpi*a0*cjd0*(2.*cjd1-chd1)/ck01-4.*ci/ck02
		xgm=xgm/(xpi*a0*a0)
		!c--------------------------      

		!---- start to loop in the D domain                


		do j=1,mx2 
			do k=1,Nz
				call GreensFunction_CompA(Test_Green, xf(j),zf(k),xxs,xzs,xzb)
				Gxx_minus(j,k)= Gxx_minusF
				Gxx_plus(j,k) = Gxx_plusF		      
			enddo
		enddo
				

!c-----do forwrad FFT for the Green's function in D domain
		call conjuGreen
		call GreenForwardFFT
		call GreenForwardFFTb
		print*,'Saving Data of GREEN function'
		write(greenk_unit,*)'Gxx_plus'
		write(greenkb_unit,*)'Gxx_pluss'
		do i=1,mx+2
			do k=1,Nz
				write(greenk_unit,'(1x,2e32.14)') Gxx_plus(i,k)
				write(greenkb_unit,'(1x,2e32.14)') Gxx_pluss(i,k)
			enddo
		enddo

		write(greenk_unit,*)'Gxx_minus'
		write(greenkb_unit,*)'Gxx_minuss'
		do i=1,mx+2
		    do k=1,mz+2
		        write(greenk_unit,'(1x,2e32.14)') Gxx_minus(i,k)
		        write(greenkb_unit,'(1x,2e32.14)') Gxx_minuss(i,k)
		    enddo
		enddo

		print*,'Saving over!'
		close(greenk_unit)
		close(greenkb_unit)

		end subroutine GreenDdomain		 
		

	  
end module green_mod



	
