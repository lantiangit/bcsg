	module precision_mod
		implicit none
		integer, parameter :: sp = selected_real_kind( 6, 37)
		integer, parameter :: dp = selected_real_kind(15, 307)
		integer, parameter :: qp = selected_real_kind(30, 291)
		
		integer, parameter :: i8 = selected_int_kind( 2)
		integer, parameter :: i16 = selected_int_kind( 4)
		integer, parameter :: i32 = selected_int_kind( 9)
		integer, parameter :: i64 = selected_int_kind(15)
		
		real (dp) , parameter :: tinydp =tiny (1.0_dp ) , tinyfactor = 5.0
	end module precision_mod


	module near0dp_mod

	contains
	elemental function Near0dp (testnumber, epsilon) result (return_value)

		use precision_mod
		implicit none
		real (dp) , intent (in) :: testnumber
		real (dp) , intent (in), optional :: epsilon
		logical :: return_value
		real (kind (epsilon) ) :: local_epsilon
		local_epsilon = tinyfactor * tinydp
		if (present (epsilon) ) then
			if (abs (epsilon)>=TINYDP) &
				local_epsilon =abs (epsilon)
		end if
		return_value =abs (testnumber )<local_epsilon

	end function Near0dp

	end module near0dp_mod
	
	
	module general_constants
	use precision_mod
		implicit none
		real, parameter :: xpi = acos(-1.0)
		real, parameter :: xmu0 = 4.0*xpi*1.0e-07
		real, parameter :: c0 = 2.99792458e8
		real, parameter :: xepsilon0 = 1.0/(c0**2*xmu0)
		real :: xomega
		complex(dp), parameter :: ci=cmplx(0,1) !?虚数单位j

	end module general_constants


	!计算区域的配置参数
	module model_mod
		implicit none
		!最大的层数，？，x方向的数量，z方向的数量，x方向cell的数量，z方向cell的数量，？
		integer, parameter :: nlayermax=8, nsegmax=2000000, mx=90, mz=60,nx=2*mx+2, nz=2*mz+2, mfr=6
		!y方向接收点数量，x方向接收点的数量，z方向接收点的数量
		integer, parameter :: nrecey=31, nrechx=31, nrechz=31
		!频率点的数量
		real, parameter :: mfreq=50
		!源的数量，接收点的数量
		integer, parameter :: nsou=3, nrec=31, mrow=nsou*nrec, nuk=mx*mz
		real :: Tloc(nsou,2)
		integer :: NT
		real :: Rloc(nrec,2)
		integer :: NR
		integer :: nrey,nrhx,nrhz
		real :: Rlocey(nrecey,2), Rlochx(nrechx,2), Rlochz(nrechz,2)
		
        
    end module model_mod
	

	
    
    module Ddomain_mod
		use model_mod
		implicit none
		real :: dx,dz,dxz
        real :: x_cell(mx),z_cell(mz)
	end module Ddomain_mod



	!c The matricizes for the Green's function spliting into *_plus and *_minus
	!c according to it is a function of  z+z' or z-z'
	!c  nkrhoGreen: the number of segmments for Sommerfeld integral evaluations
	!c              for the Green's function not used in the program

	!c*************************************************************************
	!c nlayer: the number of layers of the background
	!c nobj: the layer where the target resides
	!c nobs: the layer for the observation points, this parameter is used
	!c       for calculating Green's function when the source and observation
	!c       points are not in the same layer
	!c cer: a vector of the complex dielectric constant of each layer
	!c ck: a vector of the wave number of each layer
	!c xsig: a vector of the conductivity of each layer
	!c xc, xyc, zc: the center of the scattering object
	!c xloc:  a vector of interfaces locations 
	!c pec_nlayer: A logical variable, if the nth layer is PEC, its value 
	!c             is '.true.', otherwise, 'false'
	!c**************************************************************************
	module geo_background
		use model_mod
		implicit none
		complex :: cer(nlayermax) !每一层的相对介电常数
		complex :: cks(nlayermax), ck(nlayermax) !每一层的波数k，但是ks是什么不知道
		
		real :: k_max !最大的k值,且只有赋值
		real :: xsig(nlayermax), xer(nlayermax)  !每一层的电导率
		real :: xc, zc !散射体的中心位置
		real :: xloc(nlayermax-1) !每一层的分界面
		real :: xxsp,xzsp,xxfp,xzfp  
		real :: xfreq    !频率
		integer :: nlayer !背景的层数, number of layers
		integer :: nobj, nobs   !  are the layer for the test domain, measurement layer目标所在的层，观察点所在的层，the layer where D domain resides，
		character (len=32) :: sname, fname  
		
		logical :: pec_1layer, pec_nlayer !第一次和最后一层是不是PEC，是为1，不是为0
        real :: tolerans,tolerans2
        complex :: ckai(mx,mz)
        

		
	end module geo_background






	!  Parameters to control the iterative sovler and the output
	! istepmax: the maximu number of iterations
	! ctrmax: the error criterior to terminate the iteration
	module iter_control_mod
		implicit none
		integer :: istepmax  
		real :: ctrmax
	end module iter_control_mod

	!The parameters for the scattering object and size of D domain
	module geo_scatterer_mod
		implicit none
		integer :: nscatter
		real :: xa,za,scat_rad,scat_sig, scat_eps
		complex :: scat_ceps
		complex :: econst1 !integcoeff
    end module geo_scatterer_mod


    module field_mod
        use model_mod
        implicit none
        complex :: Eyinc(mx,mz,nsou,mfr)
        complex :: Reinc_x(nrechx,mx,mz,mfr),Reinc_z(nrechz,mx,mz,mfr)    
		complex :: cei(mx,mz)
		complex :: EyRinc(mx,mz,nrecey,mfr)
		complex :: Egrinc(mx,mz,mfr)
		complex :: EyRinchz(mx,mz,nrechz,mfr)
    end module field_mod
    
	!module field_mod
	!	use model_mod
	!	implicit none
	!	complex :: EyRinc(mx,mz,nrecey,mfr)
	!	complex :: Egrinc(mx,mz,mfr)
	!	complex :: EyRinchz(mx,mz,nrechz,mfr)
	!end module field_mod








	module weakform_mod
		implicit none
		complex :: ck01, xgm 
		real :: a0
	end module weakform_mod



	
	module para_mod
	  implicit none
      complex :: xk02
	end module para_mod


	!c***********************************************************************
	!c Deciding the layer (nobs) of the observation point
	!c***********************************************************************
	!决定观察点所在的层
	module layer_no_mod
	contains
		  subroutine LayerNo(z_cell,nobsk)
		  use model_mod
		  use geo_background
		  implicit none
		  real, intent(in) :: z_cell
		  integer, intent(inout) :: nobsk
		  integer :: index

		  if (z_cell > xloc(1)) then
			nobsk=1
		  else if (z_cell < xloc(nlayer-1)) then
			nobsk=nlayer
		  else
			do index=2, nlayer-1
			  if ( (z_cell < xloc(index-1) ) .and.  (z_cell > xloc(index) ) ) then
				nobsk=index
			  end if
			end do
		  
		  end if
		  
		end subroutine LayerNo
		
	end module layer_no_mod
	
	
	module green_global_parameters
	use model_mod
	implicit none
		complex :: Gxx_minusF, Gxx_plusF !临时的格林函数
		complex :: fgreenkm(Nx,Nz), fgreenkp(Nx,Nz)	
		complex :: fgreenkmb(Nx,Nz), fgreenkpb(Nx,Nz)
		complex :: Gxx_minus(mx+2,Nz), Gxx_minuss(mx+2,Nz) ,Gxx_plus(mx+2,Nz), Gxx_pluss(mx+2,Nz) 
	end module green_global_parameters

