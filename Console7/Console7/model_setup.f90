	!c**************************************************************************
	!c Calculate the contrast and 1/epsilon matrixes associated with 
	!c the scattering object inside the computational domain.
	!c
	!c Input parameters:
	!c      x_cell, xy, z_cell: the the coordinates of the computational domain 
	!c                  in x, y and z directions and they are three vecotors.
	!c      dx, dy, dz:  increments in the x, y, z direction, respectively
	!c Output parameters:
	!c       ckai: the contrast matrix
	!c  
	!c**************************************************************************

	module kaieps_mod
	use geo_background
	use model_mod
	use general_constants
	use geo_scatterer_mod
	
	implicit none
    contains
		subroutine kaieps(x_cell,z_cell,dx,dz)        
			implicit none
			real, intent(in) :: x_cell(mx),z_cell(mz),dx,dz
			integer :: m,p,i0,k0,i1,k1,i,k
			real :: xt1,xt3,xto,xr
			complex :: cb1,temp1

			real :: dsize,dsize1
			complex :: ci0
			complex :: ci1
			!---local variables for the disk
			complex :: disk_cSer
			real :: disk_rad,xc1,zc1,xc2,zc2,disk_xSsig
			real :: xc3,zc3,xc4,zc4
			real :: tempr,tempi,sig
			integer :: percond_unit, contrast_unit, out_unit
			!c============================================================================
			!c the matrix for the constrast and 1/er(x,y,z)
			!c============================================================================
			!c
			!c Case sphere start !
			!c
			!C  output the model

			open(newunit = percond_unit, file='target_percond.dat', action='write')
			open(newunit = contrast_unit, file='target_contrast.dat', action='write')
			open(newunit = out_unit, file='target_out_f.dat', action='write')

			xto = xomega*xepsilon0

			do m=1,mx
				do p=1,mz
					cb1=scat_ceps-ci*sig/xto-cer(nobj)
					ckai(m,p)=cb1/(cer(nobj))
				enddo
			enddo






			

			close(percond_unit)
			close(contrast_unit)
			close(out_unit)
		end subroutine kaieps




	end module kaieps_mod
