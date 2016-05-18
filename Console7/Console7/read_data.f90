!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  The subroutine to read the scattered data
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
	module read_mod
	contains 
		subroutine escat_read(scat_str,snorm,nfield)
		use model_mod

		implicit none




!       passed variables:
        character(len=20), intent(in) :: scat_str
		integer, intent(in) :: nfield
		real, intent(in) :: snorm(mfr)
		
		real :: snormt, snormt1, snormt2

		real :: snormey(mfr), snormhx(mfr), snormhz(mfr)
		complex :: escat_obs(nrecey,nsou,mfr)
		complex :: hxcat_obs(nrechx,nsou,mfr), hzcat_obs(nrechz,nsou,mfr)
		!common/scatobserved/escat_obs,hxcat_obs,hzcat_obs
		!common/scatsnorm/snormey,snormhx,snormhz

!       local variable:
		integer :: is, ir, index_source, j, index_fre
		real :: qr, qi

		open(1,file=scat_str,status='old')

       
		if(nfield == 1)then
			open(2,file='hxcat.dat',status='old')
			open(3,file='hzcat.dat',status='old')
		endif
       
		write(*,*)'Reading the scattered data: ', scat_str
        
        
        do j=1,mfr
			snormt=0.0
			snormt1=0.0
			snormt2=0.0
			do is=1,nsou
			read(1,*)index_source,index_fre
        
			if(nfield == 1)then
				read(2,*)index_source,index_fre
				read(3,*)index_source,index_fre
			endif

			do ir=1,nrecey
				read(1,'(2e30.14)')qr,qi

		!----------new version-----------------------
				escat_obs(ir,is,j)=cmplx(qr,qi)
				snormt=snormt+(abs(escat_obs(ir,is,j)))**2
			enddo

			if(nfield == 1)then
			
				do ir=1,nrechx
					read(2,'(2e30.14)')qr,qi

			!----------new version-----------------------
					hxcat_obs(ir,is,j)=cmplx(qr,qi)
					snormt1=snormt1+(abs(hxcat_obs(ir,is,j)))**2
				enddo
				
				
				do ir=1,nrechz
					read(3,'(2e30.14)')qr,qi

			!----------new version-----------------------
				    hzcat_obs(ir,is,j)=cmplx(qr,qi)
					snormt2=snormt2+(abs(hzcat_obs(ir,is,j)))**2
				enddo
			
			 endif
			

			enddo
				snormey(j)=snormt
		   
			   if(nfield == 1)then
				   snormhx(j)=snormt1
				   snormhz(j)=snormt2
			   endif
       
		enddo

        close(1)
        if(nfield == 1)then
			close(2)
			close(3)
        endif
		
        end subroutine escat_read

	

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  The subroutine is to read incident field in D domain.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

         subroutine einc_read(inc_str,nfield)
		 use model_mod
		 use field_mod
		 !use field_mod
		  
		 implicit none


		 character(len=20), intent(in) :: inc_str
		 integer, intent(in) :: nfield


        

		 complex :: Reinc(nrecey,mx,mz,mfr)
      	 !common/recevincident/Reinc(nrecey,mx,mz,mfr)

      	
      	

!c       local variables
		integer :: ns,i,k,j 
		real :: qr1,qi1
		complex :: temp(mx,mz)
		


		open(3,file=inc_str,status='old')
		!c      open(3,file='einc_ys.dat',status='old')
		write(*,*)'Reading incident file: ',inc_str
	
        do j=1,mfr
			do ns=1,nsou
				do i=1,mx
					do k=1,mz
						read(3,'(1x,2e30.14)')qr1,qi1
						Eyinc(i,k,ns,j)=cmplx(qr1,qi1)
					enddo
				enddo
	        enddo
        enddo
        close(3)
        
        open(3,file='erecinc.dat',status='old')
!c        open(3,file='erecinc_ys.dat',status='old')
		write(*,*)'Reading incident file: erecinc.dat'
        
        
		do j=1,mfr
			do ns=1,nrec

				do i=1,mx
					do k=1,mz
						read(3,'(1x,2e30.14)')qr1,qi1
						temp(i,k)=cmplx(qr1,qi1)
						Reinc(ns,i,k,j)=cmplx(qr1,qi1)
						EyRinc(i,k,ns,j)=cmplx(qr1,qi1)
					enddo
				enddo
			enddo 
		enddo
		close(3)

        
        end subroutine einc_read
        




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     The subroutine is to read the  Green 's function in k-space.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine greenk_read(greenk_str,greenkb_str)
	  use model_mod
      use green_mod
	  implicit none




	  character(len=20) :: greenk_str, greenkb_str		

!     local variables
	   real :: qr,qi
	   integer :: i,k,j

      complex :: greenkm(Nx,Nz,mfr),greenkmb(Nx,Nz,mfr),greenkp(Nx,Nz,mfr),greenkpb(Nx,Nz,mfr)
      !common/greenf/greenkm(Nx,Nz,mfr),greenkp(Nx,Nz,mfr)
      !common/greenfb/greenkmb(Nx,Nz,mfr),greenkpb(Nx,Nz,mfr)
      
      complex :: tgreenkm(Nx,Nz),tgreenkp(Nx,Nz)
      complex :: tgreenkmb(Nx,Nz),tgreenkpb(Nx,Nz)


        open(3,file=greenk_str,status='old')
        open(4,file=greenkb_str,status='old')
        
   

        write(*,*)'Reading Green function in k-space: ', greenk_str
        write(*,*)'Reading Green function in k-space: ', greenkb_str


         do j=1,mfr
			 read(3,*)
			 read(4,*)
			 do i=1,mx+2
				 do k=1,Nz
				 read(3,'(1x,2e32.14)')qr,qi 
				 Gxx_plus(i,k)=cmplx(qr,qi)
				 read(4,'(1x,2e32.14)')qr,qi 
				 Gxx_pluss(i,k)=cmplx(qr,qi)

				 enddo
             enddo
      


         
			 read(3,*)
			 read(4,*)
			 do i=1,mx+2
				 do k=1,mz+2
					 read(3,'(1x,2e32.14)')qr,qi 
					 Gxx_minus(i,k)=cmplx(qr,qi)
					 read(4,'(1x,2e32.14)')qr,qi 
					 Gxx_minuss(i,k)=cmplx(qr,qi)
				 enddo
			 enddo
         

			call FillG0(Gxx_minus,tgreenkm)
			call FillG0(Gxx_minuss,tgreenkmb)
			call FillG1(Gxx_plus,mx+2,Nz,tgreenkp,Nx,Nz)
			call FillG1(Gxx_pluss,mx+2,Nz,tgreenkpb,Nx,Nz)

			greenkm(:,:,j)=tgreenkm
			greenkmb(:,:,j)=tgreenkmb
			greenkp(:,:,j)=tgreenkp
			greenkpb(:,:,j)=tgreenkpb
        enddo
        
         close(3)
         close(4)

		end subroutine greenk_read

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    The subroutine to read the transmitter and the 
!c     receiver coordinates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine transreceiv(trans_str,receiv_str)
		use model_mod
		implicit none	
        
		
!      passed variables:
        character(len=20), intent(in) :: trans_str, receiv_str
		integer :: i
       

!       NT: the number of the transmitters
!       NR: the number of the receivers
	
  
        open(1,file=trans_str,status='old')
        write(*,*)'Reading transmitter coordinates: ',trans_str
        read(1,*)
        read(1,*)NT
        read(1,*)
        do i=1,NT
			read(1,*) Tloc(i,1), Tloc(i,2)
			Tloc(i,1)=Tloc(i,1)*350
			Tloc(i,2)=749
        enddo
        close(1)

        open(1,file='Rlocey.inp',status='old')
		write(*,*)'Reading receiver coordinates: ',receiv_str
        read(1,*)
        read(1,*) nrey
        read(1,*)
        do i=1,nrey
			read(1,*)Rlocey(i,1),Rlocey(i,2)
			Rlocey(i,1)=Rlocey(i,1)*250
			Rlocey(i,2)=749
        enddo
        
        open(1,file='Rlochx.inp',status='old')
		write(*,*)'Reading receiver coordinates: ',receiv_str
        read(1,*)
        read(1,*)nrhx
        read(1,*)
        do i=1,nrhx
			read(1,*)Rlochx(i,1),Rlochx(i,2)
			Rlochx(i,1)=Rlochx(i,1)*250
			Rlochx(i,2)=749
        enddo
        close(1)
        
        open(1,file='Rlochz.inp',status='old')
		write(*,*)'Reading receiver coordinates: ',receiv_str
        read(1,*)
        read(1,*) nrhz
        read(1,*)
        do i=1,nrhz
			read(1,*) Rlochz(i,1), Rlochz(i,2)
			Rlochz(i,1)=Rlochz(i,1)*250
			Rlochz(i,2)=749
        enddo
        close(1)
        


        return      
        end subroutine transreceiv
	
		end module read_mod



