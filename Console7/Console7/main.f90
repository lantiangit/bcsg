program ch1207
	
    use fwd2D_mod
	use input_mod
    use Ddomain_mod
    use geo_background
    use general_constants
    use kaieps_mod
	use read_mod
	use fwd2D_mod
	
	implicit none
	
	character(len=20) :: model_str, inc_str, scat_str  !inc_str 入射场的文件名
	
	character(len=20) :: greenk_str, greenkb_str, greenDS_str, cppx_str
	character(len=20) :: trans_str, receiv_str
    integer :: ntest
	integer :: nfield  !表示只计算0=only E field,1=both E and H field,2=only H field
    integer :: itmax
	real, allocatable :: vfre(:) !频率数组
	integer :: nfr !频率数量
	real :: fre !频率临时变量
	integer :: i_fre !频率迭代指数
    integer :: i !迭代指数
  
    
    integer :: switch_unit, fre_unit	
	open(newunit= switch_unit, file='switch', status='old', action='read')	
	read(switch_unit,*)
	read(switch_unit,*) model_str
	read(switch_unit,*) scat_str
    read(switch_unit,*) inc_str
	read(switch_unit,*) greenk_str
	read(switch_unit,*) greenkb_str
	read(switch_unit,*) trans_str
	read(switch_unit,*) receiv_str
	read(switch_unit,*) nfield
	read(switch_unit,*) itmax
	print *,'Running a forward solver: yes=0,no=others'
	close(switch_unit)
    
    
    open(newunit = fre_unit, file='Frequency.inp', status='old', action='read')
    read(fre_unit,*)
    read(fre_unit,*) nfr
    read(fre_unit,*)
	allocate(vfre(1:nfr))
    do i_fre=1,nfr
        read(fre_unit,*) vfre(i_fre)
	enddo
	close(fre_unit)  
       
    do i_fre=1,nfr
        fre = vfre(i_fre)
		call input(x_cell,z_cell,dx,dz,model_str,fre) 
		call kaieps(x_cell,z_cell,dx,dz)
		call transreceiv(trans_str,receiv_str)	
		call fwd2D(0, fre, i_fre, nfield, inc_str)	
    enddo
	

end program ch1207
