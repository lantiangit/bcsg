    module bessel_mod
    contains
    subroutine bslsdr(n,x,cj,ch,cjp,chp,m)
	!parameter(nmax=36)
    use precision_mod, wp => dp
    !implicit none
    integer, parameter :: nmax=36
    integer :: n,i,m
	complex (wp) :: cjd,chd,xd, jpd,hpd
	complex (wp) :: jd(-nmax:nmax),hd(-nmax:nmax)
	complex (wp) :: cj(-m:m), cjp(-m:m),ci
	complex (wp) :: ch(-m:m), chp(-m:m),x
	real (wp) :: rn


	ci = cmplx(0.0,1.0)
	do i=-n,n
	   cj(i) = 0.0
	   cjp(i) = 0.0
	   ch(i) = 0.0
	   chp(i) = 0.0
	end do
	xd = x
	call bes(xd,0,cjd,chd,0)
	jd(0) = cjd
	hd(0) = chd
	
	ir=1
	do i=1,n+1
	   call bes(xd,i,cjd,chd,0)
	   jd(i) = cjd
	   hd(i) = chd
	   jd(-i) = jd(i)
	   hd(-i) = hd(i)
	   if(ir > 0) jd(-i) = -jd(i)
	   if(ir > 0) hd(-i) = -hd(i)
	   ir = -ir
	enddo
	do i=-n,n
	   rn = real(i)
	   jpd = rn*jd(i)-xd*jd(i+1)
	   hpd = rn*hd(i)-xd*hd(i+1)
	   cj(i) = jd(i)
	   ch(i) = hd(i)
	   cjp(i) = jpd
	   chp(i) = hpd
    enddo
    
	
	end subroutine bslsdr
!*******************************************************************
! This is the code to calculate the Bessel Function of any order   *
!*******************************************************************
	 subroutine bes(x,i,bej,beh,iex)
!...   x: the input variable.
!...   i: order.
!... bej: the first kind Bessel Function J0.
!... beh: the third kind Bessel Function H0.Hankel function
!... iex: the exponential index (exp(ix*iex). Suggested value: 0   
!---------------------------------------------------------------     
    use precision_mod, wp => dp
    real (wp) :: pi,gr,gi,xi,ansr,ansi,anskr,anski !8 bytes means double precision
    complex (wp) :: ci,g,x,bej,beh                !16 bytes means both the real and 
                                                  !the imaginary parts are all in 
                                                  !double precision
	 common /bexp/iexpck                        
	 ci=(0.d0,1.d0)
	 pi=acos(-1.0)
	 iexpck=iex
	 xi=dimag(x)
	 g=x
	 if(xi.lt.0.d0)g=dconjg(x)
       	 g=-ci*g
	 gr=dreal(g)
	 gi=dimag(g)
  	 call bessel(gr,gi,i,ansr,ansi,anskr,anski)
	 bej=ci**i*dcmplx(ansr,ansi)
	 beh=2.d0*(-ci)**(i+1)*dcmplx(anskr,anski)/pi
	 if(xi.ge.0.d0)return
	 if(iex.eq.1)beh=beh*dexp(-2.d0*gr)
	 beh=2.d0*bej-beh
	 bej=dconjg(bej)
	 beh=dconjg(beh)
!-----------------------------------------------------                    

	 end subroutine bes
!*********************************************************
         subroutine bessel (zr,zi,n,ansr,ansi,anskr,anski)
         implicit double precision (a-h,o-z)
!---------------------------------------------------------                    
         common /bexp/ iexpck
         data pi/3.14159265358979323846264338d0/
         data gam/.5772156649015328606d0/
         sr=dsqrt(zr*zr+zi*zi)
         if (sr-(14.+n/2.))300,310,310
  300    continue
         ffn=n
         zdr=dabs(zr)
         zdi=dabs(zi)
         iq=1
         if (zi) 121,122,122
  122    if (zr.le.0.) iq=2
         go to 123
  121    iq=3
         if (zr.ge.0.) iq=4
  123    nm=n-1
         cer=zdr*.5d0
         cei=zdi*.5d0
         wsumr=0.
         wsumi=0.
         wwr=0.
         wwi=0.
         ucr=1.d0
         uci=0.
         if (n.eq.0) go to 210
         ucr=cer
         uci=cei
         if (n.eq.1) go to 210
         br=ucr
         bi=uci
         do 201 i=1,nm
            temp1=ucr*br-uci*bi
            uci=ucr*bi+uci*br
            ucr=temp1
  201    continue
  210    on=1.d0
         if (n.ne.0) on=(-1.d0)**n
         ucir=ucr
         ucii=uci
         ucr=ucr*on
         uci=uci*on
         vcr=ucr*.5d0
         vci=uci*.5d0
         pn=0.
         fn=1.d0
         if (n.eq.0) go to 10
         do 11 i=1,n
            fi=i
            pn=pn+1.d0/fi
   11    fn=fn*fi
   10    uur=1.d0/fn
         uui=0.
         vvr=pn/fn
         vvi=0.
         usumr=uur
         usumi=0.
         vsumr=vvr
         vsumi=0.
         temp1=cer*cer
         temp2=cei*cei
         cesr=temp1-temp2
         cesi=2.d0*cer*cei
         if (n.eq.0) go to 100
         temp1=1.d0/(temp1+temp2)
         wcr=cer*temp1
         wci=-cei*temp1
         if (nm.eq.0) go to 221
         br=wcr
         bi=wci
         do 220 i=1,nm
            temp1=wcr*br-wci*bi
            wci=wcr*bi+wci*br
            wcr=temp1
  220    continue
  221    wcr=wcr*.5d0
         wci=wci*.5d0
         fnm=fn/ffn
         wsumr=fnm
         wsumi=0.
         if (nm.eq.0) go to 35
         wwr=fnm
         wwi=0.
         do 30 k=1,nm
            fk=k
            ft=1.d0/(fk*(ffn-fk))
            wwr=-wwr*ft
            wwi=-wwi*ft
            temp1=wwr*cesr-wwi*cesi
            wwi=wwr*cesi+wwi*cesr
            wwr=temp1
            wsumr=wsumr+wwr
            wsumi=wsumi+wwi
   30    continue
   35    temp1=wsumr*wcr-wsumi*wci
         wsumi=wsumr*wci+wsumi*wcr
         wsumr=temp1
  100    pk=0.
         do 20 k=1,33
            fk=k
            fnk=fk+ffn
            d=fk*fnk
            tr=cesr/d
            ti=cesi/d
            temp1=uur*tr-uui*ti
            uui=uur*ti+uui*tr
            uur=temp1
            usumr=usumr+uur
            usumi=usumi+uui
            pk=pk+1.d0/fk
            pn=pn+1.d0/fnk
            cr=pn+pk
            vvr=uur*cr
            vvi=uui*cr
            vsumr=vsumr+vvr
            vsumi=vsumi+vvi
            temp1=usumr
            if (temp1.eq.0.) temp1=1.d0
            temp2=usumi
            if (temp2.eq.0.) temp2=1.d0
            temp1=uur/temp1
            temp2=uui/temp2
            if (dabs(temp1).ge..1e-15) go to 20
            if (dabs(temp2).ge..1e-15) go to 20
            temp1=vsumr
            if (temp1.eq.0.) temp1=1.d0
            temp2=vsumi
            if (temp2.eq.0.) temp2=1.d0
            temp1=vvr/temp1
            temp2=vvi/temp2
            if (dabs(temp1).ge..1e-15) go to 20
            if (dabs(temp2).ge..1e-15) go to 20
            go to 15
   20    continue
   15    temp1=vsumr*vcr-vsumi*vci
         vsumi=vsumr*vci+vsumi*vcr
         vsumr=temp1
         celr=.5d0*dlog(cer*cer+cei*cei)
         celi=datan2(cei,cer)
         ar=gam+celr
         ai=celi
         temp1=ar*usumr-ai*usumi
         ai=ar*usumi+ai*usumr
         ar=ucr*temp1-uci*ai
         ai=ucr*ai+uci*temp1
         anskr=wsumr-ar+vsumr
         anski=wsumi-ai+vsumi
         if (sr.ge.5.+n/4.) call zk24 (zdr,zdi,n,anskr,anski)
         ansr=usumr*ucir-usumi*ucii
         ansi=usumr*ucii+usumi*ucir
         if (iq.eq.1) go to 500
         if (iq.eq.4) go to 60
         go to 70
   60    ansi=-ansi
         anski=-anski
         go to 500
   70    cesr=1+2*(2*(n/2)-n)
         if (iq.eq.3) go to 80
         usumr=ansr
         usumi=-ansi
         br=-pi*usumi
         bi=pi*usumr
         anski=-anski
         anskr=cesr*anskr-br
         anski=cesr*anski-bi
         ansr=usumr*cesr
         ansi=usumi*cesr
         go to 500
   80    cr=-pi*ansi
         ci=pi*ansr
         anskr=anskr*cesr+cr
         anski=anski*cesr+ci
         ansr=ansr*cesr
         ansi=ansi*cesr
  500    if (iexpck.eq.0) return
         zrxp=dexp(zr)
         anskr=anskr*zrxp
       	 anski=anski*zrxp
         ansr=ansr/zrxp
         ansi=ansi/zrxp
         return
  310    call iln (zr,zi,n,ansr,ansi,anskr,anski)
!-------------------------------------------------------                    
       
         end subroutine bessel
!*******************************************************
         subroutine iln (zr,zi,n,ansr,ansi,anskr,anski)
         implicit double precision (a-h,o-z)
!-------------------------------------------------------                    
         common /bexp/ iexpck
         data pi/3.14159265358979323846264338d0/
         argr=2.d0*pi*zr
         argi=2.d0*pi*zi
!        complex square root
         if (argr.ge.0.) go to 101
         if (dabs(argi/argr).lt.4.d-7) go to 102
  101    pr=dsqrt((dsqrt(argr*argr+argi*argi)+argr)*.5d0)
         pii=argi/(2.d0*pr)
         go to 112
  102    axi=dsqrt(dabs(argr))
         if (argr) 103,104,104
  103    pr=0.
         pii=axi
         if (argi.lt.0.) pii=-pii
         go to 112
  104    pr=axi
         pii=0.
  112    fn=n
         fn2=4.d0*(fn**2)
         ar=1.d0
         if (iexpck.eq.0) ar=dexp(zr)
         br=dcos(zi)
         bi=dsin(zi)
         ai=1.d0/(pr*pr+pii*pii)
         cor=ar*(br*pr+bi*pii)*ai
         coi=ar*(-br*pii+bi*pr)*ai
         c2r=(br*pr-bi*pii)*ai/ar
         c2i=(-br*pii-bi*pr)*ai/ar
         ai=1+2*(2*(n/2)-n)
         if (zi.lt.0.) ai=-ai
         co2r=-ai*c2i
         co2i=ai*c2r
         if (iexpck.eq.0) go to 12
         tes0=dexp(-2.d0*zr)
         co2r=co2r*tes0
         co2i=co2i*tes0
   12    tes0=8.d0*dsqrt(zr*zr+zi*zi)
         z8r=zr*8.d0
         z8i=zi*8.d0
         nk=37
         psr=1.d0
         psi=0.
         pkr=1.d0
         pki=0.
         qsr=1.d0
         qsi=0.
         l=1
         do 10 k=1,nk
            l=-l
            qk=l
            kk=k-1
            ts=(2*kk+1)**2
            tt=k
            tes1=tt*tes0
            tes2=ts-fn2
            if (tes2.ge.tes1) go to 15
            ar=pkr*tes2
            ai=pki*tes2
            br=tt*z8r
            bi=tt*z8i
            tem1=br*br+bi*bi
            tem2=ai*br-ar*bi
            pkr=(ar*br+ai*bi)/tem1
            pki=tem2/tem1
            psr=psr+pkr
            psi=psi+pki
            qsr=qsr+qk*pkr
            qsi=qsi+qk*pki
   10    continue
   15    ar=cor*psr-coi*psi
         ai=cor*psi+coi*psr
         ansr=ar+co2r*qsr-co2i*qsi
         ansi=ai+co2r*qsi+co2i*qsr
         if (zi.eq.0.) ansi=0.
         anskr=pi*(qsr*c2r-qsi*c2i)
         anski=pi*(qsr*c2i+qsi*c2r)
!------------------------------------------                    
        
         end subroutine iln
!**********************************************
      subroutine zk24(zr,zi,n,zk1r,zk1i)
      implicit double precision (a-h,o-z)
!----------------------------------------------                    
      dimension x(12),c(12)
      data x( 1)/ .361913603606156014d+02/
      data x( 2)/ .276611087798460900d+02/
      data x( 3)/ .213967559361661098d+02/
      data x( 4)/ .164321950876753130d+02/
      data x( 5)/ .123904479638094713d+02/
      data x( 6)/ .907543423096120287d+01/
      data x( 7)/ .636997538803063507d+01/
      data x( 8)/ .419841564487841315d+01/
      data x( 9)/ .250984809723212788d+01/
      data x(10)/ .126958994010396155d+01/
      data x(11)/ .454506681563780283d+00/
      data x(12)/ .503618891172939534d-01/
      data c( 1)/ .332873699297821780d-15/
      data c( 2)/ .131692404861563402d-11/
      data c( 3)/ .609250853997510780d-09/
      data c( 4)/ .803794234988285940d-07/
      data c( 5)/ .431649140980466720d-05/
      data c( 6)/ .113773832728087596d-03/
      data c( 7)/ .164738496537683500d-02/
      data c( 8)/ .140967116201453420d-01/
      data c( 9)/ .748909410064614920d-01/
      data c(10)/ .255479243569118320d+00/
      data c(11)/ .572359070692886040d+00/
      data c(12)/ .853862327737398500d+00/
      if (n.ge.7) go to 77
      denr=1.d0
      deni=0.
      deni1=0.
      sz=dsin(-zi)
      cz=dcos(-zi)
      if (n.eq.0) go to 4
      do 3 i=1,n
         temp=denr*zr-deni*zi
         deni=denr*zi+deni*zr
         denr=temp
    3 continue
      denr1=denr/(denr*denr+deni*deni)
      if (zi.eq.0.) go to 4
      deni1=deni/(deni*deni+denr*denr)
    4 yr0=0.
      yi0=0.
      yr1=0.
      yi1=0.
      zrt=zr+zr
      zit=zi+zi
      do 10 i=1,12
         tr=zrt+x(i)
         temp=dsqrt(tr*tr+zit*zit)
         ar=dsqrt((tr+temp)*.5d0)
         if (ar) 1,2,1
    1    ai=zit/(2.d0*ar)
         go to 5
    2    ai=dsqrt(temp)
    5    if (n.ne.0) go to 6
         cc=ar*ar+ai*ai
         br0=ar/cc
         bi0=-ai/cc
         yr0=yr0+br0*c(i)
         yi0=yi0+bi0*c(i)
         go to 10
    6    br1=ar*x(i)
         bi1=ai*x(i)
         if (n.eq.1) go to 9
         br2=tr*x(i)
         bi2=zit*x(i)
         do 8 k=2,n
            temp=br1*br2-bi1*bi2
            tkm=2*k-1
            bi1=(br1*bi2+bi1*br2)/tkm
            br1=temp/tkm
    8    continue
    9    yr1=yr1+br1*c(i)
         yi1=yi1+bi1*c(i)
   10 continue
      et=dexp(-zr)
      if (n.ne.0) go to 12
      zk1r=(yr0*cz-yi0*sz)*et
      zk1i=(yi0*cz+yr0*sz)*et
      go to 20
   12 zk2r=yr1*denr1+yi1*deni1
      zk2i=yi1*denr1-yr1*deni1
      zk1r=(zk2r*cz-zk2i*sz)*et
      zk1i=(zk2r*sz+zk2i*cz)*et
   20 continue
      go to 88
   77 call zk32(zr,zi,n,zk1r,zk1i)
   88 continue
!------------------------------------------                    
     
      end  subroutine zk24
!******************************************
      subroutine zk32(zr,zi,n,zk1r,zk1i)
      implicit double precision (a-h,o-z)
!------------------------------------------                    
      dimension x(16),c(16)
      data x( 1)/ .507772238775370808d+02/
      data x( 2)/ .410816665254912019d+02/
      data x( 3)/ .337819704882261658d+02/
      data x( 4)/ .278314382113286757d+02/
      data x( 5)/ .228213006935252084d+02/
      data x( 6)/ .185377431786066933d+02/
      data x( 7)/ .148514313418012497d+02/
      data x( 8)/ .116770336739759564d+02/
      data x( 9)/ .895500133772339017d+01/
      data x(10)/ .664221517974144425d+01/
      data x(11)/ .470672670766758733d+01/
      data x(12)/ .312460105070214431d+01/
      data x(13)/ .187793150769607417d+01/
      data x(14)/ .953553155390865424d+00/
      data x(15)/ .342200156010947678d+00/
      data x(16)/ .379629145753134561d-01/
      data c( 1)/ .146213528547683240d-21/
      data c( 2)/ .184634730730365840d-17/
      data c( 3)/ .239468803418569740d-14/
      data c( 4)/ .843002042265289520d-12/
      data c( 5)/ .118665829267932772d-09/
      data c( 6)/ .819766432954179320d-08/
      data c( 7)/ .314833558509118800d-06/
      data c( 8)/ .730117025912475220d-05/
      data c( 9)/ .108331681236399652d-03/
      data c(10)/ .107253673105594410d-02/
      data c(11)/ .730978065330885620d-02/
      data c(12)/ .351068576631468600d-01/
      data c(13)/ .120916261911825228d+00/
      data c(14)/ .302539468153284960d+00/
      data c(15)/ .554916284605059800d+00/
      data c(16)/ .750476705185604780d+00/
      denr=1.d0
      deni=0.
      deni1=0.
      sz=dsin(-zi)
      cz=dcos(-zi)
      if (n.eq.0) go to 4
      do 3 i=1,n
         temp=denr*zr-deni*zi
         deni=denr*zi+deni*zr
         denr=temp
    3 continue
      denr1=denr/(denr*denr+deni*deni)
      if (zi.eq.0.) go to 4
      deni1=deni/(deni*deni+denr*denr)
    4 yr0=0.
      yi0=0.
      yr1=0.
      yi1=0.
      zrt=zr+zr
      zit=zi+zi
      do 10 i=1,16
         tr=zrt+x(i)
         temp=dsqrt(tr*tr+zit*zit)
         ar=dsqrt((tr+temp)*.5d0)
         if (ar) 1,2,1
    1    ai=zit/(2.d0*ar)
         go to 5
    2    ai=dsqrt(temp)
    5    if (n.ne.0) go to 6
         cc=ar*ar+ai*ai
         br0=ar/cc
         bi0=-ai/cc
         yr0=yr0+br0*c(i)
         yi0=yi0+bi0*c(i)
         go to 10
    6    br1=ar
         bi1=ai
         p=x(i)
         if (n.eq.1) go to 9
         br2=tr
         bi2=zit
         do 8 k=2,n
            p=p*x(i)
            temp=br1*br2-bi1*bi2
            tkm=2*k-1
            bi1=(br1*bi2+bi1*br2)/tkm
            br1=temp/tkm
    8    continue
    9    yr1=yr1+br1*(c(i)*p)
         yi1=yi1+bi1*(c(i)*p)
   10 continue
      et=dexp(-zr)
      if (n.ne.0) go to 12
      zk1r=(yr0*cz-yi0*sz)*et
      zk1i=(yi0*cz+yr0*sz)*et
      go to 20
   12 zk2r=yr1*denr1+yi1*deni1
      zk2i=yi1*denr1-yr1*deni1
      zk1r=(zk2r*cz-zk2i*sz)*et
      zk1i=(zk2r*sz+zk2i*cz)*et
   20 continue
  
!------------------------------------------                    
      
      end subroutine zk32
      
      end module bessel_mod

