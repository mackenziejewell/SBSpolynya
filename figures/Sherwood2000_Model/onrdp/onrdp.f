      program onrdp
c
c  One-dimensional (vertical) time-dependent 
c  model for currents, ice, and suspended sediment profiles.
c
c  This is the version used to generate newcases/1dcases for JGR manuscript.
c
c  Includes corrections to integration routines and ice thermo. from
c  Lars H. Smedsrud, and change in sediment flux b.c. calcs.
c
c  Chris Sherwood, USGS
c
      implicit none
      include 'onrdp.inc'              ! needed for array dimensions
      character*60 version/'Program onrdp, version of 23 March 1999'/

c     ...vertical water properties
      double complex cU(MAXZ)          ! double complex velocity arrays [m s-1]
      double complex cT(MAXZ+1)        ! internal shear stress [N m-2]
      double complex cbb(MAXZ), cdd(MAXZ), caa(MAXZ) ! working arrays
      double precision Km(MAXZ+1)      ! eddy viscosity for momentum [m2 s-1]
      double precision Kh(MAXZ+1)      ! eddy diffusivity for heat [m2 s-1]
      double precision Kr(MAXZ+1)      ! eddy diffusivity for mass [m2 s-1]
      double precision Rfn(MAXZ+1),Rin(MAXZ+1)!flux & gradient Richardson nums
      double precision Sm(MAXZ+1),Sh(MAXZ+1)! turb. closure params.
      double precision l(MAXZ+1)        ! turbulence length scale [m]
      double precision q(MAXZ+1)        ! 2*sqrt(TKE) [m/s] in MY model
      double precision tke(MAXZ+1)      ! tke = <u**^2> + <v'**2>; used in k-e model
      double precision diss(MAXZ+1)     ! dissipation rate in k-e model
      double precision Pr(MAXZ+1)       ! Prandtl no. (ratio Km/Kh)
      equivalence (q,tke),(Sm,diss),(Sh,Pr)  !q, Sm, and Sh are used with MY2.5 closure,
                                        ! tke, diss, and Pr are used with k-e model
      double precision Beta,sigt        ! ratio Kr/Km, neutral Prandtl no.  
      parameter(Beta=1.d0)
      double precision zetaff /1.d0/    ! shear fudge factor (Mellor, 1989)
      double precision Tw(MAXZ)         ! Temp [deg C]
      double precision S(MAXZ)          ! salinity [psu]
      double precision rho(MAXZ )       ! density [kg m-3]
      double precision z(MAXZ)          ! elevation of points [m ab]
      double precision zf(MAXZ+1)       ! elev. of interfaces [m ab]
      double precision xtran,ytran      ! cross-& along-shelf transports
      double precision xsflux(NSEDMAX),ysflux(NSEDMAX)
            ! cross- & along-shelf sediment transports
      double precision xiflux(NICEMAX),yiflux(NICEMAX)
            ! cross- & along-shelf ice transports
      double precision heat,salt,sedvol,icevol

c     ...time series
      double complex cTs(MAXT), cTb(MAXT)!surface &bottom shear stresses[N m-2]
      double precision t(MAXT)           ! time [s]
      double precision wu(MAXT), wv(MAXT)! east, north comp. wind speed [m s-1]
      double precision Hs(MAXT), Tm(MAXT)! wave height [m], period [s]
      double precision ub(MAXT)          ! near-bottom orbital velocity [m s-1]
      double precision Ta(MAXT)          ! air temperature (deg C)
      double precision snow(MAXT)        ! mass flux from atmos. [kg m-2 s-1]
      double precision Qw(MAXT)     ! surface heat flux [J m-2 s-1 = Watts m-2]
      double precision tausfmx(MAXT)     ! skin friction shear stress [N m-2]
      double precision Cd(MAXT)          ! bottom drag coef
      double precision dsdx(MAXT),dsdy(MAXT)
                                ! cross-& along-shelf water surf. slopes
c     ...ice and sediment
      double precision Ci(MAXZ,NICEMAX) ! concentration of ice 
      double precision Cs(MAXZ,NSEDMAX) ! concentration of sediment [0 to 1]
      double precision Cii(MAXZ)                ! working array for ice calcs
      double precision Tfri(MAXZ)               ! freezing point of ice [deg C]
      double precision ri(NICEMAX),D50(NSEDMAX)
                                ! radius of ice, diameter of sediment [m]
      double precision wi(NICEMAX)  ! rise velocity of ice (positive) [m s-1]
      double precision ws(NSEDMAX)
                          ! settling velocity for sediment (negative) [m s-1]
      double precision fri(NICEMAX)! frac. amount in each ice class [0-1]
      double precision frs(NSEDMAX)! frac. amount in each sediment class [0-1]
      double precision tauc(NSEDMAX)! critical shear stress [N m-2]
      integer nice /NICEMAX/! number of ice classes actually used
      integer nsed /NSEDMAX/! number of sed. classes actually used
      double precision Gt(MAXZ),Gc(MAXZ,NICEMAX),Gs(MAXZ)
                                ! source terms for T,ice, and S
      double precision Qi(NICEMAX) ! heat flux assoc. with ice formation
      double precision Av          ! ratio area to volume of ice xl
      double precision Nusi /6./   ! Nusselt number [], 
      double precision kwi ! ice xl heat exch. coeff in OS84 [J s-1 m-1 degC-1]
      double precision rhoi /920./ ! density of ice [kg m-3]
      double precision exl /0.02/  ! ratio disk thickness/disk diameter []
      double precision KT          ! ice heat xfer coefficient in JB95 [m2 s-1]
      double precision Li   ! latent heat of crystallization for ice [J kg-1]
      parameter(kwi=0.564,Li=3.35d5,KT=1.14d-6)
      double precision delT        ! water-ice temp. diff [deg C]
      double precision Ctot        ! temp variable for total ice density
      double precision rhos      ! sediment grain density [kg m-3]
      parameter(rhos=2650.)
      double precision Ca, xS    ! ref. conc. at zo [],excess shear stress []
      double precision Cb /.4/          ! bed concentration []
      double precision gamma0 /0.002/   ! resuspension coefficient []
      double precision c_we (3)         ! array of coefficients for erosion rate calc.
      double precision eflux(NSEDMAX),dflux(NSEDMAX),nflux(NSEDMAX)
                 ! erosion, deposition, and net sed fluxes at bed
      double precision fraci, fracs 
                 ! vol fractions of ice & sediment for density calcs
      double precision wspeed, wdir
                 ! wind speed [m s-1], wind & wave direction [deg T]
      double precision Kconst /0.08/
                 ! value for constant eddy viscosity [m2 s-1]
      double precision nu /2.d-6/ !molecular kinematic viscosity (min Km)
      double precision vk,g
                 ! von Karmans const, grav. accel.
      parameter(vk=0.41, g=9.80665)
      double precision betag           ! g/rhobar
      double precision alpha /0.2/     ! parameter in master length scale calc.
      double precision xamp,yamp,xlag,ylag ! ss slope parameters
      double complex czero /(0., 0.)/              
      double precision NaN /1.d35/                    
      double precision Cdmin /0.0025/ ! value for min bottom drag coeff [ ]
      double precision ustrs, ustrb
                 ! friction velocity at surface and bottom [m s-1]
      double precision zos, zob    ! roughness length at surface and bottom [m]
      double precision zoa, zobmin ! apparent rough. leng, min. zo [m]
      double precision h, dz       ! water depth [m] and grid spacing [m]
      double precision ks /0.0001/       ! nominal bottom grain size [m] 
      double precision f /1.d-4/         ! coriolis frequency [s-1]
      double precision rlat /60./        ! latitude [deg N]
      double precision fetch /50000./    ! effective wind fetch distance [m]
      double precision rhow, rhobar /1025./ ! density water, mean dens [kg m-3]
      double precision CH, rhoa
                  ! air/water heat exchange coef, density air [kg m-3]
      double precision Hfcoef        ! heat flux coefficient [?]
      double precision cpa, cpw, cpi
                  ! heat capacity of air, water [J kg-1 degC-1]
      parameter(rhoa=1.225,cpa=1004.6,cpw=3986.5,cpi=2086.7)
      double precision rho_off        ! offset to improve calcs in drho/dz
      parameter(rho_off=1000.)
      double precision t0,dt,dtdz ! start time [s], time step [s], dt/(dz*dz)
      integer nt,nz        ! number of time steps, num. of grid cells
      integer i,j,n                ! loop indices
      integer ifail /0/            ! error flag
      integer iverbose,ibbc,iktype ! flags for output, bot. bc, and turb. clos.
      integer itfreeze /0/     ! freeze point (0=tfr(S,P), 1=tfr(S,0), 2=const)
      integer iceflg /0/       ! 0=OS84 1=JB95 2=Thermo 3=Sherwood
      integer ixshelf /1/      ! flag for pressure gradient (0=input,1=xshelf balance, 2=modlulate with xamp, yamp)
      integer iss /0/          ! sus. sed. bbc (0=ref. conc., 1=flux bc)
      integer idensflg         ! density calculation (0=full eqn.,1=quadr.appx)
      integer isrestart /0/    ! flag that ques readinit to read extra columns
      integer nzinc            ! number of output levels
      integer ntinc /3600/     ! frequency of output
      integer ibigLflag /0/    ! 0=parabolic integral length scale (standard)
                               ! 1=bilinear 
      integer ilenflag /0/     ! 0=fixed parabolic or bilinear length scale
                               ! 1=q2l/q2 in MY25 (standard)
      integer iiwmflag /0/     ! 0=standard MY25
                               ! 1=add internal-wave mixing fudge
      integer lastz
      character*40 paramfile /'onrparam.dat'/
      character*40 initfile /'onrinit.dat'/
      character*40 forcefile /'onrforce.dat'/

c     ...function declarations
      double precision coriolis, density2, tfreeze, we_calc

      write(*,1010)version
 1010 format(1x,a)

c     ...read parameter file
      open(60,file=paramfile,status='old')
      call ireadfil3(60,'nzinc',5,11,1,nzinc)
      call ireadfil3(60,'ntinc',5,60,1,ntinc)
      call ireadfil3(60,'nt',2,3601,1,nt)
      call ireadfil3(60,'ixshelf',7,1,1,ixshelf)
      call ireadfil3(60,'iverbose',8,1,1,iverbose)
      call ireadfil3(60,'isrestart',9,0,1,isrestart)
      call ireadfil3(60,'iceflg',6,0,1,iceflg)
      call ireadfil3(60,'ibbc',4,1,2,ibbc)
      call ireadfil3(60,'ibigLflag',9,0,1,ibigLflag)
      call ireadfil3(60,'iiwmflag',8,0,1,iiwmflag)
      call dreadfil3(60,'zetaff',6,1.d0,1,zetaff)
      call dreadfil3(60,'sigt',4,1.2d0,1,sigt)
      call ireadfil3(60,'ilenflag',8,1,1,ilenflag)
      call ireadfil3(60,'iss',3,0,1,iss)
      call ireadfil3(60,'itfreeze',8,2,1,itfreeze)
      call ireadfil3(60,'idensflg',8,1,1,idensflg)
      call ireadfil3(60,'iktype',6,4,1,iktype)

      call ireadfil3(60,'nice',4,1,0,nice)
      if(nice.gt.NICEMAX)then
         write(*,*) 'nice must be <= ',NICEMAX
         stop
      endif
      do j=1,nice
         read(60,*) fri(j),ri(j),wi(j)
      enddo

      call ireadfil3(60,'nsed',4,2,0,nsed)
      if(nice.gt.NSEDMAX)then
         write(*,*) 'nsed must be <= ',NSEDMAX
         stop
      endif
      do j=1,nsed
         read(60,*) frs(j),D50(j),ws(j),tauc(j)
      enddo

      call dreadfil3(60,'nu',2,2.d-6,1,nu)
      call dreadfil3(60,'Cb',2,0.4d0,1,Cb)
      call dreadfil3(60,'gamma0',6,0.002d0,1,gamma0)
      call dreadfil3(60,'c_we(1)',7,.23d-7,1,c_we(1))
      call dreadfil3(60,'c_we(2)',7,1.d0,1,c_we(2))
      call dreadfil3(60,'c_we(3)',7,0.d0,1,c_we(3))
      call dreadfil3(60,'h',1,20.d0,1,h)
      call dreadfil3(60,'ks',2,0.0015d0,1,ks)
      call dreadfil3(60,'rlat',4,75.d0,1,rlat)
      call dreadfil3(60,'fetch',5,1.d5,1,fetch)
      call dreadfil3(60,'xamp',4,0.d0,2,xamp)
      call dreadfil3(60,'xlag',4,0.d0,2,xlag)
      call dreadfil3(60,'yamp',4,0.d0,2,yamp)
      call dreadfil3(60,'ylag',4,0.d0,2,ylag)
      call dreadfil3(60,'dz',2,0.2d0,1,dz)
      call dreadfil3(60,'dt',2,5.d0,1,dt)
      call dreadfil3(60,'Cdmin',5,0.0025d0,1,Cdmin)
      call dreadfil3(60,'rhobar',6,1025.d0,1,rhobar)
      call dreadfil3(60,'zos',3,0.05d0,1,zos)
      call dreadfil3(60,'zobmin',6,0.00005d0,1,zobmin)
      call dreadfil3(60,'Nusi',4,6.d0,1,Nusi)
      call dreadfil3(60,'Av',2,3.d0,1,Av)
      call dreadfil3(60,'exl',3,1.d0,1,exl)
      call dreadfil3(60,'rhoi',3,920.d0,2,rhoi)
      call dreadfil3(60,'alpha',5,0.2d0,1,alpha)
      if(iktype.eq.0) call dreadfil3(60,'Kconst',6,0.08d0,1,Kconst)
      close(60)

      betag = g/rhobar

      write(*,1020)h,zobmin,zos
 1020 format(1x,'Depth: ',f5.1,'m  zobmin: ',f9.6,'m  zos: ',f6.3,'m')
      write(*,1030)Cdmin
 1030 format(1x,'Cdmin: ',f7.4)
      if(Cdmin .ge. 0.5*dz ) stop 'Cdmin >= dz/2'
      f = coriolis(rlat)
      write(*,1040) rlat,f
 1040 format(1x,'Latitude: ',f5.1,'deg.  Coriolis f: ',f7.5,'s-1.')
      write(*,*) ' Sediment properties: '
      do j=1,nsed
         write(*,1045) frs(j),D50(j)*1000.,ws(j)*100.,tauc(j)
 1045    format(2x,'fr :',f5.3,' D50: ',f6.4,' mm  ws: ',
     &        f7.4,' cms-1  tauc: ',f5.2,' N m-2')
      enddo
      write(*,1049) Cb
 1049 format(1x,'Cb: ',f5.2)
      if(iss.eq.0) write(*,1052) gamma0
 1052 format(1x,'Reference boundary condition. gamma0: ',f6.4)
      if(iss.eq.1) write(*,1053) c_we(1),c_we(2),c_we(3)
 1053 format(1x,'Flux boundary condition. Coefficients c_we:',
     &     3(2x,f12.5))

      t0 = 0.
      if(nt.gt.MAXT)then
         nt = MAXT
         write(*,*) 'Warning: nt exceeds array size...nt decreased.'
      endif
      do n=2,nt
         t(n)=t0+(n-1)*dt
      enddo
      write(*,1050)nt,dt
 1050 format(1x,'nt: ',i7,'      dt:',f6.2,' seconds')
      write(*,1060) t0,t(nt)
 1060 format(1x,'t0: ',f8.0,', t(nt): ',f8.0,' seconds')

c     ...read and interpolate surface forcing file
      call readforce( forcefile,nt,t,wu,wv,Ta,snow,dsdx,dsdy,
     &    ixshelf)
      write(*,*) 'Input time series interpolated.'
      call wstress(nt, wu, wv, cTs )
      write(*,*) 'Windstress calculated.'
      call waves(nt, h, fetch, wu, wv, Hs, Tm, ub)
      write(*,*) 'Waves calculated.'
      if(ixshelf.eq.1)then
         call makeslope(nt,t,xamp,yamp,xlag,ylag,dsdx,dsdy)      
         write(*,*) 'Surface slope calculated.'
      endif

      if(iverbose.eq.1)then
         write(*,*) ' Interpolated forcing time series:'
         write(*,*)
     &' Time,Wind U,V,Wind Twx,Twy,|Tw|,AirT,Hs,Tm,ub,snow,dsdx,dsdy'
         do i=1,nt,nt/10
            write(*,1110)t(i),wu(i),wv(i),dreal(cTs(i)),dimag(cTs(i)),
     &           abs(cTs(i)),Ta(i),Hs(i),Tm(i),ub(i),snow(i),
     &           dsdx(i),dsdy(i)
 1110       format(1x,f9.1,9(1x,f7.2),1x,e9.2,2(1x,f9.5))
         enddo
         i=nt
         write(*,1110) t(i),wu(i),wv(i),dreal(cTs(i)),dimag(cTs(i)),
     &        abs(cTs(i)),Ta(i),Hs(i),Tm(i),ub(i),snow(i),
     &        dsdx(i),dsdy(i)
      endif
      
c     ...depth grid
      nz = nint(h/dz)
      dz = h/nz
      if(nz.gt.MAXZ) then
         dz = h /(MAXZ)
         nz = MAXZ
         write(*,*) 'Warning: dz too small for arrays...dz increased.'
      endif
      write(*,2010)nz,dz
 2010 format(1x,'nz:  ',i4,'      dz: ',f7.2,' meters')

      dtdz = dt/(dz*dz)
c     ...properties at points z
      z(1) = 0.5*dz
      do i=2,nz
         z(i)=z(i-1)+dz
      enddo
c     ...vertical fluxes at faces zf
      zf(1) = 0.
      do i=2,nz+1
         zf(i)=zf(i-1)+dz
      enddo
      write(*,2020)z(1),z(nz)
 2020 format(1x,' z(1): ',f6.2,',    z(nz): ',f6.2)
      write(*,2030)zf(1),zf(nz+1)
 2030 format(1x,'zf(1): ',f6.2,', zf(nz+1): ',f6.2)
      if(nzinc.gt.nz)nzinc=nz
      write(*,2033) nzinc
 2033 format(1x,i4,' output levels for profiles.')

c     ...check CFL criteria for seds and ice
      do i=1,nsed
         if( abs(ws(i))*dt .gt. dz )then
            write(*,*) 'ERROR: Sed. fall velocity exceeds CFL'
            write(*,*) 'ws(',i,') = ',ws(i)
            stop 'Punt.'
         endif
      enddo
      do i=1,nice
         if( abs(wi(i))*dt .gt. dz )then
            write(*,*) 'ERROR: Ice rise velocity exceeds CFL'
            write(*,*) 'wi(',i,') = ',wi(i)
            stop 'Punt.'
         endif
      enddo

c     ...read, interpolate, and initialize vertical profiles
      call readinit(initfile, isrestart, rho_off, nice, nsed,
     &     rhoi, rhos, sigt, zob, zos, vk, kconst, beta,
     &     nz, h, z, zf, Tw, S, cU, Ci, Cs, rho, Km, Kr, Kh,
     &     Tfri, cT, q, l, Rin, Rfn, Sm, Sh )

c     ...initialize for first time step
      n=1
      t(n)=t0
      cTb(1) = czero

      zob=zobmin
      if((ks/30.) .lt. zobmin) ks = 30.*zobmin
      if(zob.ge.z(1)) stop 'zob >= z(1)'

c     ...integrate velocity to determine transport
      call integ( nz, dz, zos, zob, cU, xtran, ytran )
      if(ixshelf.eq.1)
     &     dsdx(n) = dsdx(n)+(xtran/h)/(g*dt)

c     ... bottom stress
      if(ibbc.eq.0)then
         cTb(n) = czero
      elseif(ibbc.eq.1)then
         if(zob.le.0.) Cd(n) = 0.
         if(zob.gt.0.) Cd(n) = (vk/(log(z(1)/zob)))**2 !Actually, 2*Cd
         Cd(n) = max( Cd(n), Cdmin)
         cTb(n) = rhobar*Cd(n)*cU(1)*abs(cU(1)) ! 1/2 cancels here
      elseif(ibbc.eq.2)then
         call pcoord(wu(n),wv(n),wspeed,wdir)
         call bstress(Hs(n),Tm(n),wdir,dreal(cU(1)),dimag(cU(1)),
     &        zf(2),ks,h,rhobar,zoa,Cd(n),tausfmx(n))
         Cd(n) = max( Cd(n), Cdmin)
         cTb(n) = rhobar*Cd(n)*cU(1)*abs(cU(1)) ! 1/2 cancels here
         zob=max(zob,zoa)
         if( zob.gt.z(1) )then
            write(*,*) 'Warning: zob approaching z(1).'
            zob = 0.5*z(1)
         endif
      endif

c     ...heat flux
      wspeed = sqrt(wu(n)**2+wv(n)**2)
      CH = 1.10d-3
      if( Tw(nz).le.Ta(n) ) CH = 0.83d-3
      Hfcoef = rhoa*cpa*CH*wspeed
      Qw(n) = Hfcoef*( Ta(n)-Tw(nz) )
      
c     ...output files
      open(70,file='onrprof.dat',status='unknown')
      open(71,file='onrts.dat',status='unknown')

c        ...write output profiles
      do i=nz,1,-nz/nzinc
         write(70,1301) z(i),Tw(i),S(i),Km(i),
     &        dreal(cU(i)),dimag(cU(i)),
     &        (max(Ci(i,j)-1.d-10,0.d0),j=1,nice),
     &        (max(Cs(i,j)-1.d-10,0.d0),j=1,nsed),
     &        rho(i),Tfri(i),abs(cT(i)),
     &        q(i),l(i),Rin(i),Sm(i),Sh(i)
         lastz=i
      enddo
      if(lastz.ne.1)then
         i=1
         write(70,1301) z(i),Tw(i),S(i),Km(i),
     &        dreal(cU(i)),dimag(cU(i)),
     &        (max(Ci(i,j)-1.d-10,0.d0),j=1,nice),
     &        (max(Cs(i,j)-1.d-10,0.d0),j=1,nsed),
     &        rho(i),Tfri(i),abs(cT(i)),
     &        q(i),l(i),Rin(i),Sm(i),Sh(i)
      endif
c     ...integrate velocity to determine transport
      call integ( nz, dz, zos, zob, cU, xtran, ytran )
      if(ixshelf.eq.1) dsdx(n) = dsdx(n)+(xtran/h)/(g*dt)
c     ...integrate other quantities
      call conserve(nz,dz,nice,nsed,S,Tw,Ci,Cs,Li,cpw,rhobar,
     &     heat,salt,sedvol,icevol)
      do j=1,nsed
         call sflux(nz,dz,zos,zob,cU,Cs,j,xsflux(j),ysflux(j))
      enddo
      do j=1,nice
         call iflux(nz,dz,zos,zob,cU,Ci,j,xiflux(j),yiflux(j))
      enddo
      write(71,7100) t(n),abs(cTs(n)),Hs(n),
     &     Tm(n),ub(n),Cd(n),cTb(n),tausfmx(n),Qw(n),
     &     dsdx(n),dsdy(n),
     &     (xiflux(j),yiflux(j),j=1,nice),
     &           (xsflux(j),ysflux(j),j=1,nsed),
     &     heat,salt,icevol,sedvol


c     ...main loop over time 
      do n=2,nt
         
c        ...calculate stress profile (diagnostic only)
         cT(nz+1) = cTs(n)
         do i=nz,2,-1
            cT(i)= rhobar*Km(i)*(cU(i)-cU(i-1))/dz
         enddo
         cT(1)=cTb(n-1)
         
c        ...write output profiles
         if(mod(nint(t(n)),ntinc).eq.0)then
            do i=nz,1,-nz/nzinc
               write(70,1301) z(i),Tw(i),S(i),Km(i),
     &              dreal(cU(i)),dimag(cU(i)),
     &              (max(Ci(i,j)-1.d-10,0.d0),j=1,nice),
     &              (max(Cs(i,j)-1.d-10,0.d0),j=1,nsed),
     &              rho(i),Tfri(i),abs(cT(i)),
     &              q(i),l(i),Rin(i),Sm(i),Sh(i)
 1301          format(1x,f6.2,1x,g10.4,1x,f6.3,1x,f8.4,
     &              2(1x,f7.3),3(1x,e9.2),
     &              1x,f8.3,1x,f5.2,1x,g9.3,3(1x,g10.4),2(1x,g10.4))
            lastz=i
            enddo
            if(lastz.ne.1)then
               i=1
               write(70,1301) z(i),Tw(i),S(i),Km(i),
     &              dreal(cU(i)),dimag(cU(i)),
     &              (max(Ci(i,j)-1.d-10,0.d0),j=1,nice),
     &              (max(Cs(i,j)-1.d-10,0.d0),j=1,nsed),
     &              rho(i),Tfri(i),abs(cT(i)),
     &              q(i),l(i),Rin(i),Sm(i),Sh(i)
            endif
         endif

c        ...Ice (see Omstedt and Svensson)
c        ...snow into top layer
         do j=1,nice
            Ci(nz,j)=Ci(nz,j)+fri(j)*snow(n)*dt/(rhoi*dz)
         enddo

c        ...freezing point of water
         if(itfreeze.eq.0)then
            do i=1,nz
               Tfri(i)=tfreeze( S(i), 0.1*(h-z(i)) )
            enddo
         elseif(itfreeze.eq.1)then
            Tfri(nz)=tfreeze( S(nz),0.1*( h-z(nz) ) )
            do i=nz-1,1,-1
               Tfri(i)=Tfri(nz)
            enddo
         elseif(itfreeze.eq.2)then
c           ...no changes...freezing point is constant
         else
            stop 'bad value of itfreeze'
         endif

c        ...ice calculations
         if(iceflg.eq.0)then
c           ...according to Omstedt and Svensson
            do i=1,nz
               Gt(i)= 0.0
               Gs(i)= 0.0
               do j=1,nice
                  delT=max(0.d0, Tfri(i)-Tw(i)) !freezing
                  if(Tw(i) .gt. 0.) delT= -Tw(i) !melting
                  Qi(j) = delT*Ci(i,j)*(Av/ri(j))*Nusi*kwi
     &                 /(2.*ri(j))
                  Gt(i)=Gt(i)+Qi(j)/(rhobar*cpw)
                  Gc(i,j)=Qi(j)/(rhoi*Li)
                  Gs(i)=Gs(i)+S(i)*Gc(i,j)
               enddo
            enddo
         elseif(iceflg.eq.1)then
c           ...according to Jenkins and Bombosch
            do i=1,nz
               Gt(i)= 0.0
               Gs(i)= 0.0
               do j=1,nice
                  delT=max(0.d0,Tfri(i)-Tw(i)) !freezing
                  if(Tw(i) .gt. 0.) delT= -Tw(i) !melting
                  Qi(j) = delT*Ci(i,j)*(Av/ri(j))*Nusi*KT*cpw
     &                 /(exl*ri(j))
                  Gt(i)=Gt(i)+Qi(j)/(rhobar*cpw)
                  Gc(i,j)=Qi(j)/(rhobar*Li)
                  Gs(i)=Gs(i)+S(i)*Gc(i,j)
               enddo
            enddo
         elseif(iceflg.eq.2)then
c           ...instant thermal eq
c           ...only works for one size class
c           ...work on volumes (factor out dt for a while)
            do i=1,nz
               if( Tw(i) .gt. 0.d0 )then
c                 ...melting
                  Gc(i,1) = -Tw(i)*cpw/Li
                  if( Ci(i,1)+Gc(i,1) .lt. 0. )then
                     Gc(i,1) = -Ci(i,1)
                     Gt(i) = Tw(i) - Ci(i,1)*rhobar*dz*Li/cpw
                  else
                     Gc(i,1) = -Tw(i)*cpw/Li
                     Gt(i) = -Tw(i)
                  endif
               elseif( Tw(i) .lt. Tfri(i) )then
c                 ...freezing
                  Gc(i,1) = (Tfri(i)-Tw(i))*cpw/Li
                  Gt(i) = (Tfri(i)-Tw(i))
               else
c                 ...no change
                  Gc(i,1) = 0.
                  Gt(i) = 0.
               endif
c              ...convert volumes to rates
               Gc(i,1) = Gc(i,1)/dt
               Gt(i) = Gt(i)/dt
               Gs(i)=S(i)*Gc(i,1)
            enddo
         elseif(iceflg.eq.3)then
c           ...according to Sherwood
            do i=1,nz
               Gt(i)= 0.0
               Gs(i)= 0.0
               Ctot = 0.0
               do j=1,nice
                  Ctot = Ctot+Ci(i,j)
               enddo
               do j=1,nice
                  delT=max(0.d0, Tfri(i)-Tw(i)) !freezing
                  if(Tw(i) .gt. 0.d0) delT= -Tw(i) !melting
                  Qi(j) = delT*Ci(i,j)*(Av/ri(j))*(rhobar/rhoi)
     &                 *Nusi*kwi/(exl*2.*ri(j))
                  if( delT .gt. 0. )then
                     Gc(i,j)=min(Qi(j)/(Li*rhobar*(1.-Ctot)),
     &                    delT*cpw/(Li*dt))
                     Gt(i) = Gt(i)+
     &                    min(Qi(j)/(cpw*rhobar*(1.-Ctot)),delT/dt)
                  else
                     Gc(i,j)=max(Qi(j)/(Li*rhobar*(1.-Ctot)),
     &                    delT*cpw/(Li*dt))
                     Gt(i) = Gt(i)+
     &                    max(Qi(j)/(cpw*rhobar*(1.-Ctot)),delT/dt)
                  endif
                  Gs(i)=Gs(i)+S(i)*Gc(i,j)/(1.-(Ctot+Gc(i,j)))
               enddo
            enddo
         else
            stop 'Bad value of iceflg.'
         endif
         
         do j=1,nice
c           ...copy ice conc. to 1d array w/ fluxes
            do i=1,nz
               Cii(i)=Ci(i,j)+dt*Gc(i,j)
            enddo
            call advect( nz, dt, dz, wi(j), Cii )
            call diffus(1,nz,z,zf,Kr,Cii,dt,0.d0,0.d0)
 
c           ...copy answer back to ice array
            do i=1,nz
               Ci(i,j)=Cii(i)
            enddo
         enddo    

c        ...Sediment
         do j=1,nsed
            xS = max(0.d0 , (tausfmx(n-1)-tauc(j))/tauc(j))
            if(iss.eq.0)then
c              ...ref. concentration
               Ca = frs(j)*Cb*gamma0*xS/(1+gamma0*xS)
               Cs(1,j)=0.d0
               if(xS.gt.0.d0 .and. ustrb .gt. 0.d0)
     &            Cs(1,j)=Ca*(z(1)/zob)**(ws(j)/(vk*ustrb))
            elseif(iss.eq.1)then
c              ...sediment flux
c*********************************************************
c               eflux(j)=dt*we*frs(j)*Cb*xS**eta
c               dflux(j)=dt*ws(j)*Cs(1,j)
c               Cs(1,j)=Cs(1,j)+eflux(j)+dflux(j)
c*********************************************************
               dflux(j)=ws(j)*Cs(1,j)
               eflux(j)=we_calc(frs(j),tauc(j),tausfmx(n-1),c_we)
c               if((eflux(j)+dflux(j)).gt.0.)then
c                 ...net erosion limited by material in bed
c                  nflux(j)=min(zbs(1,j)/dt,eflux(j)+dflux(j))
c               else
c                 ...net deposition limited by material in wc
                  nflux(j)=max(-Cs(1,j)*dz/dt,eflux(j)+dflux(j))
c               endif
c              ...material to wc
               Cs(1,j)= (Cs(1,j)*dz+(dt*nflux(j)))/dz
            else
               stop 'bad value of iss'
            endif
c           ...copy sed conc. to 1d array w/ bb flux
            do i=1,nz
               Cii(i)=Cs(i,j)
            enddo
            call advect( nz, dt, dz, ws(j), Cii)
            call diffus(1,nz,z,zf,Kr,Cii,dt,0.d0,0.d0)
c           ...copy answer back to sed array
            do i=1,nz
               Cs(i,j)=Cii(i)
            enddo
         enddo

 
c        ...Heat flux (this could be made implicit)
         wspeed = sqrt(wu(n)**2+wv(n)**2)
         CH = 1.10d-3
         if( Tw(nz).le.Ta(n))CH = 0.83d-3
         Hfcoef = rhoa*cpa*CH*wspeed
         Qw(n) = Hfcoef*( Ta(n)-Tw(nz) )

c        ...Temperature
         do i=1,nz-1             
            Tw(i) = Tw(i)+dt*Gt(i)
         enddo
         Tw(nz) = Tw(nz)+dt*Gt(nz)+dt*Qw(n)/(dz*rhobar*cpw)
         call diffus(1,nz,z,zf,Kh,Tw,dt,0.d0,0.0d0)

c        ...Salinity
         do i=1,nz  
            S(i) = S(i)+dt*Gs(i)
         enddo
         call diffus(1,nz,z,zf,Kr,S,dt,0.d0,0.d0)

c        ...Density
         do i=1,nz
c           ...quadratic eqn. for density
            rhow = density2(S(i),Tw(i))
            fraci = 0.
            do j=1,nice
               fraci = fraci + Ci(i,j)
            enddo
            fracs = 0.
            do j=1,nsed
               fracs = fracs + Cs(i,j)
            enddo
            rho(i)=fraci*(rhoi-rho_off)
     &           +fracs*(rhos-rho_off)
     &           +(1.-(fraci+fracs))*(rhow-rho_off)
c            if(i.eq.1.or.i.eq.nz)
c     &           write(*,1140) fraci, fracs,
c     &           rhow, rho(i)
 1140       format(1x,'Ice: ',f8.5,' Sed: ',f8.5,
     &           ' rhow: ',f5.2,' rho(i): ',f5.2)
         enddo

c        ...Velocity
         do i=2,nz-1             
            cbb(i) = cmplx( -dtdz*Km(i), 0. )
            cdd(i) = cmplx( dtdz*(Km(i+1)+Km(i))+1,0.5*f*dt )
            caa(i) = cmplx( -dtdz*Km(i+1), 0. )
            cU(i) =  cU(i)*cmplx(1.,-0.5*f*dt )
     &           -dt*g*cmplx(dsdx(n-1),dsdy(n-1))
         enddo
         
c        ..bbc
         if(ibbc.eq.0)then
c           ...no-slip
            cbb(1) = NaN
            cdd(1) = 1.
            caa(1) = 0.
            cU(1) = 0.
         elseif((ibbc.eq.1).or.(ibbc .eq.2))then
c           ...bottom stress specified
            cbb(1) = NaN
            cdd(1) = cmplx(  dtdz*Km(2)+1, 0.5*f*dt)
            caa(1) = cmplx( -dtdz*Km(2), 0. )
            cU(1)  = cU(1)*cmplx( 1, -0.5*f*dt)-dt*(cTb(n-1)/rhobar)/dz
     &           -dt*g*cmplx(dsdx(n-1),dsdy(n-1))
         else
            stop 'bad value of ibbc'
         endif

c     ...surf bc
         cbb(nz) = cmplx( -dtdz*Km(nz), 0. )
         cdd(nz) = cmplx(  dtdz*Km(nz)+1, 0.5*f*dt )
         caa(nz) = NaN
         cU(nz)  = cU(nz)*cmplx(1,-0.5*f*dt)+dt*(cTs(n)/rhobar)/dz
     &        -dt*g*cmplx(dsdx(n-1),dsdy(n-1))
         if(.FALSE.)then
            do i=nz,1,-1
               write(*,1100) z(i),Km(i),cbb(i),cdd(i),caa(i),cU(i)
 1100          format(1x,f5.1,1x,f8.3,1x,4(1x,2(1x,f9.4)))
            enddo
         endif
         
         call csy(1,nz,cbb,cdd,caa,cU)
         
c        ...update bottom stress
         if(ibbc.eq.0)then
            cTb(n) = czero
         elseif(ibbc.eq.1)then
            if(zob.le.0.) Cd(n) = 0.
            if(zob.gt.0.) Cd(n) = (vk/(log(z(1)/zob)))**2 !Actually, 2*Cd
            Cd(n) = max( Cd(n), Cdmin)
            cTb(n) = rhobar*Cd(n)*cU(1)*abs(cU(1)) ! 1/2 cancels here
c           ...calculate wave/current bottom stress and zoa
         elseif(ibbc.eq.2)then
            call pcoord(wu(n),wv(n),wspeed,wdir)
            call bstress(Hs(n),Tm(n),wdir,dreal(cU(1)),dimag(cU(1)),
     &           zf(2),ks,h,rhobar,zoa,Cd(n),tausfmx(n))
            Cd(n) = max( Cd(n), Cdmin)
            cTb(n) = rhobar*Cd(n)*cU(1)*abs(cU(1)) ! 1/2 cancels here
            zob=max(zob,zoa)
            if( zob.ge.z(1) )then
               write(*,*) 'Warning: zob >=  z(1).'
               zob = 0.5*z(1)
            endif
         else
            stop 'bad value of ibbc'
         endif
            
c        write(*,1200) t(n),cU(1),abs(cTb(n)),Cd(n)
 1200    format(1x,'t:',f8.0,' B.vel.:',2(1x,f7.4),
     &        ' Tb:',1x,f7.4,' Cd:',f7.4)

c     ...new eddy viscosity
         if(iktype.eq.0)then
c           ...constant
            do i=1,nz+1
               Km(i)=Kconst
               Kh(i)=Km(i)
               Kr(i)=Beta*Km(i)
            enddo
         elseif(iktype.eq.1)then
c           ...bilinear (Madsen, 1977)
            ustrb = sqrt(abs(cTb(n))/rhobar)
            ustrs = sqrt(abs(cTs(n))/rhobar)
            do i=1,nz+1
               if( (zf(i)/h) .le. ustrb/(ustrb+ustrs) )
     &              Km(i)=ustrb*vk*zf(i)+nu
               if( ((h-zf(i))/h) .le. ustrs/(ustrb+ustrs) )
     &              Km(i)=ustrs*vk*(h-zf(i))+nu
               Kh(i)=Km(i)
               Kr(i)=Beta*Km(i)
            enddo
         elseif(iktype.eq.2)then
c           ...cubic (Signell et al, 1990)
            ustrb = sqrt(abs(cTb(n))/rhobar)
            ustrs = sqrt(abs(cTs(n))/rhobar)
            do i=1,nz+1
               Km(i)= -vk*ustrs*(zf(i)-h)
     &              + (vk*(ustrb-2.*ustrs)/h)*(zf(i)-h)**2
     &              + (vk*(ustrb-ustrs)/(h*h))*(zf(i)-h)**3
     &              + nu
               Kh(i)=Km(i)
               Kr(i)=Beta*Km(i)
            enddo
         elseif(iktype.eq.3)then
c           ...k-e closure
            ustrb = sqrt(abs(cTb(n))/rhobar)
            ustrs = sqrt(abs(cTs(n))/rhobar)
            call kepsilon(nz,dt,dz,zf,tke,diss,l,Km,Kh,Kr,cU,
     &           rho,Pr,betag,Beta,sigt,vk,nu,h,
     &           ustrs,ustrb,zos,zob,
     &           zetaff,ilenflag,iiwmflag)
         elseif(iktype.eq.4)then
c           ...MY2.5 closure
            ustrb = sqrt(abs(cTb(n))/rhobar)
            ustrs = sqrt(abs(cTs(n))/rhobar)
            call MY25(nz,dt,dz,zf,q,l,Km,Kh,Kr,Rin,Sm,Sh,cU,
     &      rho,betag,Beta,vk,nu,h,ustrs,ustrb,zos,zob,
     &           zetaff,ilenflag,ibigLflag,iiwmflag)
         else
            stop 'bad value of iktype'
         endif

c        ...integrate velocity to determine transport
         call integ( nz, dz, zos, zob, cU, xtran, ytran )
         if(ixshelf.eq.1) dsdx(n) = dsdx(n)+(xtran/h)/(g*dt)
c        ...integrate other quantities
         call conserve(nz,dz,nice,nsed,S,Tw,Ci,Cs,Li,cpw,rhobar,
     &     heat,salt,sedvol,icevol)

         if( mod(nint(t(n)),600).eq.0 )then
            do j=1,nsed
               call sflux(nz,dz,zos,zob,cU,Cs,j,xsflux(j),ysflux(j))
            enddo
            do j=1,nice
               call iflux(nz,dz,zos,zob,cU,Ci,j,xiflux(j),yiflux(j))
            enddo
            write(71,7100) t(n),abs(cTs(n)),Hs(n),
     &           Tm(n),ub(n),Cd(n),cTb(n),tausfmx(n),Qw(n),
     &           dsdx(n),dsdy(n),
     &           (xiflux(j),yiflux(j),j=1,nice),
     &           (xsflux(j),ysflux(j),j=1,nsed),
     &           heat,salt,icevol,sedvol
 7100       format(1x,f8.0,f7.4,2(1x,f4.1),1x,f4.2,
     &           1x,f7.5,3(1x,f8.5),1x,f12.4,
     &           2(1x,e12.4),6(1x,e12.4),
     &           4(1x,e14.4))
         endif
      enddo

      close(70)
      close(71)
     
      end
c***********************************************************************
      subroutine MY25(nz,dt,dz,zf,q,l,Km,Kh,Kr,GH,Sm,Sh,cU,
     &     rho,betag,Beta,vk,nu,h,ustrs,ustrb,zos,zob,
     &     zetaff,ilenflag,ibigLflag,iiwmflag)
c***********************************************************************
      implicit none
      include 'onrdp.inc'
      integer  k,nz
      double precision rho(MAXZ),zf(MAXZ+1)
      double complex   cU(MAXZ)
      double precision ustrs,ustrb,vk,nu,betag,Beta
      double precision dt,dz,h,zos,zob
      double precision q(MAXZ+1),l(MAXZ+1)
      double precision Ps(MAXZ+1),Pb(MAXZ+1)
      double precision bigW(MAXZ+1),bigL(MAXZ+1)
      double precision q2(MAXZ+1),q2l(MAXZ+1),N2(MAXZ+1)
      double precision Km(MAXZ+1),Kr(MAXZ+1),Kh(MAXZ+1)
      double precision Kq(MAXZ)
      double precision Sm(MAXZ+1),Sh(MAXZ+1),GH(MAXZ+1)
      double precision bb(MAXZ+1),dd(MAXZ+1),aa(MAXZ+1)
      double precision zetaff,dwm
      integer ilenflag, ibigLflag, iiwmflag

c     ...constants from POM users manual and POM code
c     ...they match Galperin et al., 1988 J. Atmos. Sci. 45:55-62.
      double precision Sq
      parameter(Sq = 0.2) !POM code has 0.41, MY84 suggests 0.2 or 0.5*Sm 16 May 1996

      double precision dtdz2
      double precision du_dz,dv_dz,drho,dUdz
      double precision rlx
      parameter(rlx = 0.8d0)
      double precision NaN,small,smalll
      parameter(NaN=1.d35,small=1.d-8)

      double precision     A1,A2,B1,B2,C1,E1,E2,E3
      parameter(A1=0.92,B1=16.6,A2=0.74,B2=10.1,C1=0.08)
      parameter(E1=1.8,E2=1.33,E3=1.0)
      double precision B123
      parameter(B123=6.50736837) !B1**2/3
      double precision coef1,coef2,coef3,coef4,coef5
      parameter(coef1 = 0.49392771084337D0) ! A2*(1.-6.*A1/B1) (also neutral value)
      parameter(coef2 = 34.6764D0)          ! 3.*A2*B2+18.*A1*A2
      parameter(coef3 = 0.39327228915663D0) ! A1*(1.-3.*C1-6.*A1/B1) (also neutral value)
      parameter(coef4 = 21.3624D0)          ! 18.*A1*A1+9.*A1*A2
      parameter(coef5 = 6.1272D0)           ! 9.*A1*A2
c
c     ...dwm is set to one if the correction for internal-wave shear is
c        requested (by setting iiwmflag = 1).  See Mellor (1989), p. 267.
      dwm = 0.d0
      if( iiwmflag.eq.1) dwm = 1.d0

c     ...calculate Kq at points
      do k=1,nz
         Kq(k) = 0.5d0*Sq*(Km(k)+Km(k+1))
      enddo
      
c     ...q2 equation
      dtdz2 = dt/(dz*dz)
      do k=2,nz
         N2(k) = 0.d0
         Pb(k) = 0.d0
         drho = rho(k)-rho(k-1)
         if( abs(drho) .gt. 1.d-4 )then
            N2(k) =  -betag*(drho)/dz
            Pb(k) =  -Kr(k)*N2(k)
         endif
         du_dz = (dreal(cU(k)-cU(k-1))/dz)**2
         dv_dz = (dimag(cU(k)-cU(k-1))/dz)**2
         dUdz = zetaff*(du_dz+dv_dz+dwm*0.7*N2(k))
         Ps(k) = 0.
         if( abs(dUdz) .gt. 1.d-4 ) Ps(k) =  Km(k)*(dUdz)

         bb(k) = -dtdz2*Kq(k-1)
         dd(k) =  1.+2*dt*q(k)/(B1*l(k))+dtdz2*(Kq(k-1)+Kq(k))
c         dd(k) =  1.+dtdz2*(Kq(k-1)+Kq(k))
         aa(k) = -dtdz2*Kq(k)
         q2(k)  = q(k)**2 + 2.d0*dt*(Ps(k)+Pb(k))
c         q2(k)  = q(k)**2 + 2.d0*dt*(Ps(k)+Pb(k)-q(k)/B1*l(k))
      enddo
c     ...use analytical reference conditions at boundaries
      bb(1)= NaN
      dd(1)= 1.d0
      aa(1)= 0.d0
      q2(1)=B123*ustrb*ustrb
      bb(nz+1)=0.d0
      dd(nz+1)= 1.d0
      aa(nz+1)= NaN
      q2(nz+1)=B123*ustrs*ustrs
      call sy(1,nz+1,bb,dd,aa,q2)

c     ...get rid of small and negative numbers
      q2(1)   = max( q2(1), small )
      q2(nz+1)= max( q2(nz+1), small )
      do k=2,nz
         if(q2(k).lt.small)then
            q2(k)=0.5d0*(q2(k-1)+max(q2(k+1),small))
         endif
      enddo

c     ...master length scale
      bigL(1)   =zob
      bigL(nz+1)=zos
      if(ibigLflag.eq.0)then
c        ...parabolic 
         do k=2,nz
            bigL(k) = 1./( 1./(zf(k)+zob) + 1./(h-(zf(k)+zos)))   !16 May 1996
         enddo
      elseif(ibigLflag.eq.1)then
c        ...bilinear
         do k=2,nz
            bigL(k) = min( h-(zf(k)+zos), zf(k)+zob )             !16 May 1996
         enddo
      else
         stop 'Bad value of ibigLflag in MY25'
      endif

      if(ilenflag.eq.1)then
c     ...q2l equation
         do k=2,nz
            bigW(k) = 1.+E2*(l(k)/(vk*bigL(k)))**2
            bb(k) = -dtdz2*Kq(k-1)
            dd(k) =  1.+(dt*q(k)/(B1*l(k)))*bigW(k)
     &           +dtdz2*(Kq(k-1)+Kq(k))
            aa(k) = -dtdz2*Kq(k)
            q2l(k)  = l(k)*q(k)**2 + dt*E1*l(k)*(Ps(k) + E3*Pb(k)) 
         enddo
c        ...see Melsom (1995) for discussion of these b.c.s
         bb(2)= NaN
         dd(2)= 1.d0
         aa(2)= 0.d0
         q2l(2)=q(2)**2*vk*(zf(2))
         bb(nz)=0.d0
         dd(nz)= 1.
         aa(nz)= NaN
         q2l(nz)=q(nz)**2*vk*(h-zf(nz))
         call sy(2,nz,bb,dd,aa,q2l)
         
         q2l(1)   = max(q(1)**2*vk*zob,small*vk*zob)
         q2l(nz+1)= max(q(nz+1)**2*vk*zos,small*vk*zos)
      endif
      
      if(ilenflag.eq. 0)then
c        ...fixed master length scale (parabolic or bilinear)
         l(1)=vk*zob
         l(nz+1)=vk*zos
         smalll=min(l(1),l(nz+1))
         do k=2,nz
            l(k) = vk*bigL(k)
         enddo
      elseif(ilenflag.eq. 1)then
         l(1)=vk*zob
         l(nz+1)=vk*zos
         smalll=min(l(1),l(nz+1))
         do k=2,nz
            if(q2(k).lt.small)then
               q2(k)=small
c              l(k)=vk*bigL(k)!!changed 11/18/97 to:
               l(k)=smalll
            else
               l(k)=q2l(k)/q2(k)
            endif
c           ...depth scale limit on l
            l(k)=min( l(k), vk*bigL(k) )
c           ...limit l for stable strat (Galperin et al, 1988, eqn. 22)
            if( N2(k).gt. 0.)then
               l(k)=min( 0.53*q(k)/sqrt(N2(k)), l(k) )
            endif
c           ...prevent l = 0
            l(k)=max( l(k), smalll )
         enddo
      else
         stop 'Bad ilenflag in MY25.'
      endif
         
      do k=2,nz
         q(k)=sqrt(q2(k))
c        ...Calculate Richardson number GH at faces
         GH(k) = -(l(k)**2./q2(k))*N2(k)
         GH(k) = min( GH(k), 0.0233d0 )
         GH(k) = max( GH(k),-.28d0)
c        ...calc Sm and Sh according to POM code
         Sh(k)=coef1/(1.-coef2*GH(k))
         Sm(k)=coef3+Sh(k)*coef4*GH(k)
         Sm(k)=Sm(k)/(1.-coef5*GH(k))
c        ...calc new eddy diffusivities
         Kq(k)=Kq(k)+rlx*(nu+Sq*l(k)*q(k)    -Kq(k))
         Km(k)=Km(k)+rlx*(nu+Sm(k)*l(k)*q(k) -Km(k))
         Kh(k)=Kh(k)+rlx*(nu+Sh(k)*l(k)*q(k) -Kh(k))
         Kr(k)=Beta*Km(k)
      enddo
      q(1) = sqrt(q2(1))
      q(nz+1) = sqrt(q2(nz+1))

      Km(1) = Km(1)+rlx*(vk*ustrb*zob+nu-Km(1))
      Kh(1) = Km(1)
      Kr(1) = Beta*Km(1)
      Km(nz+1) = Km(nz+1)+rlx*(vk*ustrs*zos+nu-Km(nz+1))
      Kh(nz+1) = Km(nz+1)
      Kr(nz+1) = Beta*Km(nz+1)
      GH(1)= 0.
      GH(nz+1)= 0.
      Sm(1)=coef1
      Sm(nz+1)=coef1
      Sh(1)=coef3
      Sh(nz+1)=coef3

      return
      end

c***********************************************************************
      subroutine kepsilon(nz,dt,dz,zf,tke,diss,l,Km,Kh,Kr,cU,
     &     rho,Pr,betag,Beta,sigt,vk,nu,h,
     &     ustrs,ustrb,zos,zob,
     &     zetaff,ilenflag,iiwmflag)
c
c     Turbulence closure using k-e model following Rodi (1984)
c     and standard model of Burchard and Baumert (1995)
c     with changes according to Burchard et al (in press)
c 
c***********************************************************************
      implicit none
      include 'onrdp.inc'

      integer  k,nz
      double precision rho(MAXZ),zf(MAXZ+1)
      double complex   cU(MAXZ)
      double precision ustrs,ustrb,vk,nu,betag,Beta
      double precision dt,dz,h,zos,zob
      double precision Ps(MAXZ+1),Pb(MAXZ+1),N2(MAXZ+1)
      double precision Km(MAXZ+1),Kr(MAXZ+1),Kh(MAXZ+1),Pr(MAXZ+1)
      double precision tkeold(MAXZ+1),tke(MAXZ+1)
      double precision diss(MAXZ+1),dissold(MAXZ+1)
      double precision l(MAXZ+1)
      double precision Kp(MAXZ)
      double precision bb(MAXZ+1),dd(MAXZ+1),aa(MAXZ+1)
      double precision zetaff,dwm
      integer ilenflag, iiwmflag

c     ...standard constants from Omstedt, 1983; Rodi, 1984;
c     ...except cmo and c3e from Burchard et al (1998)
      double precision cmo,cmu,c1e,c2e,c3e,sigk,sige,sigt
      parameter(cmo=0.5562,c1e=1.44,c2e=1.92)
      parameter(sigk=1.d0,sige=1.08d0)
      double precision rsqrtcmu, rsigk, rsige
      double precision hb, Pplus,Pminus,Pnet, Ri
      double precision dtdz2
      double precision du_dz,dv_dz,drho,dUdz
c     ...rlx is used to relax from old eddy diffusivities to new values
c     (set to 1.0 to disable)
      double precision rlx,minDiss,minTKE,minD,minDup,minKm
      parameter(rlx = 0.8d0, minTKE=7.6d-6,minDiss=5.d-10)
      double precision maxPr
      parameter(maxPr = 1.6d0) ! approx. equal to Ri = 0.25 when sigt = 1.2
      double precision NaN
      parameter(NaN=1.d35)

c     minKm is the minimum mixing coefficient
      minKm = max(nu,1.3d-6)

c     ...calculate reciprocal of some constants
      cmu = cmo*cmo*cmo*cmo
      rsqrtcmu = 1.d0/(sqrt(cmu))
      rsigk = 1.d0/sigk
      rsige = 1.d0/sige

c     ...hb is a representation of bl thickness in diss. bc
c     (see Rodi, 1984, p. 46)
      hb = min(0.5*h ,5.d0)

c     ...dwm is set to one if the correction for internal-wave shear is
c     requested (by setting iiwmflag = 1).  See Mellor (1989), p. 267.
      dwm = 0.d0
      if( iiwmflag.eq.1) dwm = 1.d0

c     ...convert q to tke (k) for convenience in comparing w/ MY2.5 closure
      do k=1,nz+1
         tke(k)=0.5*tke(k)*tke(k)
      enddo

      do k=1,nz+1
         tkeold(k)=max(tke(k),minTKE)
         dissold(k)=max(diss(k),minDiss)
      enddo

c     ...average eddy viscosity at points
      do k=1,nz
         Kp(k)=0.5*(Km(k)+Km(k+1))
      enddo

c     ...k equation
      dtdz2 = dt/(dz*dz)
      do k=2,nz-1
         N2(k) = 0.d0
         Pb(k) = 0.d0
         Ps(k) = 0.d0
         Pr(k) = sigt
         drho = rho(k)-rho(k-1)
         if( abs(drho) .gt. 1.d-5 )then
            N2(k) =  -betag*(drho)/dz
            Pb(k) =  -Kr(k)*N2(k)
         endif
         du_dz = (dreal(cU(k)-cU(k-1))/dz)**2
         dv_dz = (dimag(cU(k)-cU(k-1))/dz)**2
         dUdz = du_dz+dv_dz
         if( dUdz .gt. 1.d-5 )then
            Ps(k) =  Kp(k)*zetaff*(dUdz+dwm*0.7*N2(k))/sigk
c           ...Prandtl number is f(Ri) (Munk and Anderson, 1948)
c           ...(see Rodi, 1984, eqn. 2.31 or Burchard and Baumaert, 1995)
            Ri = N2(k)/dUdz
            Pr(k) = sigt
            if(Ri .gt. 1.e-4)then
               Pr(k)=min(maxPr,
     &              sigt*((1.+3.333*Ri)**(1.5))/sqrt(1.+10.*Ri))
            endif
         endif

c        ...guarantee positive coefficients (Patankar, 1980, p. 145)
         if(Ps(k)+Pb(k).gt.0.d0)then
            Pplus  = Ps(k)+Pb(k)
            Pminus = dissold(k)
         else
            Pplus  = Ps(k)
            Pminus = dissold(k)-Pb(k)
         endif
         bb(k)=-dtdz2*rsigk*Kp(k-1)
         dd(k)= 1.+dtdz2*rsigk*(Kp(k-1)+Kp(k))+Pminus*dt/tkeold(k)
         aa(k)=-dtdz2*rsigk*Kp(k)
         tke(k)= tkeold(k)+dt*Pplus
      enddo

c     ...surface boundary condition
      if( tkeold(nz).gt.rsqrtcmu*ustrs*ustrs )then
c        ...no-flux of tke
         bb(nz)= -dtdz2*rsigk*Kp(nz)
         dd(nz)= 1.+dtdz2*rsigk*Kp(nz)
         aa(nz)= NaN
         tke(nz)=tkeold(nz)
      else
c        ...Rodi(1984) eqn. 2.88
         bb(nz)=0.d0
         dd(nz)= 1.d0
         aa(nz)= NaN
         tke(nz)=rsqrtcmu*ustrs*ustrs
      endif
c     ...bottom boundary condition
      bb(1)= NaN
      dd(1)= 1.d0
      aa(1)= 0.d0
      tke(1)=rsqrtcmu*ustrb*ustrb

c     ...tridiagonal solver
      call sy(1,nz,bb,dd,aa,tke)

c     ...dissipation (epsilon) equation
      do k=2,nz-1

c        ...c3e from Bauchart et al (in press)
c        (this value is more suspect than most)
         if(Pb(k).gt.0.)then
            c3e = 1.d0
         else
            c3e =-0.4d0
         endif
                                                                     
c        ...guarantee positive coefficients (Patankar, 1980, p. 145,
c           or Burchard et al p. 11).
         if(c3e*Pb(k).gt.0.d0)then
            Pplus =(dissold(k)/tkeold(k))*(c1e*Ps(k)+c3e*Pb(k))
            Pminus=(dissold(k)/tkeold(k))*c2e*dissold(k)
         else
            Pplus =(dissold(k)/tkeold(k))*c1e*Ps(k)
            Pminus=(dissold(k)/tkeold(k))*(c2e*dissold(k)-c3e*Pb(k))
         endif

         bb(k) = -dtdz2*rsige*Kp(k-1)
         dd(k) =  1.+dtdz2*rsige*(Kp(k-1)+Kp(k))+Pminus*dt/dissold(k)
         aa(k) = -dtdz2*rsige*Kp(k)
         diss(k)  = dissold(k) + dt*Pplus
      enddo

c     ...bottom bc from Rodi (1984), eqn. 2.87
      bb(1)= NaN
      dd(1)= 1.d0
      aa(1)= 0.d0
      diss(1) = ustrb**3/(vk*zob)

      bb(nz)=0.d0
      dd(nz)= 1.
      aa(nz)= NaN
      diss(nz) = ustrs**3/(vk*zos)
c____________________________________________________________________
c     ...free surface bc from Rodi (1984), eqn. 2.89 causes probs
c     (note that value of hb is not clear, and large values decrease diss)
c      kscmu=tkeold(nz)*sqrt(cmu)
c      diss(nz)= kscmu**(3./2.)
c     &     /(vk*((h-zf(nz))+0.07*hb*(1.-ustrb**2/kscmu)))
c____________________________________________________________________
c     ...tridiagonal solver
      call sy(1,nz,bb,dd,aa,diss)

c     ...get rid of small and negative numbers
      tke(1) = max( tke(1), minTKE )
      tke(nz+1)= max( tke(nz+1), minTKE )
      do k=2,nz
         if(tke(k).lt.minTKE)then
            tke(k)=0.5d0*(tke(k-1)+max(tke(k+1),minTKE))
         endif
      enddo
      diss(1) = max( diss(1), minDiss )
      diss(nz+1)= max( diss(nz+1), minDiss )
      do k=2,nz
         minD   = max( minDiss, 0.045*tke(k)*tke(k)*N2(k) )
         minDup = max( minDiss, 0.045*tke(k+1)*tke(k+1)*N2(k+1) )
         if(diss(k).lt.minD)then
            diss(k)=0.5d0*(diss(k-1)+max(diss(k+1),minDup))
         endif
      enddo

c     ...calculate diffusivities
      Pr(1) = Pr(2)
      Km(1) = Km(1)+rlx*(max(vk*ustrb*zob,minKm)-Km(1))
      Kh(1) = Km(1)/Pr(1)
      Kr(1) = Beta*Km(1)
      l(1)  = vk*zob
      do k=2,nz
         Km(k)=Km(k)+rlx*(max(minKm,cmu*tke(k)*tke(k)/diss(k))-Km(k))
         Kh(k)=Kh(k)+rlx*(Km(k)/Pr(k) - Kh(k))
         Kr(k)=Beta*Km(k)
         l(k) =cmo**3 * tke(k)**(3./2.)/diss(k)
      enddo
      Pr(nz+1) = Pr(nz)
      Km(nz+1) = Km(nz+1)+rlx*(max(vk*ustrs*zos,minKm)-Km(nz+1))
      Kh(nz+1) = Km(nz+1)/Pr(nz+1)
      Kr(nz+1) = Beta*Km(nz+1)
      l(nz+1)  = vk*zos

      if(.false.)then
         do k=nz+1,1,-5
            Pnet = Pb(k)+Ps(k)-diss(k)
            write(80,1000)
     &           zf(k),Km(k),tke(k),diss(k),
     &           Pb(k),Ps(k),Pnet,l(k),Pr(k)
 1000       format(1x,f4.1,8(1x,f12.4))
         enddo
      endif

c     ...convert tke (k) to q for convenience in comparing w/ MY2.5 closure
      do k=1,nz+1
         tke(k)=sqrt(2.*tke(k))
      enddo
 
      return
      end
c***********************************************************************

      subroutine advect( nz, dt, dz, u, P)
c
c     Implementation of van Leer advection scheme for
c     constant velocity and uniform grid spacing
c     Assumes no-flux boundary conditions
c     Advects mass p with velocity u.
c     Updates p with new values.
c
c     Allen, Douglass, Rood, and Guthrie (1991)
c     Eqns. 2.2 - 2.4
c
      implicit none
      include 'onrdp.inc'
      integer nz
      double precision dt, dz, u, P(MAXZ)
      double precision flux(MAXZ+1)! flux into point i at z(i)-dz/2
      double precision Pold(MAXZ)
      double precision dtdz,slw,sle,delCw,delCe
      integer i

      if(u .eq. 0.) return

      do i=1,nz
         Pold(i)=P(i)
      enddo

      dtdz = dt/dz

      flux(1)=0.                ! no flux lower bc

      if(u.gt.0.)then
         flux(2) = u*P(1)       ! upwind approx.
         do i=3,nz-1
            slw = (P(i-1)-P(i-2))*(P(i)-P(i-1))
            delCw = 0.
            if(slw.gt.0. .and. abs(P(i)-P(i-2)).gt. 0.0001)
     &           delCw = 2.*slw/(P(i)-P(i-2))
            flux(i)=
     &           u*(P(i-1)+0.5*(1.-u*dtdz)*delCw)
         enddo
         flux(nz)=u*P(nz-1)
      elseif(u.lt.0.)then
         do i=2,nz-1
            sle = (P(i)-P(i-1))*(P(i+1)-P(i))
            delCe = 0.
            if(sle.gt.0. .and. abs(P(i+1)-P(i-1)).gt. 0.0001)
     &           delCe = 2.*sle/(P(i+1)-P(i-1))
            flux(i)=
     &           u*(P(i)-0.5*(1.+u*dtdz)*delCe)
         enddo
         flux(nz) = u*P(nz)     ! upwind approx.
      else
         stop 'huh?'
      endif
      flux(nz+1) = 0. ! no flux upper bc

      do i=1,nz
         P(i)=Pold(i)-dtdz*(flux(i+1)-flux(i))
      enddo
      return
      end
c***********************************************************************
c
      subroutine diffus(il,ih, xp, xf, D, P, dt, Fl, Fh )
c
c     Diffusion with specified flux boundary conditions
c     Flux formulation, fully implicit solution
c
c     Supply values of C at points listed in zp,
c     an array with length n=(kh-kl)+1, and values of D at
c     interfaces listed in zf, an array with length n+1.
c     zf(kl) = lower boundary; zf(kh+1) = upper boundary.
c
c     Specify flux Fl at bottom and Fh at top
c
c     Result is returned in C
c
c***********************************************************************
      implicit none
      include 'onrdp.inc'
      integer il, ih
      double precision dt, xp(*), xf(*), D(*), P(*), Fl, Fh
      double precision bb(MAXZ), dd(MAXZ), aa(MAXZ)

      integer i
      double precision dxh,dxl,dxm
      double precision NaN /1.e35/

c     ...bottom bc
      i=il
      dxh = xp(i+1)-xp(i)
      dxm = xf(i+1)-xf(i)
      bb(i) =  NaN
      dd(i) =  dt*D(i+1)/(dxm*dxh) + 1.
      aa(i) = -dt*D(i+1)/(dxm*dxh)
      P(i) = P(i)+2.*dt*Fl/(dxm)

c     ...interior points
      do i=il+1,ih-1
         dxl = xp(i)-xp(i-1)
         dxh = xp(i+1)-xp(i)
         dxm = xf(i+1)-xf(i)
         bb(i) = -dt*D(i)/(dxm*dxl)
         dd(i) =  dt*D(i)/(dxm*dxl)+dt*D(i+1)/(dxm*dxh)+1
         aa(i) = -dt*D(i+1)/(dxm*dxh)
      enddo

c     ...upper bc
      i=ih
      dxl = xp(i)-xp(i-1)
      dxm = xf(i+1)-xf(i)
      bb(i) = -dt*D(i)/(dxm*dxl)
      dd(i) =  dt*D(i)/(dxm*dxl) + 1
      aa(i) =  NaN
      P(i) = P(i)+2.*dt*Fh/(dxm)

c     ...tridiagonal solver
      call sy(il,ih,bb,dd,aa,P)
      return
      end
c***********************************************************************
c
      double precision function we_calc(fr,tauc,taucsfmx,c_we)
c
c     Calculate maximum erosion rate [m/s]
c
c*********************************************************************
      implicit none
      double precision fr, tauc, taucsfmx, c_we(3)
      double precision xS
      xS = max(0.d0,(taucsfmx-tauc)/tauc)
      we_calc = c_we(1)*fr*xS**c_we(2)
      return
      end
c***********************************************************************

      subroutine integ( nz, dz, zos, zob, cU, xtran, ytran )
c
c     Calculate depth-integrated transport
c
      implicit none
      include 'onrdp.inc'
      integer nz,i
      double precision dz
      double precision zos, zob
      double complex cU(MAXZ)
      double precision xtran, ytran
      xtran=0.
      ytran=0.
      xtran = xtran+(dz-zob)*dreal(cU(1))
      ytran = ytran+(dz-zob)*dimag(cU(1))
      do i=2,nz-1
         xtran = xtran+dz*dreal(cU(i))
         ytran = ytran+dz*dimag(cU(i))
      enddo
      xtran = xtran+(dz-zos)*dreal(cU(nz))
      ytran = ytran+(dz-zos)*dimag(cU(nz))
      return
      end
c***********************************************************************

      subroutine conserve(nz,dz,nice,nsed,S,Tw,Ci,Cs,Li,cpw,rhobar,
     &     heat,salt,sedvol,icevol)
c
c     Calculate depth-integrated quantities
c
      implicit none
      include 'onrdp.inc'
      integer nz,nice,nsed
      double precision dz
      double precision S(MAXZ)          ! salinity [psu]
      double precision Tw(MAXZ)         ! water temp
      double precision Ci(MAXZ,NICEMAX) ! concentration of ice 
      double precision Cs(MAXZ,NSEDMAX) ! concentration of sediment [0 to 1]
      double precision Li, cpw, rhobar

      double precision heat,salt,sedvol,icevol
      integer j,k
      double precision totC,totS
      
      heat = 0.d0
      salt = 0.d0
      sedvol = 0.d0
      icevol = 0.d0
      do k=1,nz
         totC = 0.d0
         do j=1,nice
            totC = totC+Ci(k,j)
         enddo
         totS = 0.d0
         do j=1,nsed
            totS = totS+Cs(k,j)
         enddo
         icevol = icevol+dz*totC
         sedvol = sedvol+dz*totS
         salt = salt+dz*S(k)*(1.-totC)
         heat = heat+dz*(Tw(k)+273.)*cpw*rhobar*(1.-totC)
     &        + dz*totC*rhobar*Li
      enddo
      return
      end
c***********************************************************************

      subroutine iflux( nz, dz, zos, zob, cU, Ci, j, xtran, ytran )
c
c     Calculate depth-integrated transport of ice class j
c
      implicit none
      include 'onrdp.inc'
      integer nz,i,j
      double precision dz
      double precision zos, zob
      double complex cU(MAXZ)
      double precision Ci(MAXZ,NICEMAX)
      double precision xtran, ytran
      xtran=0.
      ytran=0.
      xtran = xtran+(dz-zob)*dreal(cU(1))*Ci(1,j)
      ytran = ytran+(dz-zob)*dimag(cU(1))*Ci(1,j)
      do i=2,nz-1
         xtran = xtran+dz*dreal(cU(i))*Ci(i,j)
         ytran = ytran+dz*dimag(cU(i))*Ci(i,j)
      enddo
      xtran = xtran+(dz-zos)*dreal(cU(nz))*Ci(nz,j)
      ytran = ytran+(dz-zos)*dimag(cU(nz))*Ci(nz,j)
      return
      end
c***********************************************************************

      subroutine sflux( nz, dz, zos, zob, cU, Cs, j, xtran, ytran )
c
c     Calculate depth-integrated transport of sed class j
c
      implicit none
      include 'onrdp.inc'
      integer nz,i,j
      double precision dz
      double precision zos, zob
      double complex cU(MAXZ)
      double precision Cs(MAXZ,NSEDMAX)
      double precision xtran, ytran
      xtran=0.
      ytran=0.
      xtran = xtran+(dz-zob)*dreal(cU(1))*Cs(1,j)
      ytran = ytran+(dz-zob)*dimag(cU(1))*Cs(1,j)
      do i=2,nz-1
         xtran = xtran+dz*dreal(cU(i))*Cs(i,j)
         ytran = ytran+dz*dimag(cU(i))*Cs(i,j)
      enddo
      xtran = xtran+(dz-zos)*dreal(cU(nz))*Cs(nz,j)
      ytran = ytran+(dz-zos)*dimag(cU(nz))*Cs(nz,j)
      return
      end
c***********************************************************************
c     
      subroutine bstress(Hs,Td,wdir,u,v,zf,ks,h,rhow,
     &    zoa,Cd,tausfmx)
c
c     Calculate bottom shear stress
c
c     Input:
c        Hs - significant wave height [m]
c        Td - dominant wave period [s]
c        u,v - velocity components at elevation zr [m/s]
c        zf - height of bottom box [m]
c        ks - bottom grain roughness [m]
c        h - water depth (positive number) [m]
c        rhow - water density [kg/m3]
c
c     Output:
c        Cd - drag coefficient
c        tausfmx - magnitude of max. skin friction shear stress [N m-2]
c
      implicit none
      double precision Hs,Td,wdir,u,v,zf,ks,Cd,h,rhow,tausfmx
      
      double precision ub_qf
      integer iverbose
      double precision ustrc,ustrw,ustrcw,fwc
      double precision ub, ur, phic, udir, wr
      double precision rtd
      double precision vk, zo, zoa, zr,zrold
      parameter(vk = 0.41)
      parameter(rtd = 57.29577951)
      integer MAXIT,nit
      parameter(MAXIT=5)
      iverbose = 0

c     ...use previous zoa or zo, whichever is larger
      zo = ks/30.
      zoa = max(zo,zoa)

c     ...zr is chosen at the level where ur = <u>
      zr = zf*0.3679*exp(zoa/zf)
      Cd = (vk/(log(zr/zoa)))**2

      call pcoord( u, v, ur, udir )
      tausfmx = rhow*Cd*ur*ur
      ub = ub_qf(Td,Hs,h)
      if (ub .le. 0.01) return
      
      phic = wdir-udir
      wr = 6.2831853/Td
      nit = 0
      zrold = 999.
      do while( abs(zrold-zr).gt. 0.1 .and. nit .lt. MAXIT )
         nit = nit+1
         call madsen94(ub, wr, ur, zr, phic,
     &        ks, iverbose,
     &        ustrc, ustrw, ustrcw, fwc, zoa )
         zrold = zr
         zr = zf*0.3679*exp(zoa/zf)
         tausfmx = rhow*(ustrcw)**2
         Cd = (vk/(log(zr/zoa)))**2
      enddo
      end
c***********************************************************************
c
      subroutine madsen94( ubr, wr, ucr,
     &    zr, phiwc, kN, iverbose,
     &    ustrc, ustrwm, ustrr, fwc, zoa )
c
c   Grant-Madsen model
c   from Madsen(1994)
c   input:
c      ubr = rep. wave-orbital velocity amplitude outside wbl [m/s]
c      wr = rep. angular wave frequency = 2pi/T [rad/s]
c      ucr = current velocity at height zr [m/s]
c      zr = reference height for current velocity [m]
c      phiwc = angle between currents and waves at zr (radians)
c      kN = bottom roughness height (e.q. Nikuradse k) [m]
c      iverbose = switch; when 1, extra output
c   returned:
c      ustrc  = current friction velocity         u*c [m/s]
c      ustrr  = w-c combined friction velocity    u*r [m/s]
c      ustrwm = wave max. friction velocity      u*wm [m/s]
c      fwc = wave friction factor [ ]
c      zoa = apparent bottom roughness [m]
c
c**********************************************************************
      implicit none
      double precision ubr, ucr
      double precision zr, phiwc, kN
      double precision ustrc, ustrwm, phicwc
      double precision ustrr, fwc
      integer iverbose
      integer MAXIT
      parameter(MAXIT = 20)
      double precision PI, VK
      parameter(PI = 3.14159265)
      parameter(VK = 0.4)
      double precision zo, zoa, wr
      double precision rmu(MAXIT), Cmu(MAXIT), fwci(MAXIT), dwc(MAXIT)
      double precision ustrwm2(MAXIT), ustrr2(MAXIT), ustrci(MAXIT)
      double precision cosphiwc,lnzr,lndw,lnln,bigsqr,diff
      integer i, nit

      double precision fwc94 ! function

c     ...junk return values
      ustrc = 99.99
      ustrwm = 99.99
      ustrr = 99.99
      fwc = .4
      zoa = kN/30.
      phicwc = phiwc

c     ...some data checks
      if( wr .le. 0. ) then
	 write(*,1000) wr
 1000	 format(1x,
     &   'WARNING: Bad value for frequency in Madsen94: wr=',f6.4)
         return
      endif
      if( ubr .lt. 0. ) then
	 write(*,1010) ubr
 1010	 format(1x,
     &   'WARNING: Bad value for orbital vel. in Madsen94: ub=',f6.2)
         return
      endif
      if( kN .lt. 0. ) then
         write(*,1020) kN
 1020	 format(1x,
     &   'WARNING: Wierd value for roughness in Madsen94: kN=',f7.3)
         return
      endif
      if( (zr.lt.zoa .or. zr.lt.0.05).and. iverbose .eq. 1)then
         write(*,1030) zr
 1030	 format(1x,
     &   'WARNING: Low value for ref. level in Madsen94: zr=',f8.3)
      endif

      zo = kN/30.
  
      if(ubr .le. 0.01)then
         if(ucr .le. 0.01)then
c           ...no waves or currents
            ustrc = 0.
            ustrwm = 0.
            ustrr = 0.
            return
         endif
c        ...no waves
         ustrc = ucr * VK / log(zr/zo) 
         ustrwm = 0.
         ustrr = ustrc
         return
      endif
  
      cosphiwc =  abs(cos(phiwc))
      rmu(1) = 0.
      Cmu(1) = 1.
      fwci(1) = fwc94( Cmu(1), (Cmu(1)*ubr/(kN*wr)) ) !Eqn. 32 or 33
      ustrwm2(1)= 0.5*fwci(1)*ubr*ubr                 !Eqn. 29
      ustrr2(1) = Cmu(1)*ustrwm2(1)                   !Eqn. 26
      ustrr = sqrt( ustrr2(1) )
      dwc(1) = kN
      if ((Cmu(1)*ubr/(kN*wr)) .ge. 8.) dwc(1)= 2.*VK*ustrr/wr
      lnzr = log(zr/dwc(1))
      lndw = log(dwc(1)/zo)
      lnln = lnzr/lndw
      bigsqr = (-1.+sqrt(1+ ((4.*VK*lndw)/(lnzr*lnzr))*ucr/ustrr ))
      ustrci(1) = 0.5*ustrr*lnln*bigsqr
      nit = 1
      diff = 1.
      do i=2,MAXIT
         rmu(i) = ustrci(i-1)*ustrci(i-1)/ustrwm2(i-1)
         Cmu(i) = sqrt(1.+2.*rmu(i)*cosphiwc+rmu(i)*rmu(i))!Eqn 27
         fwci(i) = fwc94( Cmu(i), (Cmu(i)*ubr/(kN*wr)) )   !Eqn. 32 or 33
         ustrwm2(i)= 0.5*fwci(i)*ubr*ubr                   !Eqn. 29
         ustrr2(i) = Cmu(i)*ustrwm2(i)                     !Eqn. 26
         ustrr = sqrt( ustrr2(i) )
         dwc(i) = kN
         if ((Cmu(1)*ubr/(kN*wr)) .ge. 8.) dwc(i)= 2.*VK*ustrr/wr !Eqn.36
         lnzr = log( zr/dwc(i) )
         lndw = log( dwc(i)/zo )
         lnln = lnzr/lndw
         bigsqr = (-1.+sqrt(1+ ((4.*VK*lndw)/(lnzr*lnzr))*ucr/ustrr ))
         ustrci(i) = 0.5*ustrr*lnln*bigsqr                  !Eqn. 38
         diff = abs( (fwci(i)-fwci(i-1))/fwci(i) )
	 if(diff .lt. 0.0005) goto 100
         nit = nit+1
      enddo
 100  ustrwm = sqrt( ustrwm2(nit) )
      ustrc = ustrci(nit)
      ustrr = sqrt( ustrr2(nit) )
      phicwc = phiwc
      zoa = exp( log(dwc(nit))-(ustrc/ustrr)*log(dwc(nit)/zo) ) !Eqn. 11
      fwc = fwci(nit)
      if(iverbose .eq. 1)then
        do i=1,nit
           write(*,1050)
     &        i,fwci(i),dwc(i),ustrci(i),sqrt(ustrwm2(i)),
     &        sqrt(ustrr2(i))
 1050	   format(1x,'i=',i2,' fwc=',f9.6,', dwc=',f9.6,', u*c=',f9.4,
     &        ', u*wm=',f9.4,', u*r=',f9.4)
         enddo
      endif
      return
      end
c***********************************************************************

      double precision function fwc94( cmu, cukw )
c
c     Wave-current friction factor
c     Equations 32 and 33 in Madsen, 1994
c
      double precision cmu, cukw
      double precision fwc
      fwc = .00999 !meaningless (small) return value
      if( cukw .le. 0. )then
	 write(*,1000) cukw
 1000	 format(1x,'ERROR: cukw too small in fwc94: ',f9.4)
         return
      endif
      if( cukw .lt. 0.2 )then
        fwc = exp( 7.02*0.2**(-0.078) - 8.82 )
	write(*,1010)cukw
 1010	format(1x,'WARNING: cukw very small in fwc94: ',f9.4)
      endif
      if( cukw .gt. 0.2 .and. cukw .le. 100. )
     &   fwc = exp( 7.02*cukw**(-0.078)-8.82 )
      if( cukw .gt. 100. .and. cukw .le. 10000. )
     &   fwc = exp( 5.61*cukw**(-0.109)-7.30 )
      if( cukw .gt. 10000.)then
        fwc = exp( 5.61*10000.**(-0.109)-7.30 )
	write(*,1020) cukw
 1020	format(1x,'WARNING: cukw very large in fwc94: ',f9.4)
      endif
      fwc94 = cmu*fwc
      end
c***********************************************************************

      double precision function coriolis( rlat )
c
c   Coriolis parameter
c     Input: latitude [decimal degrees]
c     Returned: Coriolis frequency [s-1]
c 
      double precision omega, rlat
      omega = 2.*3.14159625/86164             ! 2*pi/sideral day (s)
      coriolis = 2.*omega*sin( rlat/57.29577951 )   ! f [radians/s]
      return
      end
c***********************************************************************
c
      subroutine csy(il,iu,bb,dd,aa,cc)
c
c  Solves the tridiagonal system of equations using Thomas' algorithm
c  (Anderson et al. pp. 549-550)
c
c  Double Complex version
c  csy solves by elimination:
c  il = subscript of first equation
c  iu = subscript of last equation
c  bb = coefficient behind diagonal (unchanged)
c  dd = coefficient on diagonal (changed on return)
c  aa = coefficient ahead of diagonal (unchanged)
c  cc = element of constant vector on input, solution vector on return
c
      integer il,iu
      double complex aa(*),bb(*),cc(*),dd(*)
      integer i,j,lp
      double complex r
c...establish upper triangular matrix
      lp = il+1
      do 10 i = lp,iu
         r = bb(i)/dd(i-1)
         dd(i) = dd(i) -r*aa(i-1)
         cc(i) = cc(i) -r*cc(i-1)
 10   continue
c...back substitution
      cc(iu) = cc(iu)/dd(iu)
      do 20 i = lp,iu
         j = iu-i+il
         cc(j) = (cc(j)-aa(j)*cc(j+1))/dd(j)
 20   continue
c...solution stored in cc
      return
      end
c***********************************************************************
c
      subroutine sy(il,iu,bb,dd,aa,cc)
c
c  Solves the tridiagonal system of equations using Thomas' algorithm
c  (Anderson et al. pp. 549-550)
c
c  SY solves by elimination
c
c
c  il = subscript of first equation
c  iu = subscript of last equation
c  bb = coefficient behind diagonal (unchanged)
c  dd = coefficient on diagonal (changed on return)
c  aa = coefficient ahead of diagonal (unchanged)
c  cc = element of constant vector on input, solution vector on return
c
c***********************************************************************
c-------SCCS time and version stamp:
c
c	@(#)sy.F	1.1	93/02/10
c

      integer iu,il
      double precision aa(*),bb(*),cc(*),dd(*)
      integer i,j,lp
      double precision r
c...establish upper triangular matrix
      lp = il+1
      do 10 i = lp,iu
         r = bb(i)/dd(i-1)
         dd(i) = dd(i) -r*aa(i-1)
         cc(i) = cc(i) -r*cc(i-1)
 10   continue
c...back substitution
      cc(iu) = cc(iu)/dd(iu)
      do 20 i = lp,iu
         j = iu-i+il
         cc(j) = (cc(j)-aa(j)*cc(j+1))/dd(j)
 20   continue
c...solution stored in cc
      return
      end

c********************************************************************

      subroutine readinit(initfile, isrestart, rho_off, 
     &     nice, nsed, rhoi, rhos, sigt, zob, zos, vk, kconst, beta,
     &     nz, h, z, zf, Temp, Sal, cU, Ci, Cs, rho, Km, Kr, Kh,
     &     Tfri, cT, q, l, Rin, Rfn, Sm, Sh )
c********************************************************************
      implicit none
      include 'onrdp.inc'
      character*40 initfile
      integer nz, isrestart, nice, nsed
      double precision h,z(MAXZ)
      double precision Temp(MAXZ),Sal(MAXZ)
      double complex cU(MAXZ)
      double precision rho(MAXZ)
      double precision zf(MAXZ+1),Km(MAXZ+1),Kr(MAXZ+1),Kh(MAXZ+1)
      double precision Ci(MAXZ,NICEMAX),Cs(MAXZ,NSEDMAX)
      double precision Tfri(MAXZ)
      double complex cT(MAXZ+1)
      double precision Sm(MAXZ+1),Sh(MAXZ+1),Rin(MAXZ+1),Rfn(MAXZ+1)
      double precision l(MAXZ+1),q(MAXZ+1)
      double precision rhoi,rhos,sigt,zob,zos,vk,kconst,beta
      double precision frac,fraci,fracs,rhow

      integer MAXZIN
      parameter(MAXZIN=100)
      double precision zin(MAXZIN),Tempin(MAXZIN),Salin(MAXZIN)
      double precision Kin(MAXZIN),uin(MAXZIN),vin(MAXZIN)
      double precision Ciin(MAXZIN,NICEMAX), Csin(MAXZIN,NICEMAX)
      double precision rhoin(MAXZIN),Tfrin(MAXZIN)
      double precision cTin(MAXZIN),qin(MAXZIN),lin(MAXZIN)
      double precision Rinin(MAXZIN),Smin(MAXZIN),Shin(MAXZIN)
      double precision rho_off
      integer i,j,m,n,incount,nin
      double precision u,v
c     ...functions
      double precision tfreeze,density2

      open(55,file=initfile,status='old')
      incount=0
      do i=1,MAXZIN-2
         if(isrestart.eq.0)then
            read(55,*,end=99) zin(i),Tempin(i),Salin(i),
     &           Kin(i),uin(i),vin(i),
     &           (Ciin(i,n),n=1,nice),(Csin(i,m),m=1,nsed)
         elseif(isrestart.eq.1)then
            read(55,*,end=99) zin(i), Tempin(i), Salin(i),
     &           Kin(i),uin(i),vin(i),
     &           (Ciin(i,n),n=1,nice),(Csin(i,m),m=1,nsed),
     &           rhoin(i),Tfrin(i),cTin(i),qin(i),lin(i),
     &           Rinin(i),Smin(i),Shin(i)
         else
            stop 'readinit: bad value of isrestart.'
         endif
         incount=incount+1
      enddo
 99   write(*,1000) incount,initfile
 1000 format(1x,'Read ',i6,' records from ',a)
      close(55)
      if(incount.lt.1) stop 'not enough data in initfile'

      nin=incount

c     ...move array down one index
      do j=incount+1,2,-1
         zin(j) = zin(j-1)
         Tempin(j)=Tempin(j-1)
         Salin(j)=Salin(j-1)
         Kin(j)=Kin(j-1)
         uin(j)=uin(j-1)
         vin(j)=vin(j-1)
         do m=1,nice
            Ciin(j,m)=Ciin(j-1,m)
         enddo
         do m=1,nsed
            Csin(j,m)=Csin(j-1,m)
         enddo
         if(isrestart.eq.1)then
            rhoin(j)=rhoin(j-1)
            Tfrin(j)=Tfrin(j-1)
            cTin(j)=cTin(j-1)
            qin(j)=qin(j-1)
            lin(j)=lin(j-1)
            Rinin(j)=Rinin(j-1)
            Smin(j)=Smin(j-1)
            Shin(j)=Shin(j-1)
         endif
      enddo
c     ...add surface value (may be redundant)
      zin(1)=h
      Tempin(1)=Tempin(2)
      Salin(1)=Salin(2)
      uin(1)=uin(2)
      vin(1)=vin(2)
      Kin(1)=Kin(2)
      do m=1,nice
         Ciin(1,m)=Ciin(2,m)
      enddo
      do m=1,nsed
         Csin(1,m)=Csin(2,m)
      enddo
      if(isrestart.eq.1)then
         rhoin(1)=rhoin(2)
         Tfrin(1)=Tfrin(2)
         cTin(1)=cTin(2)
         qin(1)=qin(2)
         lin(1)=lin(2)
         Rinin(1)=Rinin(2)
         Smin(1)=Smin(2)
         Shin(1)=Shin(2)
      endif

      nin=incount+1
c     ...make sure bottom value is at bottom
      if( zin(incount+1).gt. 0.) then
         zin(incount+2)=0.
         Tempin(incount+2)=Tempin(incount+1)
         Salin(incount+2)=Salin(incount+1)
         uin(incount+2)=0.
         vin(incount+2)=0.
         Kin(incount+2)=0.
         do m=1,nice
            Ciin(incount+2,m)=Ciin(incount+1,m)
         enddo
         do m=1,nsed
            Csin(incount+2,m)=Csin(incount+1,m)
         enddo
         if(isrestart.eq.1)then
            rhoin(incount+2)=rhoin(incount+1)
            Tfrin(incount+2)=Tfrin(incount+1)
            cTin(incount+2)=cTin(incount+1)
            qin(incount+2)=qin(incount+1)
            lin(incount+2)=lin(incount+1)
            Rinin(incount+2)=Rinin(incount+1)
            Smin(incount+2)=Smin(incount+1)
            Shin(incount+2)=Shin(incount+1)
         endif
         nin=nin+1
      endif
      
c     ...interpolate to model grid
      j=nin-1
      do i=1,nz
         do while( zin(j) .lt. z(i) .and. j .ge. 1)
            j=j-1
         enddo
         frac = (z(i)-zin(j+1))/( zin(j)-zin(j+1) )
         Temp(i)=Tempin(j+1)+frac*(Tempin(j)-Tempin(j+1))        
         Sal(i)=Salin(j+1)+frac*(Salin(j)-Salin(j+1))       
            u=uin(j+1)+frac*(uin(j)-uin(j+1))       
            v=vin(j+1)+frac*(vin(j)-vin(j+1))
            cU(i)=cmplx(u,v)     
            Km(i)=Kin(j+1)+frac*(Kin(j)-Kin(j+1))
            do m=1,nice
               Ci(i,m)=Ciin(j+1,m)+frac*(Ciin(j,m)-Ciin(j+1,m))
            enddo
            do m=1,nsed
               Cs(i,m)=Csin(j+1,m)+frac*(Csin(j,m)-Csin(j+1,m))
            enddo
         if(isrestart.eq.1)then
            rho(i)=rhoin(j+1)+frac*(rhoin(j)-rhoin(j+1))       
            Tfri(i)=Tfrin(j+1)+frac*(Tfrin(j)-Tfrin(j+1))       
            cT(i)=cmplx(cTin(j+1)+frac*(cTin(j)-cTin(j+1)),0.)      
            q(i)=qin(j+1)+frac*(qin(j)-qin(j+1))       
            l(i)=lin(j+1)+frac*(lin(j)-lin(j+1))       
            Rin(i)=Rin(j+1)+frac*(Rin(j)-Rin(j+1))       
            Sm(i)=Smin(j+1)+frac*(Smin(j)-Smin(j+1))       
            Sh(i)=Shin(j+1)+frac*(Shin(j)-Shin(j+1))       
         endif
      enddo

      if(isrestart.ne.1)then
c        ...velocity
         do i=1,nz
            cU(i)=cmplx(0.,0.)
         enddo

c        ...freezing point of water
         do i=1,nz
            Tfri(i)=tfreeze( Sal(i), 0.1*(h-z(i)) )
         enddo

c        ...calculate density
         do i=1,nz
            rhow = density2(Sal(i),Temp(i))
            fraci = 0.
            do j=1,nice
               fraci = fraci + Ci(i,j)
            enddo
            fracs = 0.
            do j=1,nsed
               fracs = fracs + Cs(i,j)
            enddo
            rho(i)=fraci*(rhoi-rho_off)
     &           +fracs*(rhos-rho_off)
     &           +(1.-(fraci+fracs))*(rhow-rho_off)
         enddo

c        ...eddy viscosity
         l(1)=zob
         l(nz+1)=zos
         do i=2,nz
            l(i)=vk*1./(1./zf(i)+1./(h-zf(i)))
         enddo

         do i=1,nz+1
            Km(i)=Kconst
            q(i)=0.d0
            Rin(i) = 0.d0
            Sm(i) = 0.406d0 !neutral
            Sh(i) = 0.507d0 !neutral
            Rfn(i) = 0.
            Kr(i)=Beta*Km(i)
            Kh(i)=Km(i)/sigt
         enddo
      endif

      return
      end

c***********************************************************************

      subroutine readforce(forcefile,nt,t,wu,wv,Ta,snow,dsdx,dsdy,
     &     ixshelf)
c
c     Reads input files for forcing parameters and linearly
c     interpolates to model time
c
c     The ixshelf switch is for backward compatibility with onrforce.dat
c     files that do not have slopes specified.
c
      implicit none
      include 'onrdp.inc'

      integer MAXTIN
      parameter(MAXTIN=1000)
      integer nt
      integer ixshelf
      double precision t(MAXT), Ta(MAXT), wu(MAXT), wv(MAXT)
      double precision snow(MAXT), dsdx(MAXT), dsdy(MAXT)
      character*40 forcefile
      integer i,j,incount
      double precision tin(MAXTIN),wuin(MAXTIN)
      double precision wvin(MAXTIN),Tain(MAXTIN)
      double precision snowin(MAXTIN)
      double precision dsdxin(MAXTIN),dsdyin(MAXTIN)
      double precision frac
      open(61,file=forcefile,status='old')
      incount=0
      if(ixshelf.ne.0)then 
c        ...read input files with no slope data
         do i=1,MAXTIN
            read(61,*,end=99) tin(i),wuin(i),wvin(i),Tain(i),
     &           snowin(i)
            incount=incount+1
         enddo
      endif
      if(ixshelf.eq.0)then
c        ...read input files with slope data
         do i=1,MAXTIN
            read(61,*,end=99) tin(i),wuin(i),wvin(i),Tain(i),
     &           snowin(i),dsdxin(i),dsdyin(i)
            incount=incount+1
         enddo
      endif
 99   write(*,1000) incount,forcefile
 1000 format(1x,'Read ',i6,' records from ',a)
      close(61)
      if(tin(1) .gt. t(1)) stop 'Input time starts too late'
      if(tin(incount) .lt. t(nt) ) stop 'input time ends too soon'
      
c     ...interpolate to model time
      j=2
      do i=1,nt
         do while( tin(j) .lt. t(i) )
            j=j+1
         enddo
         frac = (t(i)-tin(j-1))/( tin(j)-tin(j-1) )
         wu(i)=wuin(j-1)+frac*(wuin(j)-wuin(j-1))        
         wv(i)=wvin(j-1)+frac*(wvin(j)-wvin(j-1))        
         Ta(i)=Tain(j-1)+frac*(Tain(j)-Tain(j-1))        
         snow(i)=snowin(j-1)+frac*(snowin(j)-snowin(j-1))
         if(ixshelf.eq.0)then
            dsdx(i)=dsdxin(j-1)+frac*(dsdxin(j)-dsdxin(j-1))        
            dsdy(i)=dsdyin(j-1)+frac*(dsdyin(j)-dsdyin(j-1))        
         endif
      enddo
      return
      end

c***********************************************************************

      subroutine makeslope(nt,t,xamp,yamp,xlag,ylag,dsdx,dsdy)
      include 'onrdp.inc'
      integer nt
      double precision t(MAXT),dsdx(MAXT),dsdy(MAXT)
      double precision xamp,yamp,xlag,ylag
      double precision TWOPI
      parameter(TWOPI=6.28318531)
      integer n

      per = 12.4*3600.
      freq = TWOPI/per
      do n=1,nt
         dsdx(n) = xamp*cos(t(n)*freq + TWOPI*xlag/360.)
         dsdy(n) = yamp*cos(t(n)*freq + TWOPI*ylag/360.)
      enddo
      return
      end

c***********************************************************************

      subroutine wstress(n, windx, windy, cTs )

c..Calculate the surface wind stress components [N m-2] from the wind 
c..components provided [m s-1]

      implicit none
      integer i,n
      double precision windx(*),windy(*),tausx,tausy,airdens,Cd,speed
      double complex cTs(*)

      airdens = 1.3 ![kg m-3]

c     ...Cd from Large and Pond, 1981.

      do i=1,n
         speed = sqrt(windx(i)*windx(i) + windy(i)*windy(i))
         Cd = 1.140
         if (speed.gt.10.0) Cd = 0.49 + 0.065 * speed
         if (speed.gt.26.0) Cd = 2.18
         Cd = Cd * 0.001
         tausx = Cd * airdens * windx(i) * speed
         tausy = Cd * airdens * windy(i) * speed
         cTs(i)=cmplx( tausx, tausy )
      enddo
      return
      end
      subroutine areadfil3(nin,anot,nchar,dvar,avar,ifail)
C **************************************************************************
C *   FUNCTION    :  Reads CHARACTER*40 variable from file based on keyword. *
C *             (Modified version of areadfil, without output to operator *
C *                or log file)                                           *
C *                                                                       *
C     nin ..... File input device                               *
C     anot .... Keyword in file                                 *
C     nchar ... Number of characters in ANOT (max. 40)          *
C     dvar .... Default value for CHARACTER*40 variable         *
C     avar .... Resultant CHARACTER*40 variable
C     ifail
c       input: 0 = stop if unsuccessful
c              1 = warn if unsuccessful
c              2 = silent
c       return: 0 for successful execution, otherwise 1         *
c       
C *                                                                          *
C ****************************************************************************
      integer*4 nin,nchar,ifail
      character*(*) anot
      character*40 avar
      character*40 dvar
C
      integer*4 ios
      logical found
      character*80 buff
C
      rewind(nin)
      ios=0
      found=.false.
      do while(ios.eq.0.and..not.found)
        read(nin,1,iostat=ios) buff
    1   format(a80)
        if(buff(1:nchar).eq.anot.and.
     $     buff(nchar+1:nchar+1).eq.' ') then  ! Match has been found
          avar(1:40)=buff(nchar+2:nchar+41)
          found=.true.
        endif
      end do
      if(ios.eq.0) then
        ifail=0                                ! Match found and no error
        return
      else
        if(ifail.eq.0)then
           write(*,1000) anot
 1000      format(1x,'Fatal error reading ',a)
           stop
        elseif(ifail.eq.1)then
           ifail = 1
           write(*,1200) anot
 1200      format(1x,'Warning: Could not read ',a)
           avar = dvar
           write(*,1300) avar
 1300      format(1x,'Using default: ',a)
           return
        elseif(ifail.eq.2)then
           ifail = 1
           avar = dvar
           return
        else
           stop 'Bad value of ifail passed to areadfil3.'
        endif
      endif
      end
c******************************************************************************
c
      subroutine ireadfil3(nin,anot,nchar,dvar,ifail,ivar)
c
c     Reads INTEGER variable from file based on keyword.
c     Adapted from John Hunter's ireadfil2.
c
c     Input:
c        nin ..... File input device
c        anot .... Keyword in file
c        nchar ... Number of characters in ANOT (max. 40)
c        dvar  ... Default value to use
c        ifail ... 0 = Stop if not found
c                  1 = Warn and use default if not found
c                  2 = Silent and use default if not found
c     Returns:
c        ivar .... Resultant INTEGER variable
c        
c******************************************************************************
      integer nin,nchar,ivar,ifail,dvar
      character*(*) anot
      integer ios
      logical found
      character*80 buff
      logical isbot

      rewind(nin)
      ios=0
      found=.false.
      do while(ios.eq.0.and..not.found)
        read(nin,1000,iostat=ios) buff
 1000   format(a80)
        if(buff(1:nchar).eq.anot.and.
     &    isbot( buff(nchar+1:nchar+1) )) then     ! Match has been found
          read(buff(nchar+2:80),*,iostat=ios) ivar
          found=.true.
        endif
      end do
      if(ios.eq.0) then
         write(*,1400)anot,ivar
 1400    format(1x,a,4x,i7)
        return                  ! Match found and no error
      else
        if(ifail.eq.0)then
           write(*,1100) anot
 1100      format(1x,'Fatal error reading ',a)
           stop
        elseif(ifail.eq.1)then
           write(*,1200) anot
 1200      format(1x,'Warning: Could not read ',a)
           ivar = dvar
           write(*,1300) ivar
 1300      format(1x,'Using default value of: ',i8)
           return
        elseif(ifail.eq.2)then
           ivar = dvar
           return
        else
           stop 'Bad value of ifail passed to ireadfil3.'
        endif
      endif
      end
c******************************************************************************
c
      subroutine dreadfil3(nin,anot,nchar,dvar,ifail,var)
c
c     Reads double precision variable from file based on keyword.
c     Adapted from John Hunter's ireadfil2.
c
c     Input:
c        nin ..... File input device
c        anot .... Keyword in file
c        nchar ... Number of characters in ANOT (max. 40)
c        ifail ... 0 = Stop if not found
c                  1 = Warn and use default if not found
c                  2 = Silent and use default if not found
c     Returns:
c        var .... Resultant double precision variable
c        
c******************************************************************************
      integer nin,nchar,ifail
      double precision dvar,var
      character*(*) anot
      integer ios
      logical found
      character*80 buff
      logical isbot

      rewind(nin)
      ios=0
      found=.false.
      do while(ios.eq.0.and..not.found)
        read(nin,1000,iostat=ios) buff
 1000   format(a80)
        if(buff(1:nchar).eq.anot.and.
     &    isbot( buff(nchar+1:nchar+1) )) then     ! Match has been found
          read(buff(nchar+2:80),*,iostat=ios) var
          found=.true.
        endif
      end do
      if(ios.eq.0) then
         write(*,1400)anot,var
 1400    format(1x,a,4x,g14.8)
         return                 ! Match found and no error
      else
         if(ifail.eq.0)then
            write(*,1100) anot
 1100       format(1x,'Fatal error reading ',a)
            stop
         elseif(ifail.eq.1)then
            write(*,1200) anot
 1200       format(1x,'Warning: Could not read ',a)
            var = dvar
            write(*,1300) var
 1300       format(1x,'Using default value of: ',g12.6)
            return
         elseif(ifail.eq.2)then
            var = dvar
            return
         else
            stop 'Bad value of ifail passed to dreadfil3.'
         endif
      endif
      end
c****************************************************************
c
      logical function isbot( c )
c
c     Returns true if character is blank or tab
c
      character*1 c
      character*1 b,t
      b = char(32)
      t = char(9)
      isbot = ((c .eq. b) .or. (c .eq. t))
      return
      end
c***********************************************************************
c
      subroutine waves(nt, h, fetch, wu, wv, Hs, Tm, ub)
c
c***********************************************************************
      implicit none
      include 'onrdp.inc'
      integer nt
      double precision h, fetch
      double precision wu(MAXT), wv(MAXT), Hs(MAXT)
      double precision Tm(MAXT), ub(MAXT)
      integer i
      double precision ws, w2, Uafunc, ub_qf
      do i=1,nt
         w2 = ( wu(i)*wu(i) + wv(i)*wv(i) )
         if(w2 .gt. 0.1 .and. fetch .gt. 0. ) then
            ws = sqrt( w2 )
            call shallow( Uafunc(ws),fetch,h,Hs(i),Tm(i))
            ub(i)=ub_qf(Tm(i),Hs(i),h)
         else
            Hs(i) = 0.
            Tm(i) = 0.5
            ub(i) = 0.
         endif
      enddo
      return
      end
c***********************************************************************
c
      subroutine shallow( Ua,F,h,Hm,Tm )
c
c
c     Shallow-water wave height, period
c     Shore Protection Manual, eqns. 3-39 and 3-40
c     @(#)shallow.F	1.1 29 Apr 1996
c
c***********************************************************************
      double precision Ua,F,h,Hm,Tm
      double precision g
      parameter(g = 9.81)
      double precision a,b,gd,gf
      Tm = .5
      Hm = 0.
      if(Ua.le.0.) return

      gd = g*h/(Ua*Ua)
      gf = g*F/(Ua*Ua)
      a = 0.283*tanh(0.530*gd**(3./4.))*tanh( (0.00565*sqrt(gf))/
     &     tanh(0.530*gd**(3./4.)))
      Hm = a*Ua*Ua/g
      b = 7.54*tanh(0.833*gd**(3./8.))*tanh( (0.0379*gf**(1./3.))/
     &     tanh(0.833*gd**(3./8.)))
      Tm = b*Ua/g
      return
      end
c***********************************************************************

      double precision function Uafunc( U )
c     Calculates adjusted wind speed from wind windspeed [m/s]
c     Shore Protection Manual Equation 3-28a.
c     @(#)uafunc.F	1.1 29 Apr 1996
      double precision U
      Uafunc = 0.71 * U**(1.23)
      return
      end

c***********************************************************************
c
      double precision function ub_qf(T,Hs,h)
c
c  Quick version for calc. near-bottom orbital velocity
c  
c  Input: T  = sig. wave period
c         Hs = sig. wave height
c         h  = water depth
c
c  cf. Dyer(1986) eqn. 3.50, p. 98
c  Modified 15 May 1998 to avoid divide by huge kh
c
c  Returns: ub_qf = Hs*pi / T*sinh(kh)
c
c***********************************************************************
      implicit none
      double precision T,Hs,h,w,kh
      double precision twopi,g
      double precision khf
      parameter(g=9.80665,twopi=6.28318531)

      ub_qf = 0.

      if( T .le. 0.01 ) then
         return
      elseif( Hs .le. 0.001 ) then
         return
c      elseif( h .le. Hs*10.) then
c         write(*,*) 'Warning from ub_qf: linear theory violated.'
      elseif( h .le. Hs*1.28 ) then
         write(*,*) 'Warning from ub_qf: wave is breaking.'
      endif
          
      w = twopi/T
      kh = khf(T,h)
      if(Hs .gt. .01 .and. kh .lt. 10.)
     &     ub_qf = w * Hs/(2.*sinh(kh))      
      return
      end
c******************************************************************
c
      double precision function khf(T,h)
c
c     Quick version of kh in dispersion relationship
c     T is wave period in seconds
c     h is depth in meters
c     Dean and Dalrymple (1991), p. 72
c
      implicit none
      double precision T,h
      double precision y,kh2,w
      double precision g,twopi,d1,d2,d3,d4,d5,d6
      parameter(  g=9.80665,
     & twopi=6.28318531,
     & d1=0.6666666666,
     & d2=0.3555555555,
     & d3=0.1608465608,
     & d4=0.0632098765,
     & d5=0.0217540484,
     & d6=0.0065407983 )
      w=twopi/T
      y=(w**2)*h/g
      kh2=y**2+y/(1+d1*y+d2*y**2+d3*y**3+d4*y**4+d5*y**5+d6*y**6)
      khf=sqrt(kh2)
      return
      end
c***********************************************************************
c
        subroutine pcoord(x,y,r,az)
c
c       Converts u,v or x,y to polar coordiantes on 0-360 circle
c
c       Written by: Chris Sherwood, NORTEC, UW, and Battelle
c
c*********************************************************************
        double precision x,y,r,az
        double precision rtdeg,alpha
        double precision DPI /3.141592653589793d0/
        double precision rtemp, aztemp
        rtdeg=180.000d0/DPI
        rtemp=sqrt(x**2+y**2)
C       ...Prevent division by zero
        if( y.ne.0.) then        
           alpha=rtdeg*atan(abs(x/y))
        else if(x.ge.0.) then
           aztemp=90.d0
           goto 100 
        else 
           aztemp=270.d0
           goto 100 
        endif 
        if(x .ge. 0.d0) goto 10
          if(y .ge. 0.d0) goto 20               
              aztemp=180.+alpha                  
              goto 100 
 20        aztemp=360.-alpha                    
           goto 100 
 10    if(y .ge. 0.d0) goto 30                  
           azemp=180.-alpha   
           goto 100 
 30    aztemp=+alpha                            
 100   continue
       az=aztemp
       r=rtemp
       return
       end
c*********************************************************************

      double precision function tfreeze( S, P )
      double precision S, P
      
c   Freezing point of water (deg C) at salinity S (psu)
c   and pressure P (bars).  Millero (1978) cited in
c   Gill (1982) p. 602.
      
      tfreeze = -0.0575*S + 1.710523e-3*sqrt(S*S*S)
     &     - 2.154996e-4*S*S - 7.53e-3*P
      return
      end

c*********************************************************************
      double precision function density2(S,T)
c
c  Density of sea water using quadratic equation of Omstedt
c  and Svensson (1984) Eqn. 16.
c
c  Input:
c    S Salinity [psu]
c    T Temperature [deg C] 
c  Returns:
c    density2 Density [kg m-3]
c
      double precision S, T
      double precision rho0, alpha, beta, Tm
      parameter(rho0=1000.,alpha =5.6e-6,beta=8.0e-4,Tm = 2.9)
      density2 = rho0*(1.-alpha*(T-Tm)**2+beta*S)
      return
      end
