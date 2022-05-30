! time evolution
program tevol_ensvsa

  use read_netcdf
  
  implicit none
 
  ! Parameters
#ifdef lev6
  integer,parameter :: nlev=6 
#else
  integer,parameter :: nlev=3
#endif
  integer,parameter :: nv2d=2 !ps, PE(ps)
#ifdef moist
  integer,parameter :: nv3d=4 !u,v,T,q
#else
  integer,parameter :: nv3d=3 !u,v,T
#endif
  integer,parameter :: nvar=nv3d*nlev+nv2d-1
  double precision,parameter :: pi=atan(1.0d0)*4.0d0
  double precision,parameter :: cp=1005.7d0, R=287.04d0
#ifdef moist
  double precision,parameter :: Lh=2.5104d6
#ifndef weak
  double precision,parameter :: epsilon=1.0d0
#else
  double precision,parameter :: epsilon=1.0d-1
#endif
#endif 
  double precision,parameter :: Tr=270.0d0, pr=1.0d3
  ! Target region
  double precision, parameter :: slon=137.0d0, elon=142.0d0, slat=33.0d0, elat=37.0d0
  integer :: islon, ielon, islat, ielat
  !!integer,parameter :: islon=95, ielon=105, islat=67, ielat=75 
  !integer,parameter :: islon=275, ielon=285, islat=247, ielat=255 
  !!integer,parameter :: islon=271, ielon=281, islat=247, ielat=255 
  !!integer,parameter :: islon=273, ielon=279, islat=229, ielat=235 
  integer :: nlon, nlat, narea
       
  real :: lat, area
  integer :: i,j
  integer :: imem,id,it,irec
  integer :: ilt,ilu,ilv,ilq,ilon,ilat,ilev,ivar,ip,fday,tday
  integer :: imode
  !integer :: mem
  !integer :: idate,edate
  logical :: ex
  
  ! 2D variable
  double precision,allocatable ::  ps(:,:)
  ! 3D variable
  double precision,allocatable :: &
  &                    ug(:,:,:),vg(:,:,:),&
  &                    T(:,:,:),gh(:,:,:)
#ifdef moist
  double precision,allocatable :: q(:,:,:)!,rh(:,:,:)
#endif
  ! perturbation arrays
  real,allocatable ::  zv3(:,:,:),zv(:,:)
  double precision,allocatable ::  z0(:,:,:),zm(:,:,:)
  double precision,allocatable ::  ze(:,:,:,:)
  double precision,allocatable ::  z(:,:),zT(:,:)

  ! level
  double precision ::  sigma(nlev),ssg
  double precision :: plev(nlev),dplev(nlev-1)
#ifdef lev6
   data plev/300.0d0,500.0d0,700.0d0,850.0d0,925.0d0,1000.0d0/
#else
   data plev/300.0,500.0,850.0/
#endif
  real,allocatable :: sg(:),p(:),w(:,:)
  double precision,allocatable :: TE(:),KE(:),PE(:),LE(:)
  ! Namelist
  character(len=5) :: orig="jma"
  integer :: mem=26 
  integer :: idate=2019100912
  integer :: edate=2019101212
  integer :: smode=1
  integer :: emode=1
  namelist /sens_nml/ orig, mem, idate, edate, smode, emode

  character(len=100) rdf,rdw,wd(4)
  character dir*32!,dira*33
  character ns*1,ne*1,nmem*2,cnlev*1
  character yyyy*4,mm*2,mmddhh*6,yyyymmddhh*10,cedate*10
  character(len=3) :: vname(5)
  !character(len=4) :: vnamea(5)
  !data vname/'UGRD','VGRD','TMP','SPFH','PRES_meansealevel'/
  data vname/'u','v','t','q','msl'/
  !data vnamea/'air','uwnd','vwnd','shum','slp'/
     !|----/----/----/----/----/----/----/----/----/----| 
  !dir='/Users/nakashita/netcdf/tigge/'
  dir='/Volumes/dandelion/netcdf/tigge/'
 !dira='/Users/nakashita/netcdf/nc-reanl/'

   !sigma(1)=200.0/pr
   !sigma(2)=6.0/7.0*300.0/pr
   !sigma(3)=8.0/7.0*300.0/pr
   do ilev=1,nlev-1
      dplev(ilev) = plev(ilev+1) - plev(ilev)
   end do
   sigma = 0.0
   do ilev=1,nlev-1
      sigma(ilev) = sigma(ilev) + 0.5*dplev(ilev)/pr
      sigma(ilev+1) = sigma(ilev+1) + 0.5*dplev(ilev)/pr 
   end do
   print*,sigma
      
  !  データの設定 
   open(11,file="sens.nml")
   read(11,nml=sens_nml)
   write(*,nml=sens_nml)
   close(11)
   !call calc_steps(idate,edate,12,tday)
   call calc_steps(idate,edate,6,tday)
   print*,tday
   tday = tday + 4
   !yyyymmddhh="2019100912"
   write(yyyymmddhh,'(I10)') idate
   yyyy=yyyymmddhh(1:4)
   mm=yyyymmddhh(5:6)
   mmddhh=yyyymmddhh(5:10)
   print*,yyyy,mmddhh
   write(cedate,'(I10)') edate
   print*, cedate
   !smode=1
   !emode=1
   write(ns,'(I1)') smode
   write(ne,'(I1)') emode
   write(cnlev,'(I1)') nlev
!
   call set_bound(slon,elon,slat,elat,islon,ielon,islat,ielat)
   print*,islon,ielon,islat,ielat
   nlon=ielon-islon+1
   nlat=ielat-islat+1 
   narea=nlon*nlat
! File Open 
   if (smode==emode) then
#ifdef moist
#ifndef weak
      wd(1)='TE-moist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(2)='KE-moist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(3)='PE-moist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(4)='LE-moist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
#else
      wd(1)='TE-wmoist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(2)='KE-wmoist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(3)='PE-wmoist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(4)='LE-wmoist-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
#endif
#else
      wd(1)='TE-dry-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(2)='KE-dry-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(3)='PE-dry-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
#endif
   else
#ifdef moist
#ifndef weak
      wd(1)='TE-moist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(2)='KE-moist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(3)='PE-moist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(4)='LE-moist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
#else
      wd(1)='TE-wmoist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(2)='KE-wmoist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(3)='PE-wmoist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(4)='LE-wmoist-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
#endif
#else
      wd(1)='TE-dry-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(2)='KE-dry-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
      wd(3)='PE-dry-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'.txt'
#endif
   endif
   open(21,file=wd(1),status='replace')
   open(22,file=wd(2),status='replace')
   open(23,file=wd(3),status='replace')
#ifdef moist
   open(24,file=wd(4),status='replace')
#endif
   !mem=memn 
#ifdef moist
#ifndef weak
   rdw='./weight-TE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.grd'
#else
   rdw='./weight-wTE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.grd'
#endif
#else
   rdw='./weight-dTE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.grd'
#endif
   open(10,file=rdw,status='old',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*mem)

  ! 配列の割付
   allocate(&
   ps(imax,jmax),ug(imax,jmax,nlev),vg(imax,jmax,nlev)&
   ,T(imax,jmax,nlev)&
   )
#ifdef moist
   allocate(q(imax,jmax,nlev))!,rh(imax,jmax,nlev))
#endif
   allocate(&
   zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)&
   ,z0(nlon,nlat,nvar),zm(nlon,nlat,nvar)&
   )
   allocate(ze(nlon,nlat,nvar,mem))
   allocate(z(narea*nvar,mem))
   allocate(zT(mem,narea*nvar))
   allocate(sg(smode:emode))
   allocate(p(mem))
   allocate(w(mem,smode:emode))
   allocate(TE(mem+1),KE(mem+1),PE(mem+1))
#ifdef moist
   allocate(LE(mem+1))
#endif
  !singular value
   read(10,rec=1) p
   do imode=smode,emode
      sg(imode)=p(imode)
   enddo
   print*, sg
  
   do imode=smode,emode
      it=imode+1
      read(10,rec=it) w(:,imode)
      print*,imode,w(:,imode)
   enddo
  
   close(10)

   irec=1
   !edate=2019100912
   read(yyyymmddhh, *) edate
   print*,"edate=",edate
   do fday=0,tday !every 6 hours
      !idate=2019100912
      read(yyyymmddhh,*) idate
      call calc_steps(idate,edate,6,ip)
      ip=ip+1
      print *, "ip=",ip
      ze = 0.0d0
      do imem=1,mem
         write(nmem,'(I2.2)') imem
         !print*,nmem
         ilt=0
         ilu=0
         ilv=0
         ilq=0
         !rdf=dir//yyyy//'/jma/'//mmddhh//'_'//nmem//'.nc'
         rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_'//nmem//'.nc'
         !print*,rdf
         inquire(file=rdf, exist=ex)
         if(ex)then
            do id=1,3
               !print*,rdf
               !print*,vname(id)
               call fread3(rdf,vname(id),ip,zv3)
               !call fread3a(rdf,vnamea(id),ip,zv3,90.0d0,180.0d0,0.0d0,80.0d0)
               !print*,maxval(zv),minval(zv)
               if    (id==1)then
               !ug=zv3(:,:,1:3)
#ifdef lev6
               ug=dble(zv3(:,:,4:))
#else
               ug=dble(zv3(:,:,3:))
#endif
            elseif(id==2)then
               !vg=zv3(:,:,1:3)
#ifdef lev6
               vg=dble(zv3(:,:,4:))
#else
               vg=dble(zv3(:,:,3:))
#endif
               elseif(id==3)then
               !T=zv3(:,:,1:3)
#ifdef lev6
               T =dble(zv3(:,:,4:))
#else
               T =dble(zv3(:,:,3:))
#endif
            endif
         enddo
#ifdef moist
         call fread3(rdf,vname(4),ip,zv3)
#ifdef lev6
         q =dble(zv3(:,:,4:))
#else
         q =dble(zv3(:,:,3:))
#endif
#endif
            call fread(rdf,vname(5),ip,zv)
            ps=dble(zv)*1.0d-2     !Pa->hPa
            !print*,ps(1,1)
           
            do i=1,nlon
               do j=1,nlat
                  ilon=islon+i-1
                  ilat=islat+j-1
                  ze(i,j,1:nlev,imem)=ug(ilon,ilat,:)
                  ze(i,j,nlev+1:2*nlev,imem)=vg(ilon,ilat,:)
                  ze(i,j,2*nlev+1:3*nlev,imem)=T(ilon,ilat,:)
#ifdef moist
                  ze(i,j,3*nlev+1:4*nlev,imem)=q(ilon,ilat,:)
#endif
                  ze(i,j,nvar,imem)=ps(ilon,ilat)
                  !ze(i,j,10,imem)=ps(ilon,ilat)
               enddo
            enddo
            !! check
            write(6,'(A,I2)') 'mem ',imem
            write(6,'(A,'//cnlev//'f12.4)') 'u ',ze(1,1,1:nlev,imem)
            write(6,'(A,'//cnlev//'f12.4)') 'v ',ze(1,1,nlev+1:2*nlev,imem)
            write(6,'(A,'//cnlev//'f12.4)') 't ',ze(1,1,2*nlev+1:3*nlev,imem)
#ifdef moist
            write(6,'(A,'//cnlev//'es12.4)') 'q ',ze(1,1,3*nlev+1:4*nlev,imem)
#endif
            write(6,'(A,f12.4)') 'ps',ze(1,1,nvar,imem)
         endif
      enddo
     
      !idate=2019100900 !nomark,_a
      !call calc_steps(idate,edate,6,ip) !_a->12,nomark&_n->6
      !ip=ip+1
      !print *, "ip=",ip
      !ip=fday+2
      !print*,ip
      ilt=0
      ilu=0
      ilv=0
      ilq=0
      z0 = 0.0d0
      !rdf=dir//yyyy//'/jma/'//mmddhh//'_mean.nc'   !_n
      rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_mean.nc'   !_n
      !rdf=dir//yyyy//'/jma/100900_mean.nc'
      !rdf=dir//yyyy//'/jma/anl_sellev.nc' !_a
      inquire(file=rdf, exist=ex)
      if(ex)then
         do id=1,3
            call fread3(rdf,vname(id),ip,zv3)
            if    (id==1)then
            !ug=zv3(:,:,1:3)
#ifdef lev6
            ug=dble(zv3(:,:,4:))
#else
            ug=dble(zv3(:,:,3:))
#endif
         elseif(id==2)then
            !vg=zv3(:,:,1:3)
#ifdef lev6
            vg=dble(zv3(:,:,4:))
#else
            vg=dble(zv3(:,:,3:))
#endif
            elseif(id==3)then
            !T=zv3(:,:,1:3)
#ifdef lev6
            T =dble(zv3(:,:,4:))
#else
            T =dble(zv3(:,:,3:))
#endif
         endif
      enddo
#ifdef moist
      call fread3(rdf,vname(4),ip,zv3)
      !q=zv3(:,:,1:3)
#ifdef lev6
      q =dble(zv3(:,:,4:))
#else
      q =dble(zv3(:,:,3:))
#endif
      !rh=zv3(:,:,1:3)
      !print*,"rh",rh(1,1,1)
      !do k=1,kmax
      !   do j=1,jmax
      !      do i=1,imax
      !         call calc_q(T(i,j,k),rh(i,j,k),plev(k),q(i,j,k))
      !      enddo
      !   enddo
      !enddo
#endif
         call fread(rdf,vname(5),ip,zv)
         ps=dble(zv)*1.0d-2     !Pa->hPa
         !print*,ps(1,1)
         !print*,maxval(ps),minval(ps)
        
         do i=1,nlon
            do j=1,nlat
               ilon=islon+i-1
               ilat=islat+j-1
               z0(i,j,1:nlev)=ug(ilon,ilat,:)
               z0(i,j,nlev+1:2*nlev)=vg(ilon,ilat,:)
               z0(i,j,2*nlev+1:3*nlev)=T(ilon,ilat,:)
#ifdef moist
               z0(i,j,3*nlev+1:4*nlev)=q(ilon,ilat,:)
#endif
               z0(i,j,nvar)=ps(ilon,ilat)
               !z0(i,j,10)=ps(ilon,ilat)
            enddo
         enddo
      endif
      !! check
      write(6,'(A,'//cnlev//'f12.4)') 'u ',z0(1,1,1:nlev)
      write(6,'(A,'//cnlev//'f12.4)') 'v ',z0(1,1,nlev+1:2*nlev)
      write(6,'(A,'//cnlev//'f12.4)') 't ',z0(1,1,2*nlev+1:3*nlev)
#ifdef moist
      write(6,'(A,'//cnlev//'es12.4)') 'q ',z0(1,1,3*nlev+1:4*nlev)
#endif
      write(6,'(A,f12.4)') 'ps',z0(1,1,nvar)

      !1.calcurate perturbation
      do imem=1,mem
         ze(:,:,:,imem)=ze(:,:,:,imem)-z0(:,:,:)
      enddo

      !2.multiply by cos(lat) and layer thickness factor
      ssg=0.0
      do imode=smode,emode
         ssg=ssg+sg(imode)
      enddo

      area=0.0
      do i=1,nlon
         do j=1,nlat
            lat=(islat+j-2)*dlat - 90.0
            area=area+cos(lat*pi/180.0)
            ze(i,j,:,:)=ze(i,j,:,:)*sqrt(cos(lat*pi/180.0))
         enddo
      enddo
      area=area*nlon
      print*, area
      ze = ze / sqrt(area)
     
      do ilev=1,nlev            !300,500,700,850,925,1000hPa
         do ivar=1,3         !ug,vg,T
            ze(:,:,nlev*(ivar-1)+ilev,:)=ze(:,:,nlev*(ivar-1)+ilev,:)*sqrt(sigma(ilev))
         enddo
#ifdef moist
         ivar=4 !q
         ze(:,:,nlev*(ivar-1)+ilev,:)=ze(:,:,nlev*(ivar-1)+ilev,:)*sqrt(sigma(ilev))
#endif
      enddo
     !3.Multiply by coefficient
     !T
      ze(:,:,2*nlev+1:3*nlev,:)=ze(:,:,2*nlev+1:3*nlev,:)*sqrt(cp/Tr)
     !ps
      ze(:,:,nvar,:)=ze(:,:,nvar,:)*sqrt(R*Tr)/pr
     !ze(:,:,10,:)=ze(:,:,10,:)*sqrt(R*Tr)/pr
#ifdef moist
     !q
      ze(:,:,3*nlev+1:4*nlev,:)=ze(:,:,3*nlev+1:4*nlev,:)*Lh/sqrt(cp*Tr)*sqrt(epsilon)
#endif
      !do imem=1,mem
      !   print*,imem
      !   print*,ze(1,1,:,imem)
      !enddo

     !4.calcurate energy
     !TE:1~nvar, KE=1~6, PE=7~9+13, LE=10~12
      TE=0.0
      KE=0.0
      PE=0.0
#ifdef moist
      LE=0.0
#endif
      zm=0.0
      do imode=smode,emode
         do imem=1,mem
            zm=zm+ze(:,:,:,imem)*w(imem,imode)*sg(imode)/ssg
         enddo
      enddo
     
      do ivar=1,nvar
         do j=1,nlat
            do i=1,nlon
               TE(1)=TE(1)+zm(i,j,ivar)**2
            enddo 
         enddo
      enddo
      TE(1)=TE(1)*0.5
      do ivar=1,2*nlev
         do j=1,nlat 
            do i=1,nlon
               KE(1)=KE(1)+zm(i,j,ivar)**2
            enddo
         enddo
      enddo
      KE(1)=KE(1)*0.5
      do ivar=2*nlev+1,3*nlev
         do j=1,nlat 
            do i=1,nlon
               PE(1)=PE(1)+zm(i,j,ivar)**2
            enddo
         enddo
      enddo
      do j=1,nlat
         do i=1,nlon
            PE(1)=PE(1)+zm(i,j,nvar)**2
         enddo
      enddo
      PE(1)=PE(1)*0.5
#ifdef moist
      do ivar=3*nlev+1,4*nlev
         do j=1,nlat 
            do i=1,nlon
               LE(1)=LE(1)+zm(i,j,ivar)**2
            enddo
         enddo
      enddo
      LE(1)=LE(1)*0.5
#endif
      ! member
      do imem=1,mem
         do ivar=1,nvar
            do j=1,nlat
               do i=1,nlon
                  TE(imem+1)=TE(imem+1)+ze(i,j,ivar,imem)**2
               enddo 
            enddo
         enddo
         TE(imem+1)=TE(imem+1)*0.5
         do ivar=1,2*nlev
            do j=1,nlat 
               do i=1,nlon
                  KE(imem+1)=KE(imem+1)+ze(i,j,ivar,imem)**2
               enddo
            enddo
         enddo
         KE(imem+1)=KE(imem+1)*0.5
         do ivar=2*nlev+1,3*nlev
            do j=1,nlat 
               do i=1,nlon
                  PE(imem+1)=PE(imem+1)+ze(i,j,ivar,imem)**2
               enddo
            enddo
         enddo
         do j=1,nlat
            do i=1,nlon
               PE(imem+1)=PE(imem+1)+ze(i,j,nvar,imem)**2
            enddo
         enddo
         PE(imem+1)=PE(imem+1)*0.5
#ifdef moist
         do ivar=3*nlev+1,4*nlev
            do j=1,nlat 
               do i=1,nlon
                  LE(imem+1)=LE(imem+1)+ze(i,j,ivar,imem)**2
               enddo
            enddo
         enddo
         LE(imem+1)=LE(imem+1)*0.5
#endif
      enddo
      write(nmem,'(I2.2)') mem+1
      write(21,'(f5.1,'//nmem//'e16.8)') real(fday*6), TE
      write(22,'(f5.1,'//nmem//'e16.8)') real(fday*6), KE
      write(23,'(f5.1,'//nmem//'e16.8)') real(fday*6), PE
#ifdef moist
      write(24,'(f5.1,'//nmem//'e16.8)') real(fday*6), LE
#endif
      idate=edate
      !call calc_date(idate,2,edate)
      call calc_date(idate,1,edate)
      print*, "edate=",edate
   enddo
   close(21)
   close(22)
   close(23)
#ifdef moist
   close(24)
#endif
   deallocate(zv3,zv,z0,ug,vg,T,ps)
#ifdef moist
   deallocate(q)
#endif
   deallocate(ze,z,zT,sg,p,w)
  
   stop  
    
end program tevol_ensvsa