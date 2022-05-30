program grads_ensvsa

  use read_netcdf
  
  implicit none

  ! Parameters
#ifdef lev6
  integer,parameter :: nlev=8,nlev_te=6 
#else
  integer,parameter :: nlev=5,nlev_te=3 
#endif
  integer,parameter :: nv2d=2 !ps, PE(ps)
#ifdef moist
  integer,parameter :: nv3d=5+4 !u,v,T,q,gh + TE,KE,PE(T),LE
  integer,parameter :: nvar=(nv3d-4)*nlev+nv2d-1
#else
  integer,parameter :: nv3d=4+3 !u,v,T,gh + TE,KE,PE(T)
  integer,parameter :: nvar=(nv3d-3)*nlev+nv2d-1
#endif
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

  integer :: i,j,k,n
  integer :: imem,id,it,irec
  integer :: ilt,ilu,ilv,ilq,ilev,ivar,ip,fday,tday
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
  double precision :: plev(nlev),dplev(nlev)
#ifdef lev6
  data plev/200.0d0,250.0d0,300.0d0,500.0d0,700.0d0,850.0d0,925.0d0,1000.0d0/
#else
  data plev/200.0d0,250.0d0,300.0d0,500.0d0,850.0d0/
#endif
  ! weights
  real,allocatable :: sg(:),p(:),w(:,:)
  ! diagnostic variables
  double precision,allocatable :: TE(:,:,:),KE(:,:,:),PE(:,:,:)
#ifdef moist
  double precision,allocatable :: LE(:,:,:)
#endif
  real,allocatable :: v3d(:,:,:,:) !ug,vg,T,q,gh,TE,KE,PE(T),LE
  real,allocatable :: v2d(:,:,:)   !ps,PE(ps)
  real,allocatable :: buf4(:,:)
! Namelist
  character(len=5) :: orig="jma"
  integer :: mem=26 
  integer :: idate=2019100912
  integer :: edate=2019101212
  integer :: smode=1
  integer :: emode=1
  namelist /sens_nml/ orig, mem, idate, edate, smode, emode

  character rdf*100,rdw*100,wd*100,wds*100
  character dir*32!,dira*33
  character ns*1,ne*1,nmem*2,cnlev*1
  character yyyy*4,mm*2,mmddhh*6,yyyymmddhh*10,cedate*10
  character(len=3) :: vname(6)
  !character(len=4) :: vnamea(5)
  !data vname/'UGRD','VGRD','TMP','SPFH','PRES_meansealevel'/
  data vname/'u','v','t','gh','q','msl'/
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
  sigma = 0.0d0
  do ilev=1,nlev-1
     sigma(ilev) = sigma(ilev) + 0.5d0*dplev(ilev)/pr
     sigma(ilev+1) = sigma(ilev+1) + 0.5d0*dplev(ilev)/pr 
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
   tday = tday+1
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
   !mem=memn 
   write(cnlev,'(I1)') nlev_te
!
   call set_bound(slon,elon,slat,elat,islon,ielon,islat,ielat)
   print*,islon,ielon,islat,ielat
   nlon=ielon-islon+1
   nlat=ielat-islat+1 
   narea=nlon*nlat
! File Open 
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
   write(ns,'(I1)') smode
   write(ne,'(I1)') emode
   if (smode==emode) then
#ifdef moist
#ifndef weak
      wd='./ensvsa-TE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
     wds='./ps-ensvsa-TE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
#else
     wd='./ensvsa-wTE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
     wds='./ps-ensvsa-wTE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
#endif
#else
     wd='./ensvsa-dTE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
     wds='./ps-ensvsa-dTE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
#endif
   else
#ifdef moist
#ifndef weak
      wd='./ensvsa-TE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
     wds='./ps-ensvsa-TE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
#else
     wd='./ensvsa-wTE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
     wds='./ps-ensvsa-wTE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
#endif
#else
     wd='./ensvsa-dTE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
     wds='./ps-ensvsa-dTE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_nlev'//cnlev//'-gr'
#endif
   endif
   open(21,file=wd,status='replace',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*imax*jmax)
   open(22,file=wds,status='replace',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*imax*jmax)

  ! 配列の割付
   allocate(&
   ps(imax,jmax),ug(imax,jmax,nlev),vg(imax,jmax,nlev)&
   ,T(imax,jmax,nlev),gh(imax,jmax,nlev)&
   )
#ifdef moist
   allocate(q(imax,jmax,nlev))!,rh(imax,jmax,nlev))
#endif
   allocate(&
   zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)&
   ,z0(imax,jmax,nvar),zm(imax,jmax,nvar)&
   )
   allocate(ze(imax,jmax,nvar,mem))
   allocate(z(narea*nvar,mem))
   allocate(zT(mem,narea*nvar))
   allocate(sg(smode:emode))
   allocate(p(mem))
   allocate(w(mem,smode:emode))
   allocate(&
   TE(imax,jmax,nlev),KE(imax,jmax,nlev)&
   ,PE(imax,jmax,nlev)&
   ,v3d(imax,jmax,nlev,nv3d),v2d(imax,jmax,nv2d),buf4(imax,jmax) &
   )
#ifdef moist
   allocate(LE(imax,jmax,nlev))
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

   write(cnlev,'(I1)') nlev
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
            do id=1,4
               call fread3(rdf,vname(id),ip,zv3)
               if    (id==1)then
                  !ug=zv3(:,:,1:3)
#ifdef lev6
                  ug=zv3(:,:,2:)
#else
                  ug=zv3(:,:,:)
#endif
                  !print*,ug(1,1,1)
               elseif(id==2)then
                  !vg=zv3(:,:,1:3)
#ifdef lev6
                  vg=zv3(:,:,2:)
#else
                  vg=zv3(:,:,:)
#endif
                  !print*,vg(1,1,1)
               elseif(id==3)then
                  !T=zv3(:,:,1:3)
#ifdef lev6
                  T =zv3(:,:,2:)
#else
                  T =zv3(:,:,:)
#endif
                  !print*,T(1,1,1)
               else
#ifdef lev6
                  gh=zv3(:,:,2:)
#else
                  gh=zv3(:,:,:)
#endif
               endif
            enddo
#ifdef moist
            call fread3(rdf,vname(5),ip,zv3)
            !q=zv3(:,:,1:3)
#ifdef lev6
            q =zv3(:,:,2:)
#else
            q =zv3(:,:,:)
#endif
#endif
            call fread(rdf,vname(6),ip,zv)
            ps=dble(zv)*1.0d-2     !Pa->hPa
            !print*,ps(1,1)
           
            ze(:,:,1:nlev,imem)=ug
            ze(:,:,nlev+1:2*nlev,imem)=vg
            ze(:,:,2*nlev+1:3*nlev,imem)=T
            ze(:,:,3*nlev+1:4*nlev,imem)=gh
#ifdef moist
            ze(:,:,4*nlev+1:5*nlev,imem)=q
#endif
            ze(:,:,nvar,imem)=ps
            !ze(:,:,10,imem)=ps
           
            !! check
            write(6,'(A,I2)') 'mem ',imem
            write(6,'(A,'//cnlev//'f12.4)') 'u ',ze(1,1,1:nlev,imem)
            write(6,'(A,'//cnlev//'f12.4)') 'v ',ze(1,1,nlev+1:2*nlev,imem)
            write(6,'(A,'//cnlev//'f12.4)') 't ',ze(1,1,2*nlev+1:3*nlev,imem)
            write(6,'(A,'//cnlev//'es12.4)') 'gh',ze(1,1,3*nlev+1:4*nlev,imem)
#ifdef moist
            write(6,'(A,'//cnlev//'es12.4)') 'q ',ze(1,1,4*nlev+1:5*nlev,imem)
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
      !rdf=dir//yyyy//'/jma/'//mmddhh//'_mean.nc'   !_n
      rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_mean.nc'   !_n
      !rdf=dir//yyyy//'/jma/100900_mean.nc'
      !rdf=dir//yyyy//'/jma/anl_sellev.nc' !_a
      inquire(file=rdf, exist=ex)
      if(ex)then
         do id=1,4
            !print*,rdf
            !print*,vname(id)
            call fread3(rdf,vname(id),ip,zv3)
            !call fread3a(rdf,vnamea(id),ip,zv3,90.0d0,180.0d0,0.0d0,80.0d0)
            !print*,maxval(zv),minval(zv)
            if    (id==1)then
               !ug=zv3(:,:,1:3)
#ifdef lev6
               ug=zv3(:,:,2:)
#else
               ug=zv3(:,:,:)
#endif
               !print*,ug(1,1,1)
            elseif(id==2)then
               !vg=zv3(:,:,1:3)
#ifdef lev6
               vg=zv3(:,:,2:)
#else
               vg=zv3(:,:,:)
#endif
               !print*,vg(1,1,1)
            elseif(id==3)then
               !T=zv3(:,:,1:3)
#ifdef lev6
               T =zv3(:,:,2:)
#else
               T =zv3(:,:,:)
#endif
               !print*,T(1,1,1)
            else
               !gh=zv3(:,:,1:3)
#ifdef lev6
               gh=zv3(:,:,2:)
#else
               gh=zv3(:,:,:)
#endif
               !print*,gh(1,1,1)
            endif
         enddo
#ifdef moist
         call fread3(rdf,vname(5),ip,zv3)
         !q=zv3(:,:,1:3)
#ifdef lev6
         q =zv3(:,:,2:)
#else
         q =zv3(:,:,:)
#endif
#endif
         !rdf=dir//yyyy//'/jma/anl.nc' !_a
         call fread(rdf,vname(6),ip,zv)
         !call freada(rdf,vnamea(5),ip,zv,90.0d0,180.0d0,0.0d0,80.0d0)
         ps=dble(zv)*1.0d-2        !Pa->hPa
         !print*,ps(1,1)
         !print*,maxval(ps),minval(ps)
        
         z0(:,:,1:nlev)=ug
         z0(:,:,nlev+1:2*nlev)=vg
         z0(:,:,2*nlev+1:3*nlev)=T
         z0(:,:,3*nlev+1:4*nlev)=gh
#ifdef moist
         z0(:,:,4*nlev+1:5*nlev)=q
#endif
         z0(:,:,nvar)=ps
         !z0(:,:,10)=ps
      endif
      !! check
      write(6,'(A,'//cnlev//'f12.4)') 'u ',z0(1,1,1:nlev)
      write(6,'(A,'//cnlev//'f12.4)') 'v ',z0(1,1,nlev+1:2*nlev)
      write(6,'(A,'//cnlev//'f12.4)') 't ',z0(1,1,2*nlev+1:3*nlev)
      write(6,'(A,'//cnlev//'es12.4)') 'gh',z0(1,1,3*nlev+1:4*nlev)
#ifdef moist
      write(6,'(A,'//cnlev//'es12.4)') 'q ',z0(1,1,4*nlev+1:5*nlev)
#endif
      write(6,'(A,f12.4)') 'ps',z0(1,1,nvar)
     
     
      !1.calcurate perturbation
      do imem=1,mem
         ze(:,:,:,imem)=ze(:,:,:,imem)-z0(:,:,:)
      enddo

      !2.conserve weighted perturbation
      ssg=0.0
      do imode=smode,emode
         ssg=ssg+sg(imode)
      enddo
     
      !zm=0.0
      v3d = 0.0d0
      v2d = 0.0d0
      !prognostic variables
#ifdef moist
      do n = 1,5
#else
      do n = 1,4
#endif
         do k = 1,nlev
            ivar = (n-1)*nlev + k
            do j = 1,jmax
               do i = 1,imax      
                  do imode=smode,emode
                     do imem=1,mem
                        v3d(i,j,k,n)=v3d(i,j,k,n)+ze(i,j,ivar,imem)*w(imem,imode)*sg(imode)/ssg
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      do j=1,jmax
         do i=1,imax
            do imode = smode,emode
               do imem=1,mem
                  v2d(i,j,1) = v2d(i,j,1)+ze(i,j,nvar,imem)*w(imem,imode)*sg(imode)/ssg
               enddo
            enddo
         enddo
      enddo
      ! ps
      ps(:,:) = ps(:,:) + v2d(:,:,1)
      buf4 = real(ps, kind=4)
      write(22,rec=fday+1) buf4

      !do ilev=1,nlev            !300,500,700,850,925,1000hPa
      !   do ivar=1,4         !ug,vg,T,q
      !      ze(:,:,nlev*(ivar-1)+ilev,:)=ze(:,:,nlev*(ivar-1)+ilev,:)*sigma(ilev)
      !   enddo
      !enddo
     !3.Multiply by coefficient
     !T
      ze(:,:,2*nlev+1:3*nlev,:)=ze(:,:,2*nlev+1:3*nlev,:)*sqrt(cp/Tr)
     !ps
      ze(:,:,nvar,:)=ze(:,:,nvar,:)*sqrt(R*Tr)/pr
     !ze(:,:,10,:)=ze(:,:,10,:)*sqrt(R*Tr)/pr
#ifdef moist
     !q
      ze(:,:,4*nlev+1:5*nlev,:)=ze(:,:,4*nlev+1:5*nlev,:)*Lh/sqrt(cp*Tr)*sqrt(epsilon)
#endif
      !do imem=1,mem
      !   print*,imem
      !   print*,ze(1,1,:,imem)
      !enddo

     !4.calcurate energy
     !TE:1~4*nlev,nvar, KE=1~2*nlev, PE=2*nlev+1~3*nlev,nvar, LE=3*nvar+1~4*nlev
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
     
      ! TE(except for Ps)
      do n = 1,3
         do k = 1,nlev
            ivar = (n-1)*nlev + k
            do j = 1,jmax
               do i = 1,imax      
                  TE(i,j,k) = TE(i,j,k) + zm(i,j,ivar)**2
               enddo
            enddo
         enddo
      enddo
#ifdef moist
      do k = 1,nlev
         ivar = 4*nlev + k
         do j = 1,jmax
            do i = 1,imax      
               TE(i,j,k) = TE(i,j,k) + zm(i,j,ivar)**2
            enddo
         enddo
      enddo
#endif
      TE = TE * 0.5
      ! KE
      do n = 1,2 !u,v
         do k = 1,nlev
            ivar = (n-1)*nlev + k
            do j = 1,jmax
               do i = 1,imax      
                  KE(i,j,k) = KE(i,j,k) + zm(i,j,ivar)**2
               enddo
            enddo
         enddo
      enddo
      KE = KE * 0.5
      ! PE(T)
      do k = 1,nlev
         ivar = 2*nlev + k
         do j = 1,jmax
            do i = 1,imax      
               PE(i,j,k) = PE(i,j,k) + zm(i,j,ivar)**2
            enddo
         enddo
      enddo
      PE = PE * 0.5
      !PE(ps)
      v2d(:,:,2) = zm(:,:,nvar)**2*0.5
#ifdef moist
      !LE
      do k = 1,nlev
         ivar = 4*nlev + k
         do j = 1,jmax
            do i = 1,imax      
               LE(i,j,k) = LE(i,j,k) + zm(i,j,ivar)**2
            enddo
         enddo
      enddo
      LE = LE * 0.5
#endif
!      do ivar=1,nvar
!         TE=TE+zm(:,:,ivar)**2/2
!      enddo
!      do ivar=1,2*nlev
!         KE=KE+zm(:,:,ivar)**2/2
!      enddo
!      do ivar=2*nlev+1,3*nlev
!         PE=PE+zm(:,:,ivar)**2/2
!      enddo
!      PE=PE+zm(:,:,nvar)**2/2
!      do ivar=3*nlev+1,4*nlev
!         LE=LE+zm(:,:,ivar)**2/2
!      enddo
!         
!      !print*,"max",maxval(TE),"min",minval(TE)

!      v2d(:,:,2) = TE
!      v2d(:,:,3) = KE
!      v2d(:,:,4) = PE
!      v2d(:,:,5) = LE
#ifdef moist
      v3d(:,:,:,6) = TE
      v3d(:,:,:,7) = KE
      v3d(:,:,:,8) = PE
      v3d(:,:,:,9) = LE
#else
      v3d(:,:,:,5) = TE
      v3d(:,:,:,6) = KE
      v3d(:,:,:,7) = PE
#endif
      it=fday+1
      do n=1,nv3d
         do k=1,nlev
            buf4 = v3d(:,:,k,n)
            write(21,rec=irec) buf4
            irec = irec+1
         enddo
      enddo
      do n=1,nv2d
         buf4 = v2d(:,:,n)
         write(21,rec=irec) buf4
         irec = irec+1
      enddo
      idate=edate
      !call calc_date(idate,2,edate)
      call calc_date(idate,1,edate)
      print*, "edate=",edate
   enddo
   close(21)
   close(22)
  
   deallocate(ze,z,zT,sg,p,w)
  
   stop  
    
end program grads_ensvsa
