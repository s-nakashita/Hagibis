program ensvsa

  use read_netcdf
  
  implicit none
    
  ! Parameters
#ifdef lev6
  integer,parameter :: nlev=6 
#else
  integer,parameter :: nlev=3
#endif
#ifdef moist
  integer,parameter :: nvar=nlev*4+1 
#else
  integer,parameter :: nvar=nlev*3+1 
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

  double precision :: lat, area
  integer :: i,j,ieof,info,l,imem,id,iarea,ilon,ilat,it
  integer :: ilt,ilu,ilv,ilq,ilev,ivar,ip
  logical :: ex
  
  real ::  tmp, score, crate
  real,allocatable ::  zv3(:,:,:),zv(:,:)
  double precision,allocatable ::  ps(:,:)
  double precision,allocatable ::  ug(:,:,:),vg(:,:,:),T(:,:,:)
#ifdef moist
  double precision,allocatable ::  q(:,:,:)!,rh(imax,jmax,3)
#endif
  double precision,allocatable ::  z0(:,:,:)
  double precision,allocatable ::  ze(:,:,:,:)
  double precision ::  sigma(nlev)
  double precision :: plev(nlev),dplev(nlev-1)
#ifdef lev6
   data plev/300.0d0,500.0d0,700.0d0,850.0d0,925.0d0,1000.0d0/
#else
   data plev/300.0,500.0,850.0/
#endif
  double precision,allocatable ::  z(:,:),zT(:,:)
  double precision,allocatable :: u8(:,:),vt8(:,:)
  real,allocatable :: vt(:,:),v(:,:),vtv(:,:)
  integer :: lwork!=3*narea*nvar 
  double precision,allocatable :: sg8(:),work(:)
  real,allocatable :: sg(:),p(:)
  
  ! Namelist
  character(len=5) :: orig="jma"
  integer :: mem=26 
  integer :: idate=2019100912
  integer :: edate=2019101212
  integer :: smode=1
  integer :: emode=1
  namelist /sens_nml/ orig, mem, idate, edate, smode, emode

  character rdf*100,wd*100,logf*100
  character dir*32!,dira*33
  character nmem*2,cnlev*1
  character yyyy*4,mm*2,mmddhh*6,yyyymmddhh*10
  character(len=3) :: vname(5)
  !character(len=4) :: vnamea(5)
  !data vname/'UGRD','VGRD','TMP','SPFH','PRES_meansealevel'/
  data vname/'u','v','t','q','msl'/
  !data vnamea/'air','uwnd','vwnd','shum','slp'/
  
     !|----/----/----/----/----/----/----/----/----/----/----/----| 
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
   sigma(ilev) = sigma(ilev) + 0.5d0*dplev(ilev)/pr
   sigma(ilev+1) = sigma(ilev+1) + 0.5d0*dplev(ilev)/pr 
  end do
  print*,sigma
  
! Settings for variables
  open(11,file="sens.nml")
  read(11,nml=sens_nml)
  write(*,nml=sens_nml)
  close(11)
  !idate=2019100912
  write(yyyymmddhh,'(I10)') idate
  yyyy=yyyymmddhh(1:4)
  mm=yyyymmddhh(5:6)
  mmddhh=yyyymmddhh(5:10)
  print*,yyyy,mmddhh
  write(cnlev,'(I1)') nlev
!
  call set_bound(slon,elon,slat,elat,islon,ielon,islat,ielat)
  print*,islon,ielon,islat,ielat
  nlon=ielon-islon+1
  nlat=ielat-islat+1 
  narea=nlon*nlat
  lwork=3*narea*nvar 
  !mem=memn
  it=1
  !edate=2019101212
  call calc_steps(idate,edate,6,ip)
  print*,ip
  ip = ip+1
  !ip=13
#ifdef moist
#ifndef weak
  wd='./weight-TE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.grd'
  logf='./TE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.log'
#else
  wd='./weight-wTE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.grd'
  logf='./wTE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.log'
#endif
#else
  wd='./weight-dTE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.grd'
  logf='./dTE-'//trim(orig)//'-'//yyyymmddhh//'_nlev'//cnlev//'.log'
#endif
  open(21,file=wd,status='replace',access='direct',&
       &        convert='big_endian',&
       &        form='unformatted', recl=4*mem)
  open(30,file=logf)
  ! 配列の割付
  l=min(mem,narea*nvar)
  print*,l
  allocate(ps(imax,jmax),ug(imax,jmax,nlev),vg(imax,jmax,nlev),T(imax,jmax,nlev))
#ifdef moist
  allocate(q(imax,jmax,nlev))
#endif
  allocate(zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax))
  allocate(z0(nlon,nlat,nvar),ze(nlon,nlat,nvar,mem))
  allocate(z(narea*nvar,mem))
  allocate(zT(mem,narea*nvar))
  !allocate(a8(narea*nvar,mem))
  allocate(sg8(l))
  allocate(sg(l))
  allocate(u8(narea*nvar,narea*nvar))  
  allocate(vt8(mem,mem))
  allocate(vt(mem,mem))
  allocate(p(mem))
  allocate(v(mem,mem))
  allocate(vtv(mem,mem))
  allocate(work(lwork))
  
  ilt=0
  ilu=0
  ilv=0
  ilq=0
  !rdf=dir//yyyy//'/jma/'//mmddhh//'_mean.nc'   !_n
  rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_mean.nc'   !_n
  !rdf=dir//yyyy//'/jma/100900_mean.nc'
  !rdf=dir//yyyy//'/jma/anl_sellev.nc' !_a
  !ip = ip+2 !nomark
  !idate=2019100900
  !call calc_steps(idate,edate,12,ip)
  !print*,ip
  !ip = ip+1
  !ip = 8 !_a
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

      !rdf=dira//yyyy//'/'//mm//'/init.nc' !_a
      !rdf=dir//yyyy//'/jma/anl.nc' !_a
      call fread(rdf,vname(5),ip,zv)
      !call freada(rdf,vnamea(5),ip+2,zv,90.0d0,180.0d0,0.0d0,80.0d0)
      ps=dble(zv)*1.0d-2        !Pa->hPa
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
   do imem=1,mem
      write(nmem,'(I2.2)') imem
      !print*,nmem
      ilt=0
      ilu=0
      ilv=0
      ilq=0
      !rdf=dir//yyyy//'/jma/'//mmddhh//'_'//nmem//'.nc'
      rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_'//nmem//'.nc'
      print*,rdf
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
#ifdef lev6
         q =dble(zv3(:,:,4:))
#else
         q =dble(zv3(:,:,3:))
#endif
#endif
         call fread(rdf,vname(5),ip,zv)
         ps=dble(zv)*1.0d-2     !Pa->hPa

  !1.calcurate perturbation
         do i=1,nlon
            do j=1,nlat
               ilon=islon+i-1
               ilat=islat+j-1
               ze(i,j,1:nlev,imem)=ug(ilon,ilat,:)-z0(i,j,1:nlev)
               ze(i,j,nlev+1:2*nlev,imem)=vg(ilon,ilat,:)-z0(i,j,nlev+1:2*nlev)
               ze(i,j,2*nlev+1:3*nlev,imem)=T(ilon,ilat,:)-z0(i,j,2*nlev+1:3*nlev)
#ifdef moist
               ze(i,j,3*nlev+1:4*nlev,imem)=q(ilon,ilat,:)-z0(i,j,3*nlev+1:4*nlev)
#endif
               ze(i,j,nvar,imem)=ps(ilon,ilat)-z0(i,j,nvar)
               !ze(i,j,10,imem)=ps(ilon,ilat)
            enddo
         enddo
      endif
   enddo 
    
   deallocate(zv3,zv,z0,ug,vg,T,ps)
#ifdef moist
   deallocate(q)
#endif
  !do imem=1,mem
  !   ze(:,:,:,imem)=ze(:,:,:,imem)-z0(:,:,:)
  !enddo
  !2.Multiply by area and layer thickness factor
   area=0.0d0
   do j=1,nlat
      lat=(islat+j-2)*dlat + lat0
      area=area+cos(lat*pi/180.0)
      ze(:,j,:,:)=ze(:,j,:,:)*sqrt(cos(lat*pi/180.0))
   enddo
   area=area*nlon
   ze=ze/sqrt(area)
   !print*, ze(1,1,:,:)   
   do ilev=1,nlev         !300,500,700,850,925,1000hPa
      do ivar=1,3         !ug,vg,T
         ze(:,:,3*(ivar-1)+ilev,:)=ze(:,:,3*(ivar-1)+ilev,:)*sqrt(sigma(ilev))
      enddo
#ifdef moist
      ivar=4 !q
      ze(:,:,3*(ivar-1)+ilev,:)=ze(:,:,3*(ivar-1)+ilev,:)*sqrt(sigma(ilev))
#endif
   enddo
   !print*, ze(1,1,:,:)   
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
  !! check
   do imem=1,mem
      write(6,'(A,I2)') 'mem ',imem
      write(6,'(A,'//cnlev//'f12.4)') 'u ',ze(1,1,1:nlev,imem)
      write(6,'(A,'//cnlev//'f12.4)') 'v ',ze(1,1,nlev+1:2*nlev,imem)
      write(6,'(A,'//cnlev//'f12.4)') 't ',ze(1,1,2*nlev+1:3*nlev,imem)
#ifdef moist
      write(6,'(A,'//cnlev//'f12.4)') 'q ',ze(1,1,3*nlev+1:4*nlev,imem)
#endif
      write(6,'(A,f12.4)') 'ps',ze(1,1,nvar,imem)
   enddo

  !4.make data array Z
  
   do imem=1,mem
      iarea=0
      do i=1,nlon
         do j=1,nlat
            do ivar=1,nvar
               iarea=iarea+1
               z(iarea,imem)=ze(i,j,ivar,imem)
            enddo
         enddo
      enddo
   enddo
  !print*,z
  ! 転置行列を作る
   do j=1,mem
      do i=1,narea*nvar
         zT(j,i) = z(i,j)
      end do
   end do
  
   !a8=real(z,kind=8)
   call dgesvd('o','a',narea*nvar,mem,z,narea*nvar,sg8,u8,narea*nvar,vt8,mem,work,lwork,info)
   !特異値分解(eof.f90)と異なり、固有値の小さい順にデータが格納されていることに注意
   write(30,*) 'info=',info
   !print *,'eigen values=',s(:) !固有値
   !do i=l,1,-1
   !   sg(l-i+1)=sqrt(mem*s(i))
   !enddo
   !print*,"0"
   !call calc_svd(z,int(narea*nvar),mem,sg,vt)

   vt=real(vt8,kind=4)
   sg=real(sg8,kind=4)
  
   write(30,*) 'singular values (double)=', sg8(1:10) !特異値(sigma=sqrt(m*s))
   write(30,*) 'singular values=', sg(1:10) !特異値(sigma=sqrt(m*s))
   write(21,rec=it) sg
   it=it+1

   tmp=0.0
   do ieof=1,l 
      tmp=tmp+sg(ieof)*sg(ieof) 
   end do
      
   !s=sg
   crate=0.0
   !do ieof=l,l-9,-1          !固有値の大きい順に直して表示
   do ieof=1,10
      write(30,*) '-- PC',ieof
      write(30,*) 'proportion=',sg(ieof)*sg(ieof)/tmp
      crate=crate+sg(ieof)*sg(ieof)/tmp*100
      write(30,*) 'contribution rate(%)=',crate
      write(30,*) 'vector=',vt(ieof,:)
      write(30,*) 'points:'
      do i=1,10
         score = 0.0
         do j=1,mem
            score = score + real(z(i,j)*vt(ieof,j),kind=4)
         enddo
         write(30,*) i,score
        !     print *,i,dot_product(z(i,:),vt(ieof,:)) !<-特異値分解のu(i,ieof)*s(ieof)に相当
      end do
      p(:)=vt(ieof,:)
      write(21,rec=it) p
      it=it+1
   end do
  
   do i=1,mem
      do j=1,mem
         v(j,i)=vt(i,j)
      enddo
   enddo
  
   vtv=matmul(vt,v)
   print*, vtv(3,:)
      
  
   close(21)
   deallocate(ze,z,zT,sg8,u8,vt8,vt,sg,p,v,vtv,work)
    
end program ensvsa
