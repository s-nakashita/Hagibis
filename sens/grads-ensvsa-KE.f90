program grads_ensvsa_KE

  use read_netcdf
  
  implicit none
 
  integer,parameter :: dslon=95, delon=105, dslat=67, delat=75 
  integer,parameter :: nlon=delon-dslon+1, nlat=delat-dslat+1 
  integer,parameter :: narea=nlon*nlat 
  integer,parameter :: nv3d=2,nlev=3
  integer,parameter :: nvar=nv3d*nlev
  integer,parameter :: memo=50, memn=26 
  real,parameter :: dtheta=0.5, pi=atan(1.0)*4.0 
  real,parameter :: cp=1005.7, R=287.04, Lh=2.5104*10**6 
  real,parameter :: Tr=270.0, pr=1000.0 
      
  integer :: i,j,k,n
  integer :: imem,id,it,irec
  integer :: ilt,ilu,ilv,ilq,ilev,ivar,ip,fday,imode
  integer :: mem
  integer :: idate,edate
  logical :: ex
  
  real ::  ug(imax,jmax,nlev),vg(imax,jmax,nlev)
  real ::  zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)
  real ::  z0(imax,jmax,nvar),zm(imax,jmax,nvar)
  real,allocatable ::  ze(:,:,:,:)
  real ::  sigma(3),ssg
  real :: plev(3)
  data plev/850.0,500.0,300.0/
  real,allocatable ::  z(:,:),zT(:,:)
  real,allocatable :: sg(:),p(:),w(:,:)
  real :: KE(imax,jmax)
  real :: v3d(imax,jmax,nlev,nv3d) !ug,vg
  real :: v2d(imax,jmax) !TE
  real :: buf4(imax,jmax)
      
  character rdf*100,rdw*100,wd*100
  character dir*30,dira*33,nmem*2,yyyy*4,mm*2,mmddhh*6,yyyymmddhh*10
  character(len=17) :: vname(2)
  !character(len=4) :: vnamea(5)
  data vname/'UGRD','VGRD'/
  !data vnamea/'air','uwnd','vwnd','shum','slp'/
     !|----/----/----/----/----/----/----/----/----/----| 
  dir='/Users/nakashita/netcdf/tigge/'
 !dira='/Users/nakashita/netcdf/nc-reanl/'

   sigma(1)=8.0/7.0*300.0/pr
   sigma(2)=6.0/7.0*300.0/pr
   sigma(3)=200.0/pr
   print*,sigma
      
  !  データの設定 
   yyyymmddhh="2019100912"
   yyyy=yyyymmddhh(1:4)
   mm=yyyymmddhh(5:6)
   mmddhh=yyyymmddhh(5:10)
   print*,yyyy,mmddhh
   wd='./ensvsa-KE-m1-jma-'//yyyymmddhh//'_a-gr'
   open(21,file=wd,status='new',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*imax*jmax)
  
   mem=memn 
   rdw='./weight-KE-jma-'//yyyymmddhh//'_a.grd'
   open(10,file=rdw,status='old',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*mem)

  ! 配列の割付
   allocate(ze(imax,jmax,nvar,mem))
   allocate(z(narea*nvar,mem))
   allocate(zT(mem,narea*nvar))
   allocate(sg(10))
   allocate(p(mem))
   allocate(w(mem,10))

  !singular value
   read(10,rec=1) p
   do imode=1,10
      sg(imode)=p(imode)
   enddo
   print*, sg
  
   do imode=1,10
      it=imode+1
      read(10,rec=it) p
      w(:,imode)=p(:)
      print*,imode,w(:,imode)
   enddo
  
   close(10)

   irec=1
   edate=2019100912
   print*,"edate=",edate
   do fday=0,6 !every 12 hours
      idate=2019100912
      call calc_steps(idate,edate,6,ip)
      ip=ip+1
      print *, "ip=",ip
      do imem=1,mem
         write(nmem,'(I2.2)') imem
         print*,nmem
         ilt=0
         ilu=0
         ilv=0
         ilq=0
         rdf=dir//yyyy//'/jma/'//mmddhh//'_'//nmem//'.nc'
         !print*,rdf
         inquire(file=rdf, exist=ex)
         if(ex)then
            call fread3(rdf,vname(1),ip,zv3)
            ug=zv3(:,:,1:3)
            print*,ug(1,1,1)
            call fread3(rdf,vname(2),ip,zv3)
            vg=zv3(:,:,1:3)
            print*,vg(1,1,1)
            
            ze(:,:,1:3,imem)=ug
            ze(:,:,4:6,imem)=vg
            !print*,ze(1,1,:,imem)
         endif
      enddo
     
      idate=2019100900
      call calc_steps(idate,edate,12,ip)
      ip=ip+1
      print *, "ip=",ip
      !ip=fday+2
      !print*,ip
      ilt=0
      ilu=0
      ilv=0
      ilq=0
      !rdf=dir//yyyy//'/jma/'//mmddhh//'_mean.nc'   !_n
      !rdf=dir//yyyy//'/jma/100900_mean.nc'
      rdf=dir//yyyy//'/jma/anl_sellev.nc' !_a
      inquire(file=rdf, exist=ex)
      if(ex)then
         call fread3(rdf,vname(1),ip,zv3)
         ug=zv3(:,:,1:3)
         print*,ug(1,1,1)
         call fread3(rdf,vname(2),ip,zv3)
         vg=zv3(:,:,1:3)
         print*,vg(1,1,1)
         
         z0(:,:,1:3)=ug
         z0(:,:,4:6)=vg
      endif
      !print*,z0(1,1,:)
     
     
      !1.calcurate perturbation
      do imem=1,mem
         ze(:,:,:,imem)=ze(:,:,:,imem)-z0(:,:,:)
      enddo

      !2.conserve weighted perturbation
      ssg=0.0
      do imode=1,1
         ssg=ssg+sg(imode)
      enddo
     
      zm=0.0
      do n = 1,nv3d
         do k = 1,nlev
            ivar = (n-1)*nlev + k
            do j = 1,jmax
               do i = 1,imax      
                  do imode=1,1
                     do imem=1,mem
                        v3d(i,j,k,n)=ze(i,j,ivar,imem)*w(imem,imode)*sg(imode)/ssg
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      do ilev=1,3            !850,500,300hPa
         do ivar=1,2         !ug,vg
            ze(:,:,3*(ivar-1)+ilev,:)=ze(:,:,3*(ivar-1)+ilev,:)*sigma(ilev)
         enddo
      enddo
     
      do imem=1,mem
         print*,imem
         print*,ze(1,1,:,imem)
      enddo

     !4.calcurate KE
      KE=0.0
      zm=0.0
      do imode=1,1
         do imem=1,mem
            zm=ze(:,:,:,imem)*w(imem,imode)*sg(imode)/ssg
         enddo
      enddo
     
      do ivar=1,nvar
         KE=KE+zm(:,:,ivar)**2/2
      enddo
         
      print*,"max",maxval(KE),"min",minval(KE)

      v2d = KE

      it=fday+1
      do n=1,nv3d
         do k=1,nlev
            buf4 = v3d(:,:,k,n)
            write(21,rec=irec) buf4
            irec = irec+1
         enddo
      enddo
      buf4 = v2d
      write(21,rec=irec) buf4
      irec = irec+1
      
      idate=edate
      call calc_date(idate,2,edate)
      print*, "edate=",edate
   enddo
   close(21)
  
   deallocate(ze,z,zT,sg,p,w)
  
   stop  
    
end program grads_ensvsa_KE
