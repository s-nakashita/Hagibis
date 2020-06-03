program plot_ensvsa_TE

  use read_netcdf
  
  implicit none
    
  integer,parameter :: dslon=35, delon=45, dslat=67, delat=75 
  integer,parameter :: nlon=delon-dslon+1, nlat=delat-dslat+1 
  integer,parameter :: narea=nlon*nlat 
  integer,parameter :: nvar=3*3+1 
  integer,parameter :: memo=50, memn=26 
  real,parameter :: dtheta=0.5, pi=atan(1.0)*4.0 
  real,parameter :: cp=1005.7, R=287.04, Lh=2.5104*10**6 
  real,parameter :: Tr=270.0, pr=1000.0 
      
  integer :: i,j,ieof,info,l,imem,id,iarea,ilon,ilat,it
  integer :: ilt,ilu,ilv,ilq,ilev,ivar,ip,fday,imode
  integer :: lat,nd,mem,ios
  logical :: ex
  
  real ::  tmp, score, crate
  real ::  ps(imax,jmax)
  real ::  zv(0:imax-1,jmax)
  real ::  z0(nlon,nlat),zm(nlon,nlat),zf(nlon,nlat)
  real,allocatable ::  ze(:,:,:)
  real ::  sigma(3),area
  real,allocatable ::  z(:,:),zT(:,:)
  real,allocatable :: sg(:),p(:),w(:,:),slp(:,:)
  real :: slpm(ntime),slpf(ntime),day(ntime)
      
  character rdf*100,rdw*100,wd*100,wdm*100,wdf*100
  character dir*40,nmem*2,fd*1,date*10,yyyy*4,mmddhh*6,yyyymmddhh*10,orig*4
  character(len=17) :: vname(5)
  data vname/'TMP','UGRD','VGRD','SPFH','PRES_meansealevel'/

     !|----/----/----/----/----/----/----/----/----/----| 
  dir='/Users/nakashita/Documents/netcdf/tigge/'

  sigma(1)=8.0/7.0*300.0/pr
  sigma(2)=6.0/7.0*300.0/pr
  sigma(3)=200.0/pr
  print*,sigma
      
  !  データの設定 
  yyyymmddhh="2019100912"
  yyyy=yyyymmddhh(1:4)
  mmddhh=yyyymmddhh(5:10)
  print*,yyyy,mmddhh
  
  mem=memn 
  rdw='./weight-slp-jma-'//yyyymmddhh//'.grd'
  open(10,file=rdw,status='old',access='direct',&
       &        convert='big_endian',&
       &        form='unformatted', recl=4*mem)

  ! 配列の割付
  allocate(ze(nlon,nlat,mem))
  allocate(z(narea,mem))
  allocate(zT(mem,narea))
  allocate(sg(10))
  allocate(p(mem))
  allocate(w(mem,10))
  allocate(slp(mem,ntime))
  
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
  
  do ip=1,ntime
     day(ip)=6*(ip-1)
     do imem=1,mem
        write(nmem,'(I2.2)') imem
        !print*,nmem
        rdf=dir//yyyy//'/jma/'//mmddhh//'_'//nmem//'.nc'
        !print*,rdf
        inquire(file=rdf, exist=ex)
        if(ex)then
           call fread(rdf,vname(5),ip,zv)
           ps=zv/100     !Pa->hPa
            
           do i=1,nlon
              do j=1,nlat
                 ilon=dslon+i-1
                 ilat=dslat+j-1
                 ze(i,j,imem)=ps(ilon,ilat)
              enddo
           enddo
           !print*,ze(1,1,:,imem)
        endif
     enddo

     rdf=dir//yyyy//'/jma/100900_mean.nc'
     inquire(file=rdf, exist=ex)
     if(ex)then
        call fread(rdf,vname(5),ip,zv)
        ps=zv/100        !Pa->hPa
        !print*,maxval(ps),minval(ps)
        
        do i=1,nlon
           do j=1,nlat
              ilon=dslon+i-1
              ilat=dslat+j-1
              z0(i,j)=ps(ilon,ilat)
           enddo
        enddo
     endif
     !print*,z0(1,1,:)
         
       
     !1.calcurate perturbation
     do imem=1,mem
        ze(:,:,imem)=ze(:,:,imem)-z0(:,:)
     enddo
     !2.Multiply by cos(lat) and layer thickness factor
     area=0.0
     do i=1,nlon
        do j=1,nlat
           lat=dslat+(j-1)-1
           area=area+cos(lat*dtheta*pi/180.0)
           ze(i,j,:)=ze(i,j,:)*cos(lat*dtheta*pi/180.0)
        enddo
     enddo
      
     do imem=1,mem
        print*,imem
        print*,ze(1,1,imem)
     enddo

     !4.calcurate TE
     slp(:,ip)=0.0
     do imem=1,mem
        do i=1,nlon
            do j=1,nlat
               slp(imem,ip)=slp(imem,ip)+ze(i,j,imem)**2/area
            enddo
        enddo
     enddo
     slp = sqrt(slp)
     print*,"max",maxval(slp(:,ip)),"min",minval(slp(:,ip))
     
     !mean
     slpm(ip)=0.0
     do imem=1,mem
        slpm(ip)=slpm(ip)+slp(imem,ip)/mem
     enddo
     print*,"mean",slpm(ip)
     
     !first mode
     zm=0.0
     do imem=1,mem
        zm(:,:)=zm(:,:)+ze(:,:,imem)*w(imem,1)
     enddo
     slpf(ip)=0.0
     do i=1,nlon
        do j=1,nlat
           slpf(ip)=slpf(ip)+zm(i,j)**2/area
        enddo
     enddo
     slpf = sqrt(slpf)     
     print*,"1 mode",slpf(ip)
  end do
  
  wdm='./slp-mean-jma-'//yyyymmddhh//'.txt'
  open(21,file=wdm,status='new')
  wdf='./slp-1mode-jma-'//yyyymmddhh//'.txt'
  open(22,file=wdf,status='new')
  do ip=1,ntime
     write(21,'(f5.1,f10.4)') day(ip),slpm(ip)
     write(22,'(f5.1,f10.4)') day(ip),slpf(ip)
  enddo
  close(21)
  close(22)
  
  do imem=1,mem
     write(nmem,'(I2.2)') imem
     wd='./slp-'//nmem//'-jma-'//yyyymmddhh//'.txt'
     open(23,file=wd,status='new')
     do ip=1,ntime
        write(23,'(f5.1,f10.4)') day(ip),slp(imem,ip)
     enddo
     close(23)
  enddo
      
  deallocate(ze,z,zT,sg,p,w,slp)
  
  stop  
    
contains

end program plot_ensvsa_TE
