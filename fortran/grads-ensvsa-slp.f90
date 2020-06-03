program grads_ensvsa_TE

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
  real ::  z0(imax,jmax),zm(imax,jmax)
  real,allocatable ::  ze(:,:,:)
  real ::  sigma(3),area,ssg
  real,allocatable ::  z(:,:),zT(:,:)
  real,allocatable :: sg(:),p(:),w(:,:)
  real :: slp(imax,jmax)
      
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
  wd='./ensvsa-slp-m1-jma-'//yyyymmddhh//'-gr'
  open(21,file=wd,status='new',access='direct',&
       &        convert='big_endian',&
       &        form='unformatted', recl=4*imax*jmax*1)
  
  mem=memn 
  rdw='./weight-slp-jma-'//yyyymmddhh//'.grd'
  open(10,file=rdw,status='old',access='direct',&
       &        convert='big_endian',&
       &        form='unformatted', recl=4*mem)

  ! 配列の割付
  allocate(ze(imax,jmax,mem))
  allocate(z(narea,mem))
  allocate(zT(mem,narea))
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

  do fday=0,7 !every 12 hours
     ip=1+2*fday
     do imem=1,mem
        write(nmem,'(I2.2)') imem
        !print*,nmem
        rdf=dir//yyyy//'/jma/'//mmddhh//'_'//nmem//'.nc'
        !print*,rdf
        inquire(file=rdf, exist=ex)
        if(ex)then
           call fread(rdf,vname(5),ip,zv)
           ps=zv/100     !Pa->hPa
           ze(:,:,imem)=ps
           !print*,ze(1,1,:,imem)
        endif
     enddo
     
     rdf=dir//yyyy//'/jma/100900_mean.nc'
     inquire(file=rdf, exist=ex)
     if(ex)then
        call fread(rdf,vname(5),ip,zv)
        ps=zv/100        !Pa->hPa
        z0(:,:)=ps
     endif
     !print*,z0(1,1,:)
     
     !1.calcurate perturbation
     do imem=1,mem
        ze(:,:,imem)=ze(:,:,imem)-z0(:,:)
     enddo
     
     do imem=1,mem
        print*,imem
        print*,ze(1,1,imem)
     enddo

     !4.calcurate TE
     ssg=0.0
     do imode=1,10
        ssg=ssg+sg(imode)
     enddo
     
     slp=0.0
     zm=0.0
     do imode=1,10
        do imem=1,mem
           zm=ze(:,:,imem)*w(imem,imode)*sg(imode)/ssg
        enddo
     enddo
     
     slp=sqrt(zm**2)
         
     print*,"max",maxval(slp),"min",minval(slp)

     it=fday+1
     write(21,rec=it) slp
  enddo
  close(21)
  
  deallocate(ze,z,zT,sg,p,w)
  
  stop  
    
end program grads_ensvsa_TE
