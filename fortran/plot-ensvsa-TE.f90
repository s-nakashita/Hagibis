program plot_ensvsa_TE

  use read_netcdf
  
  implicit none
    
  integer,parameter :: dslon=35, delon=45, dslat=67, delat=75 
  integer,parameter :: nlon=delon-dslon+1, nlat=delat-dslat+1 
  integer,parameter :: narea=nlon*nlat 
  integer,parameter :: nvar=3*kmax+1 
  integer,parameter :: memo=50, memn=26 
  real,parameter :: dtheta=0.5, pi=atan(1.0)*4.0 
  real,parameter :: cp=1005.7, R=287.04, Lh=2.5104*10**6 
  real,parameter :: Tr=270.0, pr=1000.0 
      
  integer :: i,j,ieof,info,l,imem,id,iarea,ilon,ilat,it
  integer :: ilt,ilu,ilv,ilq,ilev,ivar,ip,fday,imode
  integer :: lat,nd,mem,ios
  logical :: ex
  
  real ::  tmp, score, crate
  real ::  ps(imax,jmax),ug(imax,jmax,3),vg(imax,jmax,3)
  real ::  T(imax,jmax,3),q(imax,jmax,3)
  real ::  zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)
  real ::  z0(nlon,nlat,nvar),zm(nlon,nlat,nvar),zf(nlon,nlat,nvar)
  real,allocatable ::  ze(:,:,:,:)
  real ::  sigma(3),area
  real,allocatable ::  z(:,:),zT(:,:)
  real,allocatable :: sg(:),p(:),w(:,:),TE(:,:)
  real :: TEm(ntime),TEf(ntime),day(ntime)
      
  character rdf*100,rdw*100,wd*100,wdm*100,wdf*100
  character dir*48,nmem*2,fd*1,date*10,yyyy*4,mmddhh*6,yyyymmddhh*10,orig*4
  character(len=17) :: vname(5)
  data vname/'TMP','UGRD','VGRD','SPFH','PRES_meansealevel'/

     !|----/----/----/----/----/----/----/----/----/----| 
  dir='/Users/nakashita/Documents/hagibis/netcdf/tigge/'

  sigma(1)=200.0/pr
  sigma(2)=6.0/7.0*300.0/pr
  sigma(3)=8.0/7.0*300.0/pr
  print*,sigma
      
  !  データの設定 
  yyyymmddhh="2019101000"
  yyyy=yyyymmddhh(1:4)
  mmddhh=yyyymmddhh(5:10)
  print*,yyyy,mmddhh
  
  mem=memn 
  rdw='./weight-dryTE-jma-'//yyyymmddhh//'.grd'
  open(10,file=rdw,status='old',access='direct',&
       &        convert='big_endian',&
       &        form='unformatted', recl=4*mem)

  ! 配列の割付
  allocate(ze(nlon,nlat,nvar,mem))
  allocate(z(narea*nvar,mem))
  allocate(zT(mem,narea*nvar))
  allocate(sg(10))
  allocate(p(mem))
  allocate(w(mem,10))
  allocate(TE(mem,ntime))
  
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
        ilt=0
        ilu=0
        ilv=0
        ilq=0
        rdf=dir//yyyy//'/jma/'//mmddhh//'_'//nmem//'.nc'
        !print*,rdf
        inquire(file=rdf, exist=ex)
        if(ex)then
           do id=1,4
              call fread3(rdf,vname(id),ip,zv3)
              if(mod(id,4)==1)then
                 T(:,:,:)=zv3
              elseif(mod(id,4)==2)then
                 ug(:,:,:)=zv3
              elseif(mod(id,4)==3)then
                 vg(:,:,:)=zv3
              else
                 q(:,:,:)=zv3
              endif
           enddo
           
           call fread(rdf,vname(5),ip,zv)
           ps=zv/100     !Pa->hPa
            
           do i=1,nlon
              do j=1,nlat
                 ilon=dslon+i-1
                 ilat=dslat+j-1
                 ze(i,j,1:3,imem)=ug(ilon,ilat,:)
                 ze(i,j,4:6,imem)=vg(ilon,ilat,:)
                 ze(i,j,7:9,imem)=T(ilon,ilat,:)
                 !ze(i,j,10:12,imem)=q(ilon,ilat,:)
                 !ze(i,j,13,imem)=ps(ilon,ilat)
                 ze(i,j,10,imem)=ps(ilon,ilat)
              enddo
           enddo
           !print*,ze(1,1,:,imem)
        endif
     enddo

     ilt=0
     ilu=0
     ilv=0
     ilq=0
     rdf=dir//yyyy//'/jma/'//mmddhh//'_mean.nc'
     inquire(file=rdf, exist=ex)
     if(ex)then
        do id=1,4
           !print*,rdf
           !print*,vname(id)
           call fread3(rdf,vname(id),ip,zv3)
           !print*,maxval(zv),minval(zv)
           if(mod(id,4)==1)then
              T(:,:,:)=zv3
           elseif(mod(id,4)==2)then
              ug(:,:,:)=zv3
           elseif(mod(id,4)==3)then
              vg(:,:,:)=zv3
           else
              q(:,:,:)=zv3
           endif
        enddo
        
        call fread(rdf,vname(5),ip,zv)
        ps=zv/100        !Pa->hPa
        !print*,maxval(ps),minval(ps)
        
        do i=1,nlon
           do j=1,nlat
              ilon=dslon+i-1
              ilat=dslat+j-1
              z0(i,j,1:3)=ug(ilon,ilat,:)
              z0(i,j,4:6)=vg(ilon,ilat,:)
              z0(i,j,7:9)=T(ilon,ilat,:)
              !z0(i,j,10:12)=q(ilon,ilat,:)
              !z0(i,j,13)=ps(ilon,ilat)
              z0(i,j,10)=ps(ilon,ilat)
           enddo
        enddo
     endif
     !print*,z0(1,1,:)
         
       
     !1.calcurate perturbation
     do imem=1,mem
        ze(:,:,:,imem)=ze(:,:,:,imem)-z0(:,:,:)
     enddo
     !2.Multiply by cos(lat) and layer thickness factor
     area=0.0
     do i=1,nlon
        do j=1,nlat
           lat=dslat+(j-1)-1
           area=area+cos(lat*dtheta*pi/180.0)
           ze(i,j,:,:)=ze(i,j,:,:)*cos(lat*dtheta*pi/180.0)
        enddo
     enddo
      
     do ilev=1,3            !250,500,850hPa
        do ivar=1,3!4         !ug,vg,T,q
           ze(:,:,3*(ivar-1)+ilev,:)=ze(:,:,3*(ivar-1)+ilev,:)*sigma(ilev)
        enddo
     enddo
     !3.Multiply by coefficient
     !T
     ze(:,:,7:9,:)=ze(:,:,7:9,:)*sqrt(cp/Tr)
     !ps
     !ze(:,:,13,:)=ze(:,:,13,:)*sqrt(R*Tr)/pr
     ze(:,:,10,:)=ze(:,:,10,:)*sqrt(R*Tr)/pr
     !q
     !ze(:,:,10:12,:)=ze(:,:,10:12,:)*Lh/sqrt(cp*Tr)

     do imem=1,mem
        print*,imem
        print*,ze(1,1,:,imem)
     enddo

     !4.calcurate TE
     TE(:,ip)=0.0
     do imem=1,mem
        do ivar=1,nvar
           do i=1,nlon
              do j=1,nlat
                 TE(imem,ip)=TE(imem,ip)+ze(i,j,ivar,imem)**2/2/area
              enddo
           enddo
        enddo
     enddo
     print*,"max",maxval(TE(:,ip)),"min",minval(TE(:,ip))
     
     !mean
     TEm(ip)=0.0
     do imem=1,mem
        TEm(ip)=TEm(ip)+TE(imem,ip)/mem
     enddo
     print*,"mean",TEm(ip)
     
     !first mode
     zm=0.0
     do imem=1,mem
        zm(:,:,:)=zm(:,:,:)+ze(:,:,:,imem)*w(imem,1)
     enddo
     TEf(ip)=0.0
     do ivar=1,nvar
        do i=1,nlon
           do j=1,nlat
              TEf(ip)=TEf(ip)+zm(i,j,ivar)**2/2/area
           enddo
        enddo
     enddo
     
     print*,"1 mode",TEf(ip)
  end do
  
  wdm='./TE-mean-jma-'//yyyymmddhh//'.txt'
  open(21,file=wdm,status='new')
  wdf='./TE-1mode-jma-'//yyyymmddhh//'.txt'
  open(22,file=wdf,status='new')
  do ip=1,ntime
     write(21,'(f5.1,f9.4)') day(ip),TEm(ip)
     write(22,'(f5.1,f9.4)') day(ip),TEf(ip)
  enddo
  close(21)
  close(22)
  
  do imem=1,mem
     write(nmem,'(I2.2)') imem
     wd='./TE-'//nmem//'-jma-'//yyyymmddhh//'.txt'
     open(23,file=wd,status='new')
     do ip=1,ntime
        write(23,'(f5.1,f9.4)') day(ip),TE(imem,ip)
     enddo
     close(23)
  enddo
      
  deallocate(ze,z,zT,sg,p,w,TE)
  
  stop  
    
contains

end program plot_ensvsa_TE
