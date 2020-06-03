program ensvsa_TE

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
  integer :: ilt,ilu,ilv,ilq,ilev,ivar,ip,fday
  integer :: lat,day,nd,mem,ios
  logical :: ex
  
  real ::  tmp, score, crate
  real ::  ps(imax,jmax),ug(imax,jmax,3),vg(imax,jmax,3)
  real ::  T(imax,jmax,3),q(imax,jmax,3)
  real ::  zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)
  real ::  z0(nlon,nlat,nvar)
  real,allocatable ::  ze(:,:,:,:)
  real ::  sigma(3)
  real,allocatable ::  z(:,:),zT(:,:)
  double precision,allocatable :: a8(:,:),u8(:,:),vt8(:,:)
  real,allocatable :: vt(:,:),v(:,:),vtv(:,:)
  integer,parameter :: lwork=3*narea*nvar 
  double precision,allocatable :: sg8(:),work(:)
  real,allocatable :: sg(:),p(:)
  
  character rdf*100,wd*100,rdt*100
  character dir*40,nmem*2,fd*1,yyyy*4,mmddhh*6,yyyymmddhh*10,orig*4
  character(len=17) :: vname(5)
  data vname/'TMP','UGRD','VGRD','SPFH','PRES_meansealevel'/
  
     !|----/----/----/----/----/----/----/----/----/----/----/----| 
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
  it=1
  ip=13
  wd='./weight-dryTE-jma-'//yyyymmddhh//'.grd'
  open(21,file=wd,status='new',access='direct',&
       &        convert='big_endian',&
       &        form='unformatted', recl=4*mem)
  
  ! 配列の割付
  l=min(mem,narea*nvar)
  print*,l
  allocate(ze(nlon,nlat,nvar,mem))
  allocate(z(narea*nvar,mem))
  allocate(zT(mem,narea*nvar))
  allocate(a8(narea*nvar,mem))
  allocate(sg8(l))
  allocate(sg(l))
  allocate(u8(narea*nvar,narea*nvar))  
  allocate(vt8(mem,mem))
  allocate(vt(mem,mem))
  allocate(p(mem))
  allocate(v(mem,mem))
  allocate(vtv(mem,mem))
  allocate(work(3*narea*nvar))
  
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
              T=zv3(:,:,1:3)
              print*,T(1,1,1)
           elseif(mod(id,4)==2)then
              ug=zv3(:,:,1:3)
              print*,ug(1,1,1)
           elseif(mod(id,4)==3)then
              vg=zv3(:,:,1:3)
              print*,vg(1,1,1)
           else
              q=zv3(:,:,1:3)
              print*,q(1,1,1)
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
        print*,imem!ze(1,1,:,imem)
     endif
  enddo
  
  ilt=0
  ilu=0
  ilv=0
  ilq=0
  rdf=dir//yyyy//'/jma/100900_mean.nc'
  inquire(file=rdf, exist=ex)
  if(ex)then
     do id=1,4
        !print*,rdf
        !print*,vname(id)
        call fread3(rdf,vname(id),ip,zv3)
        !print*,maxval(zv),minval(zv)
        if(mod(id,4)==1)then
           T=zv3(:,:,1:3)
           print*,T(1,1,1)
        elseif(mod(id,4)==2)then
           ug=zv3(:,:,1:3)
           print*,ug(1,1,1)
        elseif(mod(id,4)==3)then
           vg=zv3(:,:,1:3)
           print*,vg(1,1,1)
        else
           q=zv3(:,:,1:3)
           print*,q(1,1,1)
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
  do j=1,nlat
     lat=dslat+(j-1)-1
     ze(:,j,:,:)=ze(:,j,:,:)*sqrt(cos(lat*dtheta*pi/180.0))
  enddo
      
  do ilev=1,3            !850,500,300hPa
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
  
  a8=real(z,kind=8)
  call dgesvd('o','a',narea*nvar,mem,a8,narea*nvar,sg8,u8,narea*nvar,vt8,mem,work,lwork,info)
  !特異値分解(eof.f90)と異なり、固有値の小さい順にデータが格納されていることに注意
  print *,'info=',info
  !print *,'eigen values=',s(:) !固有値
  !do i=l,1,-1
  !   sg(l-i+1)=sqrt(mem*s(i))
  !enddo
  !print*,"0"
  !call calc_svd(z,int(narea*nvar),mem,sg,vt)

  vt=real(vt8,kind=4)
  sg=real(sg8,kind=4)
  
  print *,'singular values (double)=', sg8(1:10) !特異値(sigma=sqrt(m*s))
  print *,'singular values=', sg(1:10) !特異値(sigma=sqrt(m*s))
  write(21,rec=it) sg
  it=it+1

  tmp=0
  do ieof=1,l 
     tmp=tmp+sg(ieof)*sg(ieof) 
  end do
      
  !s=sg
  crate=0.0
  !do ieof=l,l-9,-1          !固有値の大きい順に直して表示
  do ieof=1,10
     print *,'-- PC',ieof
     print *,'proportion=',sg(ieof)*sg(ieof)/tmp
     crate=crate+sg(ieof)*sg(ieof)/tmp*100
     print *,'contribution rate(%)=',crate
     print *,'vector=',vt(ieof,:)
     print *,'points:'
     do i=1,10
        score = 0.0
        do j=1,mem
           score = score + z(i,j)*vt(ieof,j)
        enddo
        print *,i,score
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
  deallocate(ze,z,zT,a8,sg8,u8,vt8,vt,sg,p,v,vtv,work)
    
end program ensvsa_TE
