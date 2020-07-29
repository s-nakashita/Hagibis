program reg_prtb

  use read_netcdf
  use stat
  
  implicit none
 
  integer,parameter :: dslon=95, delon=105, dslat=67, delat=75 
  integer,parameter :: nlon=delon-dslon+1, nlat=delat-dslat+1 
  integer,parameter :: narea=nlon*nlat 
  integer,parameter :: nv3d=4,nv2d=1,nlev=3
  integer,parameter :: nvar=nv3d*nlev+nv2d
  integer,parameter :: mem=26 
  real,parameter :: dtheta=0.5, pi=atan(1.0)*4.0 
  real,parameter :: cp=1005.7, R=287.04, Lh=2.5104*10**6 
  real,parameter :: Tr=270.0, pr=1000.0 
  real,parameter :: thres=1.782 !<=t((26-13-1),0.05)
      
  integer :: i,j,k,l,n
  integer :: imem,id,irec
  integer :: ilt,ilu,ilv,ilq,ivar,ip,fday
  logical :: ex
  
  real ::  ps(imax,jmax),ug(imax,jmax,nlev),vg(imax,jmax,nlev)
  real ::  T(imax,jmax,nlev),q(imax,jmax,nlev),rh(imax,jmax,nlev)
  real ::  zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)
  real ::  z0(imax,jmax,nvar)
  real ::  ze(imax,jmax,nvar,0:mem)
  real ::  plev(nlev)
  data plev/850.0,500.0,350.0/
  
  real :: err(0:mem)
  integer :: em
  real :: a,b(nvar),cor(nvar)
  real :: thr(nvar)
  real :: R2
  real :: x(nvar,mem),y(mem)

  real :: v3d(imax,jmax,nvar,3) !b,cor,thr(ug,vg,T,q,slp)
  real :: v2d(imax,jmax) !R2
  real :: buf4(imax,jmax)
      
  character rdf*100,rdw*100,wd*100
  character dir*30,dira*33,nmem*2,yyyy*4,mm*2,mmddhh*6,yyyymmddhh*10
  character(len=17) :: vname(5)
  character(len=4) :: vnamea(5)
  data vname/'UGRD','VGRD','TMP','SPFH','PRES_meansealevel'/
  data vnamea/'air','uwnd','vwnd','shum','slp'/
     !|----/----/----/----/----/----/----/----/----/----| 
  dir='/Users/nakashita/netcdf/tigge/'
 dira='/Users/nakashita/netcdf/nc-reanl/'

!  データの設定 
   yyyymmddhh="2019100912"
   yyyy=yyyymmddhh(1:4)
   mm=yyyymmddhh(5:6)
   mmddhh=yyyymmddhh(5:10)
   print*,yyyy,mmddhh
   rdw='../pytrack/error'//yyyymmddhh//'-jma.txt'
   open(10,file=rdw,status='old')
   do i=0,mem
      read(10,'(I2.2,f8.3)') em,err(i)
      print *, em,err(i)
   enddo
   y = err(1:mem)

   close(10)

   irec=1
   !do fday=0,6 !every 12 hours
   ip=1
   ilt=0
   ilu=0
   ilv=0
   ilq=0
   rdf=dir//yyyy//'/jma/'//mmddhh//'_mean.nc'
   !print*,rdf
   inquire(file=rdf, exist=ex)
   if(ex)then
      do id=1,4
         call fread3(rdf,vname(id),ip,zv3)
         if(mod(id,4)==1)then
            ug=zv3(:,:,1:3)
            !print*,T(1,1,1)
         elseif(mod(id,4)==2)then
            vg=zv3(:,:,1:3)
            !print*,ug(1,1,1)
         elseif(mod(id,4)==3)then
            T=zv3(:,:,1:3)
            !print*,vg(1,1,1)
         else
            q=zv3(:,:,1:3)
            do k=1,nlev
               do j=1,jmax
                  do i=1,imax
                     call calc_rh(T(i,j,k),q(i,j,k),plev(k),rh(i,j,k))
                  enddo
               enddo
            enddo
            !print*,q(1,1,1)
         endif
      enddo
           
      call fread(rdf,vname(5),ip,zv)
      ps=zv/100     !Pa->hPa
           
      ze(:,:,1:3,0)=ug
      ze(:,:,4:6,0)=vg
      ze(:,:,7:9,0)=T
      ze(:,:,10:12,0)=rh
      ze(:,:,13,0)=ps
      !ze(:,:,10,imem)=ps
   endif
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
               ug=zv3(:,:,1:3)
               !print*,T(1,1,1)
            elseif(mod(id,4)==2)then
               vg=zv3(:,:,1:3)
               !print*,ug(1,1,1)
            elseif(mod(id,4)==3)then
               T=zv3(:,:,1:3)
               !print*,vg(1,1,1)
            else
               q=zv3(:,:,1:3)
               do k=1,nlev
                  do j=1,jmax
                     do i=1,imax
                        call calc_rh(T(i,j,k),q(i,j,k),plev(k),rh(i,j,k))
                     enddo
                  enddo
               enddo
               !print*,q(1,1,1)
            endif
         enddo
           
         call fread(rdf,vname(5),ip,zv)
         ps=zv/100     !Pa->hPa
           
         ze(:,:,1:3,imem)=ug
         ze(:,:,4:6,imem)=vg
         ze(:,:,7:9,imem)=T
         ze(:,:,10:12,imem)=rh
         ze(:,:,13,imem)=ps
         !ze(:,:,10,imem)=ps
           
         !print*,ze(1,1,:,imem)
      endif
   enddo
     
   ip=3 !analysis start from 2019-10-09T00
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
      do id=1,4
         !print*,rdf
         !print*,vname(id)
         call fread3(rdf,vname(id),ip,zv3)
         !call fread3a(rdf,vnamea(id),ip,zv3,90.0d0,180.0d0,0.0d0,80.0d0)
         !print*,maxval(zv),minval(zv)
         if(mod(id,4)==1)then
            ug=zv3(:,:,1:3)
            !print*,T(1,1,1)
         elseif(mod(id,4)==2)then
            vg=zv3(:,:,1:3)
            !print*,ug(1,1,1)
         elseif(mod(id,4)==3)then
            T=zv3(:,:,1:3)
            !print*,vg(1,1,1)
         else
            q=zv3(:,:,1:3)
            !rh=zv3(:,:,1:3)
            !print*,"rh",rh(1,1,1)
            !do k=1,kmax
            !   do j=1,jmax
            !      do i=1,imax
            !         call calc_q(T(i,j,k),rh(i,j,k),plev(k),q(i,j,k))
            !      enddo
            !   enddo
            !enddo
            !print*,q(1,1,1)
         endif
      enddo
     
      rdf=dir//yyyy//'/jma/anl.nc' !_a
      call fread(rdf,vname(5),ip,zv)
      !call freada(rdf,vnamea(5),ip,zv,90.0d0,180.0d0,0.0d0,80.0d0)
      ps=zv/100        !Pa->hPa
      !print*,ps(1,1)
      !print*,maxval(ps),minval(ps)
        
      z0(:,:,1:3)=ug
      z0(:,:,4:6)=vg
      z0(:,:,7:9)=T
      z0(:,:,10:12)=q
      z0(:,:,13)=ps
      !z0(:,:,10)=ps
   endif
   !print*,z0(1,1,:)
     
     
   !1.calcurate perturbation
   do imem=1,mem
      ze(:,:,:,imem)=ze(:,:,:,imem)-ze(:,:,:,0)!z0(:,:,:)
   enddo

   !2.regression one by one point
   do j = 1, jmax
      do i = 1, imax
         x(:,:) = ze(i,j,:,1:mem)
         call reg_mul(nvar,mem,x,y,.TRUE.,a,b,cor,R2,thr)
         if(maxval(cor)>1.0.or.minval(cor)<-1.0) print *, maxval(cor),minval(cor)
         v3d(i,j,:,1) = b
         v3d(i,j,:,2) = cor
         do ivar = 1, nvar
            if(abs(thr(ivar)) .lt. thres) thr(ivar)=0.0
         end do
         v3d(i,j,:,3) = thr
             v2d(i,j) = R2
      end do
   end do
  
   wd='./reg-mul-rh-jma-'//yyyymmddhh//'-gr'
   open(21,file=wd,status='new',access='direct',&
      &        convert='big_endian',&
      &        form='unformatted', recl=4*imax*jmax)
      
   irec = 1
   do l = 1,3
      do n = 1,nv3d
         do k = 1,nlev
            ivar = (n-1)*nlev + k
            buf4 = v3d(:,:,ivar,l)
            write(21,rec=irec) buf4
            irec = irec+1
         enddo
      enddo
      buf4 = v3d(:,:,nvar,l)
      write(21,rec=irec) buf4
      irec = irec+1
   enddo

   buf4 = v2d
   write(21,rec=irec) buf4
   close(21)
   
   stop  
    
end program reg_prtb
