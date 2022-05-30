program grads_ensvsa_dTE

  use read_netcdf
  
  implicit none
 
  !integer,parameter :: dslon=95, delon=105, dslat=67, delat=75 
  integer,parameter :: dslon=275, delon=285, dslat=247, delat=255 
  !integer,parameter :: dslon=271, delon=281, dslat=247, delat=255 
  integer,parameter :: nlon=delon-dslon+1, nlat=delat-dslat+1 
  integer,parameter :: narea=nlon*nlat 
  integer,parameter :: nv3d=3,nv2d=4,nlev=3
  integer,parameter :: nvar=nv3d*nlev+nv2d-3 
  !integer,parameter :: memo=50, memn=26 
  real,parameter :: dtheta=0.5, pi=atan(1.0)*4.0 
  real,parameter :: cp=1005.7, R=287.04, Lh=2.5104*10**6 
  real,parameter :: Tr=270.0, pr=1000.0 
      
  integer :: i,j,k,n
  integer :: imem,id,it,irec
  integer :: ilt,ilu,ilv,ilev,ivar,ip,fday,tday
  integer :: imode
  !integer :: mem
  !integer :: idate,edate
  logical :: ex
  
  real ::  ps(imax,jmax),ug(imax,jmax,nlev),vg(imax,jmax,nlev)
  real ::  T(imax,jmax,nlev)
  real ::  zv3(0:imax-1,jmax,kmax),zv(0:imax-1,jmax)
  real ::  z0(imax,jmax,nvar),zm(imax,jmax,nvar)
  real,allocatable ::  ze(:,:,:,:)
  real ::  sigma(3),ssg
  real :: plev(3)
  data plev/300.0,500.0,850.0/
  real,allocatable ::  z(:,:),zT(:,:)
  real,allocatable :: sg(:),p(:),w(:,:)
  real :: TE(imax,jmax),KE(imax,jmax),PE(imax,jmax)
  real :: v3d(imax,jmax,nlev,nv3d) !ug,vg,T
  real :: v2d(imax,jmax,nv2d) !ps,TE
  real :: buf4(imax,jmax)
      
  character(len=5) :: orig="jma"
  integer :: mem=26 
  integer :: idate=2019100912
  integer :: edate=2019101212
  integer :: smode=1
  integer :: emode=1
  namelist /sens_nml/ orig, mem, idate, edate, smode, emode

  character rdf*100,rdw*100,wd*100
  character dir*32,dira*33,ns*1,ne*1,nmem*2
  character yyyy*4,mm*2,mmddhh*6,yyyymmddhh*10,cedate*10
  character(len=3) :: vname(4)
  !character(len=4) :: vnamea(5)
  !data vname/'UGRD','VGRD','TMP','PRES_meansealevel'/
  data vname/'u','v','t','msl'/
  !data vnamea/'air','uwnd','vwnd','shum','slp'/
     !|----/----/----/----/----/----/----/----/----/----| 
  dir='/Volumes/dandelion/netcdf/tigge/'
 !dira='/Users/nakashita/netcdf/nc-reanl/'

   sigma(1)=200.0/pr
   sigma(2)=6.0/7.0*300.0/pr
   sigma(3)=8.0/7.0*300.0/pr
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
   print *, cedate
   !smode=1
   !emode=1
   write(ns,'(I1)') smode
   write(ne,'(I1)') emode
   if (smode==emode) then
      wd='./ensvsa-dTE-m'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_n-gr'
   else
      wd='./ensvsa-dTE-m'//ns//'-'//ne//'-'//trim(orig)//'-'//yyyymmddhh//'-'//cedate//'_n-gr'
   endif
   open(21,file=wd,status='new',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*imax*jmax)
  
   !mem=memn 
   rdw='./weight-dTE-'//trim(orig)//'-'//yyyymmddhh//'_n.grd'
   open(10,file=rdw,status='old',access='direct',&
          &        convert='big_endian',&
          &        form='unformatted', recl=4*mem)

  ! 配列の割付
   allocate(ze(imax,jmax,nvar,mem))
   allocate(z(narea*nvar,mem))
   allocate(zT(mem,narea*nvar))
   allocate(sg(smode:emode))
   allocate(p(mem))
   allocate(w(mem,smode:emode))

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
   !edate=2019100900
   read(yyyymmddhh,*) edate
   print*,"edate=",edate
   do fday=0,tday !every 6 hours
      !idate=2019100900
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
         rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_'//nmem//'.nc'
         !print*,rdf
         inquire(file=rdf, exist=ex)
         if(ex)then
            do id=1,3
               call fread3(rdf,vname(id),ip,zv3)
               if(mod(id,3)==1)then
                  !ug=zv3(:,:,1:3)
                  ug=zv3(:,:,3:)
                  !print*,ug(1,1,1)
               elseif(mod(id,3)==2)then
                  !vg=zv3(:,:,1:3)
                  vg=zv3(:,:,3:)
                  !print*,vg(1,1,1)
               else
                  !T=zv3(:,:,1:3)
                  T=zv3(:,:,3:)
                  !print*,T(1,1,1)
               endif
            enddo
           
            call fread(rdf,vname(4),ip,zv)
            ps=zv/100     !Pa->hPa
            !print*,ps(1,1)
           
            ze(:,:,1:3,imem)=ug
            ze(:,:,4:6,imem)=vg
            ze(:,:,7:9,imem)=T
            ze(:,:,10,imem)=ps
           
            !print*,ze(1,1,:,imem)
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
      rdf=dir//yyyy//'/'//trim(orig)//'/glb_'//yyyymmddhh//'_mean.nc'   !_n
      !rdf=dir//yyyy//'/jma/100900_mean.nc'
      !rdf=dir//yyyy//'/jma/anl_sellev.nc' !_a
      inquire(file=rdf, exist=ex)
      if(ex)then
         do id=1,3
            !print*,rdf
            !print*,vname(id)
            call fread3(rdf,vname(id),ip,zv3)
            !call fread3a(rdf,vnamea(id),ip,zv3,90.0d0,180.0d0,0.0d0,80.0d0)
            !print*,maxval(zv),minval(zv)
            if(mod(id,3)==1)then
               !ug=zv3(:,:,1:3)
               ug=zv3(:,:,3:)
               !print*,ug(1,1,1)
            elseif(mod(id,3)==2)then
               !vg=zv3(:,:,1:3)
               vg=zv3(:,:,3:)
               !print*,vg(1,1,1)
            else
               !T=zv3(:,:,1:3)
               T=zv3(:,:,3:)
               !print*,T(1,1,1)
            endif
         enddo
     
         !rdf=dir//yyyy//'/jma/anl.nc' !_a
         call fread(rdf,vname(4),ip,zv)
         !call freada(rdf,vnamea(5),ip,zv,90.0d0,180.0d0,0.0d0,80.0d0)
         ps=zv/100        !Pa->hPa
         !print*,ps(1,1)
         !print*,maxval(ps),minval(ps)
        
         z0(:,:,1:3)=ug
         z0(:,:,4:6)=vg
         z0(:,:,7:9)=T
         z0(:,:,10)=ps
      endif
      !print*,z0(1,1,:)
     
     
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
      do n = 1,nv3d
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
               
     !area=0.0
     !do i=1,nlon
     !   do j=1,nlat
     !      lat=dslat+(j-1)-1
     !      area=area+cos(lat*dtheta*pi/180.0)
     !      ze(i,j,:,:)=ze(i,j,:,:)*cos(lat*dtheta*pi/180.0)
     !   enddo
     !enddo
     
      do ilev=1,3            !850,500,300hPa
         do ivar=1,3         !ug,vg,T
            ze(:,:,3*(ivar-1)+ilev,:)=ze(:,:,3*(ivar-1)+ilev,:)*sigma(ilev)
         enddo
      enddo
     !3.Multiply by coefficient
     !T
      ze(:,:,7:9,:)=ze(:,:,7:9,:)*sqrt(cp/Tr)
     !ps
      ze(:,:,10,:)=ze(:,:,10,:)*sqrt(R*Tr)/pr
     
      !do imem=1,mem
      !   print*,imem
      !   print*,ze(1,1,:,imem)
      !enddo

     !4.calcurate energy
      TE=0.0
      KE=0.0
      PE=0.0
      zm=0.0
      do imode=smode,emode
         do imem=1,mem
            zm=zm+ze(:,:,:,imem)*w(imem,imode)*sg(imode)/ssg
         enddo
      enddo

      do ivar=1,nvar
         TE=TE+zm(:,:,ivar)**2/2
      enddo
      do ivar=1,6
         KE=KE+zm(:,:,ivar)**2/2
      enddo
      do ivar=7,nvar
         PE=PE+zm(:,:,ivar)**2/2
      enddo
      !print*,"max",maxval(TE),"min",minval(TE)

      v2d(:,:,2) = TE
      v2d(:,:,3) = KE
      v2d(:,:,4) = PE

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
  
   deallocate(ze,z,zT,sg,p,w)
  
   stop  
    
end program grads_ensvsa_dTE
