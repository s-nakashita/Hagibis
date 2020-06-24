program test_sp
    USE read_netcdf

    implicit none
    integer,parameter :: nlat=jmax,nlon=imax,nlev=kmax,nnt=nlev
    integer,parameter :: mmdab=(nlon+2)/2,mmdc=(nlon+1)/2
    integer,parameter :: lleng = 5*nlat*nlat*nlon,llsav=5*nlat*nlat*nlon
    double precision :: work(lleng),wsave(llsav)
    integer,parameter :: lldwork = 4*nlat*nlat
    double precision :: dwork(lldwork)
    double precision,parameter :: pi = 4.0d0*atan(1.0d0)
    double precision,parameter :: re = 6371.0d3
    double precision :: thetag(nlat),dtheta(nlat),dwts(nlat) 
    double precision :: br(mmdc,nlat,nnt),bi(mmdc,nlat,nnt)
    double precision :: cr(mmdc,nlat,nnt),ci(mmdc,nlat,nnt)
    double precision :: vt(nlat,nlon,nnt),dv(nlat,nlon,nnt)
    double precision :: sf(nlat,nlon,nnt),vp(nlat,nlon,nnt)
    double precision :: v(nlat,nlon,nnt),w(nlat,nlon,nnt)
    double precision :: pertrb(nnt)
    integer :: isym,ityp,ier
    integer :: nt,nmax,mdab,mdc,lwork,lsave,ldwork
    integer :: i,j
                          !0    5   10   15   20   25   30   35   40
                          !|----/----/----/----/----/----/----/----/
    character(100) :: infile="../../netcdf/gsm/gl/2019/10/init.nc"
    character(7) :: outfile="vor.grd"
    character(len=18) :: vname(5)
    data vname/'TMP','UGRD','VGRD','SPFH','PRMSL_meansealevel'/
    real :: ug(0:nlon-1,nlat,nlev),vg(0:nlon-1,nlat,nlev)
    double precision :: b3d(nlon,nlat,nlev,5) !vor,div,strm,vpot,Tbar
    real :: buf(nlon,nlat),Tbar(nlon,nlat,nlev)
    integer :: iolen,k,n,irec

    !call fread(infile,"PRMSL_meansealevel",0,buf)
    !print*,buf(1,1)
    call fread3(infile,vname(2),1,ug)
    call fread3(infile,vname(3),1,vg)
    print*, "ug =",ug(12,1:10,2)
    print*, "vg =",vg(12,1:10,2)

    nt = nnt
    nmax = max(nlat,nlon)
    mdab = mmdab
    mdc = mmdc
    lwork = lleng
    lsave = llsav
    ldwork = lldwork
    ityp=0
    isym=0
    print*,nlat
    
    do k=1,nlev
        call geo2mathv(0,nlon,nlat,real(ug(:,:,k),8),real(vg(:,:,k),8),v(:,:,k),w(:,:,k),work)
        !do j=1,nlat
        !    do i=1,nlon
        !        v(j,i,k)=-v3d(i,nlat+1-j,k,iv3d_v)
        !        w(j,i,k)=v3d(i,nlat+1-j,k,iv3d_u)
        !    enddo
        !enddo
    enddo
    print*, "vm =",v(12,1:10,2)
    print*, "wm =",w(12,1:10,2)

    call vhaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ier)
    print*, "ierr=",ier
    
    call vhaes(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,mdc,nlat,wsave,lsave,work,lwork,ier)
    print*, "ierr=",ier

    call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ier)
    print*, "ierr=",ier  

    call vrtes(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdc,nlat,wsave,lsave,work,lwork,ier)
    print*, "ierr=",ier
    vt = vt/re
    print*, "Vort=",vt(12,1:10,2)

    call dives(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdc,nlat,wsave,lsave,work,lwork,ier)
    print*, "ierr=",ier
    dv = dv/re
    print*, "Div=",dv(12,1:10,2)

    call sfvpes(nlat,nlon,isym,nt,sf,vp,nlat,nlon,br,bi,cr,ci,mdc,nlat,wsave,lsave,work,lwork,ier)
    print*, "ierr=",ier
    do j=1,nlat
        sf(j,:,:) = sf(j,:,:)*re!*fcori(j)
        vp(j,:,:) = vp(j,:,:)*re!*fcori(j)
    enddo
    print*, "strm=",sf(12,1:10,2)
    print*, "vpot=",vp(12,1:10,2)

    do k=1,nlev
        call math2geos(0,nlat,nlon,vt(:,:,k),b3d(:,:,k,1),work)
        call math2geos(0,nlat,nlon,dv(:,:,k),b3d(:,:,k,2),work)
        call math2geos(0,nlat,nlon,sf(:,:,k),b3d(:,:,k,3),work)
        call math2geos(0,nlat,nlon,vp(:,:,k),b3d(:,:,k,4),work)
        !do i=1,nlon
        !    do j=1,nlat
        !        v3d(i,j,k,1) = vt(nlat+1-j,i,k)/re
        !        v3d(i,j,k,2) = dv(nlat+1-j,i,k)/re
        !        v3d(i,j,k,3) = sf(nlat+1-j,i,k)/re/fcori(j)
        !        v3d(i,j,k,4) = vp(nlat+1-j,i,k)/re/fcori(j)
        !    enddo
        !enddo
    enddo

    OPEN(55,FILE=outfile,FORM='unformatted',ACCESS='direct',RECL=4*nlon*nlat)

    irec=1
    DO n=1,5
        DO k=1,nlev
            buf = REAL(b3d(:,:,k,n),4)
            WRITE(55,REC=irec) buf
            irec = irec + 1
        ENDDO
    ENDDO

    CLOSE(55)

end program