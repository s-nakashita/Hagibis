module read_netcdf
  implicit none
  include '/opt/local/include/netcdf.inc'

  !integer,parameter :: imax=720,jmax=361,kmax=5,ntime=21
#ifdef lev6
  integer,parameter :: imax=720,jmax=361,kmax=9,ntime=21
#else
  integer,parameter :: imax=720,jmax=361,kmax=5,ntime=21
#endif
  integer,parameter :: imaxa=720,jmaxa=361,kmaxa=3,ntimea=15
  double precision, parameter :: dlon=0.5d0, dlat=0.5d0
  double precision, parameter :: lon0=0.0d0, lat0=-90.0d0

contains

  subroutine fread(fname,vname,ip,z)
    implicit none
    integer :: istart(3),icount(3)
    integer :: ic,ncid,idvar,idlon,idlat,idtime
    integer,intent(in) :: ip
    real :: rlon(0:imax-1),rlat(jmax)
    integer :: rtime(ntime)
    real,intent(out) :: z(0:imax-1,jmax)
    character,intent(in) :: fname*100,vname*3

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'longitude',idlon)
    ic=NF_INQ_VARID(ncid,'latitude',idlat)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_REAL(ncid,idlon,rlon)
    ic=NF_GET_VAR_REAL(ncid,idlat,rlat)
    ic=NF_GET_VAR_INT(ncid,idtime,rtime)

    istart(1) = 1
    istart(2) = 1
    istart(3) = ip

    icount(1) = imax
    icount(2) = jmax
    icount(3) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine fread

  subroutine fread3(fname,vname,ip,z)
    implicit none
    integer :: istart(4),icount(4)
    integer :: ic,ncid,idvar,idlon,idlat,idlev,idtime
    integer,intent(in) :: ip
    real :: rlon(0:imax-1),rlat(jmax),rlev(kmax)
    integer :: rtime(ntime)
    real,intent(out) :: z(0:imax-1,jmax,kmax)
    character,intent(in) :: fname*100,vname*2

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'longitude',idlon)
    ic=NF_INQ_VARID(ncid,'latitude',idlat)
    ic=NF_INQ_VARID(ncid,'level',idlev)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_REAL(ncid,idlon,rlon)
    ic=NF_GET_VAR_REAL(ncid,idlat,rlat)
    ic=NF_GET_VAR_REAL(ncid,idlev,rlev)
    ic=NF_GET_VAR_INT(ncid,idtime,rtime)

    !print *, rlev(:)
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = ip

    icount(1) = imax
    icount(2) = jmax
    icount(3) = kmax
    icount(4) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine fread3
!-----------------------------------------------------------------------
! Calculate boundary indexes
!-----------------------------------------------------------------------
  subroutine set_bound(slon,elon,slat,elat,islon,ielon,islat,ielat)
    implicit none 
    double precision,intent(in) ::  slon, elon, slat, elat
    integer,intent(out)         :: islon,ielon,islat,ielat
    islon = CEILING((slon-lon0)/dlon) + 1
    ielon = CEILING((elon-lon0)/dlon) + 1
    islat = CEILING((slat-lat0)/dlat) + 1
    ielat = CEILING((elat-lat0)/dlat) + 1
    return
  end subroutine set_bound
!-----------------------------------------------------------------------
! Subroutine for NCEP-NCAR reanalysis data (regrid 720x361)
!-----------------------------------------------------------------------
  subroutine freada(fname,vname,ip,z,slon,elon,slat,elat)
    implicit none
    integer :: istart(3),icount(3)
    integer :: ic,ncid,idvar,idlon,idlat,idtime
    integer,intent(in) :: ip
    double precision,intent(in) :: slon,elon,slat,elat
    integer :: islon,ielon,islat,ielat
    double precision :: rlon(0:imaxa-1),rlat(jmaxa),rtime(ntimea)
    real,intent(out) :: z(0:imax-1,jmax)
    character,intent(in) :: fname*100,vname*4

    islon = int(slon/0.5d0) + 1
    ielon = int(elon/0.5d0) + 1
    islat = int((slat+90.0d0)/0.5d0) + 1
    ielat = int((elat+90.0d0)/0.5d0) + 1
    print *, islon,ielon,islat,ielat

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'lon',idlon)
    ic=NF_INQ_VARID(ncid,'lat',idlat)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_DOUBLE(ncid,idlon,rlon)
    ic=NF_GET_VAR_DOUBLE(ncid,idlat,rlat)
    ic=NF_GET_VAR_DOUBLE(ncid,idtime,rtime)
    !print*,rlon(islon:ielon)
    !print*,rlat(islat:ielat)

    istart(1) = islon
    istart(2) = islat
    istart(3) = ip

    icount(1) = (ielon-islon)+1
    icount(2) = (ielat-islat)+1
    icount(3) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine freada

  subroutine fread3a(fname,vname,ip,z,slon,elon,slat,elat)
    implicit none
    integer :: istart(4),icount(4)
    integer :: ic,ncid,idvar,idlon,idlat,idlev,idtime
    integer,intent(in) :: ip
    double precision,intent(in) :: slon,elon,slat,elat
    integer :: islon,ielon,islat,ielat
    double precision :: rlon(0:imaxa-1),rlat(jmaxa),rtime(ntimea)
    real :: rlev(kmaxa)
    real,intent(out) :: z(0:imax-1,jmax,kmax)
    character,intent(in) :: fname*100,vname*4

    islon = int(slon/0.5d0) + 1
    ielon = int(elon/0.5d0) + 1
    islat = int((slat+90.0d0)/0.5d0) + 1
    ielat = int((elat+90.0d0)/0.5d0) + 1
    print *, islon,ielon,islat,ielat

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'lon',idlon)
    ic=NF_INQ_VARID(ncid,'lat',idlat)
    ic=NF_INQ_VARID(ncid,'level',idlev)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_DOUBLE(ncid,idlon,rlon)
    ic=NF_GET_VAR_DOUBLE(ncid,idlat,rlat)
    ic=NF_GET_VAR_REAL(ncid,idlev,rlev)
    ic=NF_GET_VAR_DOUBLE(ncid,idtime,rtime)
    !print *, rlon


    istart(1) = islon
    istart(2) = islat
    istart(3) = 1
    istart(4) = ip

    icount(1) = (ielon-islon)+1
    icount(2) = (ielat-islat)+1
    icount(3) = kmax
    icount(4) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine fread3a  
!-----------------------------------------------------------------------
! Calculate time steps between 2 dates
!-----------------------------------------------------------------------
  subroutine calc_steps(idate,edate,dt,nt)
    implicit none
    integer,intent(in) :: idate,edate !format:yyyymmddhh
    integer,intent(in) :: dt ! timestep
    integer,intent(out) :: nt
    integer :: diff,ny,nm,nd,nh
    integer :: iy,im,id,ih
    integer :: ey,em,ed,eh

    print *, dt
    diff = idate
    iy = diff/1000000
    diff = diff - iy*1000000
    im = diff/10000
    diff = diff - im*10000
    id = diff/100
    diff = diff - id*100
    ih = diff
    print*, ih,id,im,iy

    diff = edate
    ey = diff/1000000
    diff = diff - ey*1000000
    em = diff/10000
    diff = diff - em*10000
    ed = diff/100
    diff = diff - ed*100
    eh = diff
    print*, eh,ed,em,ey

    ny = ey - iy
    if(em.ge.im)then
      nm = em - im
    else
      nm = em+12 - im
      ny = ny - 1
    endif
    if(ed.ge.id)then
      nd = ed - id
    else
      if(em.eq.3)then
        if(mod(ey,4).ne.0)then
          nd = ed+28 - id
          nm = nm - 1
        else
          nd = ed+29 - id
          nm = nm -1
        endif
      elseif(em.eq.1.or.em.eq.2.or.em.eq.4.or.em.eq.6.or.em.eq.9.or.em.eq.11)then
        nd = ed+31 - id
        nm = nm - 1
      else
        nd = ed+30 - id
        nm = nm - 1
      endif
    endif
    if(eh.ge.ih)then
      nh = ceiling(real(eh - ih))/dt
    else
      nh = ceiling(real(eh+24 - ih))/dt
      nd = nd - 1
    endif
    print*, nh,nd,nm,ny
    nt = nh + nd*(24/dt) + nm*(24/dt)*30 + ny*(24/dt)*365
    return
  end subroutine calc_steps
!-----------------------------------------------------------------------
! Calculate date after specified timesteps
!-----------------------------------------------------------------------
  subroutine calc_date(idate,nt,edate)
    implicit none
    integer,intent(in) :: idate, nt !format:yyyymmddhh
    integer,intent(out) :: edate
    integer :: dt=6 !default timestep(6h)
    integer :: ey,em,ed,eh

    edate = idate + dt*nt
    eh = mod(edate,100)
    if (eh==24) then
      eh = 0
      edate = edate + 100 - 24
    endif
    ed = mod(edate,10000)/100
    em = mod(edate,1000000)/10000
    ey = (edate - em*10000 - ed*100 - eh)/1000000
    if (ed.gt.31.and.&
    &(em.eq.1.or.em.eq.3.or.em.eq.5.or.em.eq.7&
    &.or.em.eq.8.or.em.eq.10.or.em.eq.12))then
      edate = edate + 10000 - 3100
    elseif (ed.gt.30.and.&
      &(em.eq.4.or.em.eq.6.or.em.eq.9.or.em.eq.11))then
      edate = edate + 10000 - 3000
    elseif (em==2.and.ed.gt.28) then
      if (mod(ey,4)==0.and.ed.gt.29) then
        edate = edate + 10000 - 2900
      elseif(mod(ey,4).ne.0) then
        edate = edate + 10000 - 2800
      endif
    endif
    em = mod(edate,1000000)/10000
    if (em.gt.12) then
      edate = edate + 1000000 - 120000
    endif

    return
  end subroutine calc_date
!-----------------------------------------------------------------------
! Calculate specific humidity from relative humidity
!-----------------------------------------------------------------------
  SUBROUTINE calc_q(t,rh,p,q)
    IMPLICIT NONE
    double precision,PARAMETER :: t0=273.15d0
    double precision,PARAMETER :: e0c=6.11d0
    double precision,PARAMETER :: al=17.3d0
    double precision,PARAMETER :: bl=237.3d0
    double precision,PARAMETER :: e0i=6.1121d0
    double precision,PARAMETER :: ai=22.587d0
    double precision,PARAMETER :: bi=273.86d0
    double precision,INTENT(IN) :: t,rh,p !rh=%,p=hPa
    double precision,INTENT(OUT) :: q
    double precision :: e,es,tc

    tc = t-t0
    IF(tc >= 0.0d0) THEN
      es = e0c * exp(al*tc/(bl+tc))
    ELSE IF(tc <= -15.d0) THEN
      es = e0i * exp(ai*tc/(bi+tc))
    ELSE
      es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
        + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
    END IF

    e = rh * 0.01 * es
  
    q = e * 0.622d0 / (p - e * 0.378d0)
   
    RETURN
  END SUBROUTINE calc_q
!-----------------------------------------------------------------------
! Calculate relative humidity (RH)
!-----------------------------------------------------------------------
  SUBROUTINE calc_rh(t,q,p,rh)
    IMPLICIT NONE
    DOUBLE PRECISION,PARAMETER :: t0=273.15d0
    DOUBLE PRECISION,PARAMETER :: e0c=6.11d0
    DOUBLE PRECISION,PARAMETER :: al=17.3d0
    DOUBLE PRECISION,PARAMETER :: bl=237.3d0
    DOUBLE PRECISION,PARAMETER :: e0i=6.1121d0
    DOUBLE PRECISION,PARAMETER :: ai=22.587d0
    DOUBLE PRECISION,PARAMETER :: bi=273.86d0
    DOUBLE PRECISION,INTENT(IN) :: t,q,p
    DOUBLE PRECISION,INTENT(OUT) :: rh
    DOUBLE PRECISION :: e,es,tc
  
    e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)
  
    tc = t-t0
    IF(tc >= 0.0d0) THEN
      es = e0c * exp(al*tc/(bl+tc))
    ELSE IF(tc <= -15.d0) THEN
      es = e0i * exp(ai*tc/(bi+tc))
    ELSE
      es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
         + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
    END IF
  
    rh = e/es
  
    RETURN
  END SUBROUTINE calc_rh
end module read_netcdf
