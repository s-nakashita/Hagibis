module read_netcdf
  implicit none
  include '/opt/local/include/netcdf.inc'

  integer,parameter :: imax=181,jmax=161,kmax=3,ntime=21
  integer,parameter :: imaxa=720,jmaxa=361,kmaxa=3,ntimea=8

contains

  subroutine fread(fname,vname,ip,z)
    implicit none
    integer :: istart(3),icount(3)
    integer :: ic,ncid,idvar,idlon,idlat,idtime
    integer,intent(in) :: ip
    double precision :: rlon(0:imax-1),rlat(jmax),rtime(ntime)
    real,intent(out) :: z(0:imax-1,jmax)
    character,intent(in) :: fname*100,vname*17

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'longitude',idlon)
    ic=NF_INQ_VARID(ncid,'latitude',idlat)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_DOUBLE(ncid,idlon,rlon)
    ic=NF_GET_VAR_DOUBLE(ncid,idlat,rlat)
    ic=NF_GET_VAR_DOUBLE(ncid,idtime,rtime)

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
    double precision :: rlon(0:imax-1),rlat(jmax),rtime(ntime)
    real :: rlev(kmax)
    real,intent(out) :: z(0:imax-1,jmax,kmax)
    character,intent(in) :: fname*100,vname*17

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'longitude',idlon)
    ic=NF_INQ_VARID(ncid,'latitude',idlat)
    ic=NF_INQ_VARID(ncid,'level',idlev)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_DOUBLE(ncid,idlon,rlon)
    ic=NF_GET_VAR_DOUBLE(ncid,idlat,rlat)
    ic=NF_GET_VAR_REAL(ncid,idlev,rlev)
    ic=NF_GET_VAR_DOUBLE(ncid,idtime,rtime)


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
!
! subroutine for GSM initial file
!
  subroutine freada(fname,vname,ip,z,slon,elon,slat,elat)
    implicit none
    integer :: istart(3),icount(3)
    integer :: ic,ncid,idvar,idlon,idlat,idtime
    integer,intent(in) :: ip
    double precision,intent(in) :: slon,elon,slat,elat
    integer :: islon,ielon,islat,ielat
    double precision :: rlon(0:imaxa-1),rlat(jmaxa),rtime(ntimea)
    real,intent(out) :: z(0:imax-1,jmax)
    character,intent(in) :: fname*100,vname*18

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
    character,intent(in) :: fname*100,vname*17

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
!
! calcurate time steps between 2 dates
!
  subroutine calc_steps(idate,edate,nt)
    implicit none
    integer,intent(in) :: idate,edate !format:yyyymmddhh
    integer,intent(out) :: nt
    integer,parameter :: dt=6 ! timestep=6h
    integer :: diff,ny,nm,nd,nh,maxday
    integer :: iy,im,id,ih
    integer :: ey,em,ed,eh

    diff = idate
    iy = diff/1000000
    diff = diff - iy*1000000
    im = diff/10000
    diff = diff - im*100000
    id = diff/100
    diff = diff - id*100
    ih = diff/dt

    diff = edate
    ey = diff/1000000
    diff = diff - ey*1000000
    em = diff/10000
    diff = diff - em*100000
    ed = diff/100
    diff = diff - ed*100
    eh = diff/dt

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
      nh = (eh - ih)/dt
    else
      nh = (eh+24 - ih)/dt
      nd = nd - 1
    endif
    nt = nh + nd*4 + nm*4*30 + ny*4*365
    return
  end subroutine calc_steps
!
! calcurate specific humidity from relative humidity
!
  SUBROUTINE calc_q(t,rh,p,q)
    IMPLICIT NONE
    double precision,PARAMETER :: t0=273.15d0
    double precision,PARAMETER :: e0c=6.11d0
    double precision,PARAMETER :: al=17.3d0
    double precision,PARAMETER :: bl=237.3d0
    double precision,PARAMETER :: e0i=6.1121d0
    double precision,PARAMETER :: ai=22.587d0
    double precision,PARAMETER :: bi=273.86d0
    real,INTENT(IN) :: t,rh,p !rh=%,p=hPa
    real,INTENT(OUT) :: q
    real :: e,es,tc

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
end module read_netcdf
