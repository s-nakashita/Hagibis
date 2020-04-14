module read_netcdf
  implicit none
  include '/opt/local/include/netcdf.inc'

  integer,parameter :: imax=61,jmax=161,kmax=5,ntime=29

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

end module read_netcdf
