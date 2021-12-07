module netcdf_write

use grid
use mpi_interface
use netcdf

implicit none

integer :: ncid
integer :: u_varid,w_varid,bu_varid,bs_varid,p_varid,ql_varid

contains

subroutine init_par_vars(filename,cmode,comm,info,ncid)

implicit none

character :: filename
integer :: comm,cmode,info,ncid

call nf_create_par(filename,cmode,comm,info,ncid)

end subroutine 


subroutine init_netcdf_vars

implicit none

integer,parameter :: ndims=3
integer :: u_dimids(ndims),w_dimids(ndims),t_dimids(ndims),xm_dimid,zm_dimid,xt_dimid,zt_dimid,rec_dimid,nc_rec,x_dim
!character (len=*), parameter :: file_name="file_name"
!character (len=*), parameter :: u_name="u"
!character (len=*) :: w_name="w"
!character (len=*) :: bu_name="bu"
!character (len=*) :: bs_name="bs"

nc_rec = nf90_create(file_name,NF90_CLOBBER,ncid)
nc_rec = nf90_def_dim(ncid,"xm",N_x,xm_dimid)
nc_rec = nf90_def_dim(ncid,"zm",N_vert,zm_dimid)
nc_rec = nf90_def_dim(ncid,"xt",N_x,xt_dimid)
nc_rec = nf90_def_dim(ncid,"zt",N_vert,zt_dimid)
nc_rec = nf90_def_dim(ncid,"time",nf90_unlimited,rec_dimid)

u_dimids=(/zt_dimid,xm_dimid,rec_dimid/) 
w_dimids=(/zm_dimid,xt_dimid,rec_dimid/) 
t_dimids=(/zt_dimid,xt_dimid,rec_dimid/)

nc_rec = nf90_def_var(ncid,"u",nf90_double,u_dimids,u_varid)
nc_rec = nf90_def_var(ncid,"w",nf90_double,w_dimids,w_varid)
nc_rec = nf90_def_var(ncid,"bu",nf90_double,t_dimids,bu_varid)
nc_rec = nf90_def_var(ncid,"bs",nf90_double,t_dimids,bs_varid)
nc_rec = nf90_def_var(ncid,"p",nf90_double,t_dimids,p_varid)
nc_rec = nf90_def_var(ncid,"ql",nf90_double,t_dimids,ql_varid)

nc_rec = nf90_enddef(ncid)

end subroutine init_netcdf_vars


subroutine write_netcdf_vars(U,W,Bu,Bs,P,ql,k)

implicit none

integer :: start(3), ucounts(3), wcounts(3),tcounts(3),nc_rec,k,varid
real*8 :: time
real*8, dimension(:,:) :: U,W,Bu,Bs,P,ql
ucounts = (/N_vert,N_x,1/)
wcounts = (/N_vert,N_x,1/)
tcounts = (/N_vert,N_X,1/)
start = (/1,1,k/)
nc_rec = nf90_inq_varid(ncid,"u",varid)
nc_rec= nf90_put_var(ncid,varid,U,start,ucounts)
nc_rec = nf90_inq_varid(ncid,"w",varid)
nc_rec= nf90_put_var(ncid,varid,W,start,wcounts)
nc_rec = nf90_inq_varid(ncid,"bu",varid)
nc_rec = nf90_put_var(ncid,varid,Bu,start,tcounts)
nc_rec = nf90_inq_varid(ncid,"bs",varid)
nc_rec = nf90_put_var(ncid,varid,Bs,start,tcounts)
nc_rec = nf90_inq_varid(ncid,"p",varid)
nc_rec = nf90_put_var(ncid,varid,P,start,tcounts)
nc_rec = nf90_inq_varid(ncid,"ql",varid)
nc_rec = nf90_put_var(ncid,varid,ql,start,tcounts)
!nc_rec = nf90_close(ncid)


end subroutine write_netcdf_vars

subroutine end_netcdf

implicit none

integer :: nc_rec

nc_rec = nf90_close(ncid)

end subroutine end_netcdf

end module
