module read_qse
use hdf5
use hdf5_utils
use iso_fortran_env, only: dp => real64
use linspace_mod

implicit none

contains

subroutine data_qse(y)
  ! This subroutine reads params.jld
  ! add path to filename
  implicit none

  real(kind=8), intent(inout)      :: y(75,1,1,75,5371)
  integer, parameter               :: ydims = 5
  integer(HID_T)                   :: yshape(ydims)
  integer                          :: errorflag
  integer(HID_T)                   :: file_id,group_id,root_id, d_id
  character(len=60)                :: dataname, rootname, filename

  filename = "QSE_table.jld"
  rootname = "/"
  dataname = "/"//"data"
  call open_hdf5file(filename,file_id,errorflag)
  call h5gopen_f(file_id, rootname, root_id, errorflag)
  if (errorflag .lt. 0) then
          print*, "Error root group"
          return
  end if
  call h5dopen_f(root_id, dataname, d_id, errorflag)
  if (errorflag .lt. 0) then
          print*, "Error data open"
          return
  end if
  call H5DREAD_F(d_id, H5T_NATIVE_DOUBLE, y, yshape, errorflag)
  if (errorflag .lt. 0) then
          print*, "Error data open"
          return
  end if
  call h5dclose_f(d_id, errorflag)
  call h5gclose_f(root_id,errorflag)
  call h5fclose_f(file_id,errorflag)
  call h5close_f(errorflag)
  if (errorflag .lt. 0) then
          print*, "Error closing files"
          return
  end if
end subroutine data_qse

subroutine range_params(yrange, dataname_y)
  ! This subroutine reads params.jld
  ! add path to filename
  ! can be called for srange as well
  real(kind=8), intent(inout) :: yrange(75)
  character(len=7), intent(in) :: dataname_y
  integer, parameter          :: ydims = 1
  integer(HID_T)              :: yshape(ydims)
  integer                     :: errorflag
  integer(HID_T)              :: file_id,group_id,root_id, d_id
  character(len=60)           :: filename, rootname

  yshape = shape(yrange)
  filename = "params.jld"
  rootname = "/"
  !dataname_y = dataname_y
  call open_hdf5file(filename,file_id,errorflag)
  call h5gopen_f(file_id, rootname, root_id, errorflag)
  if (errorflag .lt. 0) then
          print*, "Error root group"
          return
  end if
  call h5dopen_f(root_id, dataname_y, d_id, errorflag)
  if (errorflag .lt. 0) then
          print*, "Error data open"
          return
  end if
  call H5DREAD_F(d_id, H5T_NATIVE_DOUBLE, yrange, yshape, errorflag)
  if (errorflag .lt. 0) then
          print*, "Error data open"
          return
  end if
  call h5dclose_f(d_id, errorflag)
  call h5gclose_f(root_id,errorflag)
  call h5fclose_f(file_id,errorflag)
  call h5close_f(errorflag)
  if (errorflag .lt. 0) then
          print*, "Error closing files"
          return
  end if
end subroutine range_params

subroutine qse(rho,temp,ye,x,x_cl,enbyrst)
      implicit none
      real(kind=8), intent(in)    :: rho, temp, ye, x_cl,enbyrst
      real(kind=8), intent(inout) :: x(5371)
      real(kind=8)                :: y(75,1,1,75,5371)
      real(kind=8)                :: yrange(75), rrange(64),trange(75),srange(75)
      !call range_params(yrange, "/yrange")
      !call range_params(srange, "/yrange")
      !call data_qse(y)
      !print*, yrange
      trange = linspace(3e9_dp,5e9_dp,75) !trange
      rrange =  linspace(1e8_dp,5e9_dp,64) !rrange
      print*, trange

end subroutine qse

subroutine interpolate4D(rho,tem,ye,cl)
! f11: rho(i_r),T(i_t)
! f12: rho(i_r), T(i_t+1)
! f21: rho(i_r+1),T(i_t)
! f22: rho(i_r+1),T(i_t+1)
! i1: ye_i, i2: rho_i, i3: x_cl_i, i4: t_i
implicit none
real(kind=8), intent(in) :: rho,tem,ye,cl
real(kind=8)     :: yrange(75), rrange(64),trange(75),srange(75)
integer, intent(in) :: i_ye, i_cl
real(kind=8), intent(out) :: x_inter
real(kind=8) :: x1,x2,x3,x4,&
              & x5,x6,x7,x8,&
              & x9,x10,x11,x12,&
              & x13,x14,x15,x16
real(kind=8) :: w10,w20,w30,w40,&
              & w01,w21,w31,w41
!convert dir name to rho,tem
! find i_ye, i_cl in those 4 files
rrange = linspace(1e8_dp,5e9_dp,64)
trange = linspace(3e9_dp,5e9_dp,75)
call range_params(yrange, "/yrange")
call range_params(srange, "/srange")
minloc(abs(a-rho / b(idx2)), 64,mask=())

! (rho,temp,ye,xcl)
! xi is a 1D value of element i:
x1 = f11(i_ye,i_cl)      ! 0000
x2 = f22(i_ye+1,i_cl+1)  ! 1111
x3 = f21(i_ye,i_cl)      ! 1000
x4 = f12(i_ye,i_cl)      ! 0100
x5 = f11(i_ye+1,i_cl)    ! 0010
x6 = f11(i_ye,i_cl+1)    ! 0001
x7 = f22(i_ye,i_cl)      ! 1100
x8 = f11(i_ye+1,i_cl+1)  ! 0011
x9 = f21(i_ye+1,i_cl)    ! 1010
x10 = f12(i_ye,i_cl+1)   ! 0101
x11 = f12(i_ye+1,i_cl)   ! 0110
x12 = f21(i_ye,i_cl+1)   ! 1001
x13 = f22(i_ye+1,i_cl)   ! 1110
x14 = f22(i_ye,i_cl+1)   ! 1101
x15 = f21(i_ye+1,i_cl+1) ! 1011
x16 = f12(i_ye+1,i_cl+1) ! 0111
! 12 weights:
w10 = (rrange(i+1) - r) / (rrange(i+1) - rrange(i))
w21 = 1. - w10
w20 = (trange(i+1) - t) / (trange(i+1) - trange(i))
w21 = 1. - w20
w30 = (yrange(i+1) - yel) / (yrange(i+1) - yrange(i))
w31 = 1. - w30
w40 = (srange(i+1) - cl) / (srange(i+1) - srange(i))
w41 = 1. - w40

x_inter
end subroutine interpolate4D

subroutine binary_search(x,x1,i1,i2)
  real,              intent(in) :: x(:)
 real,              intent(in) :: x1
 integer                       :: i1,i2
 integer :: i0, i1

 ! starting interval is the whole array
 i0 = 1
 i1 = size(x)

 ! test if x1 is outside of interval [x(1), x(end)]
 if (x1 <= x(i0)) then
   i = i0
   return
 end if
 if (x1 >= x(i1)) then
   i = i1
   return
 end if

 ! binary search
 do while (i1 > (i0 + 1))
   i = (i1 + i0) / 2
   if      (x(i) < x1) then
     i0 = i
   else if (x(i) > x1) then
     i1 = i
   else
     return
   end if
 end do

 ! pick index of value that is closer
 i = merge(i0, i1, ((2 * x1) < (x(i0) + x(i1))))
end subroutine binary_search

end module read_qse
