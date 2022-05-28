module read_qse
use hdf5
use hdf5_utils
use iso_fortran_env, only: dp => real64
use linspace_mod

implicit none

contains

subroutine data_qse(y,path)
  ! This subroutine reads params.jld
  ! add path to filename
  implicit none

  real(kind=8), intent(inout)      :: y(75,1,1,75,5371)
  character(len=60), intent(in)    :: path
  integer, parameter               :: ydims = 5
  integer(HID_T)                   :: yshape(ydims)
  integer                          :: errorflag
  integer(HID_T)                   :: file_id,group_id,root_id, d_id
  character(len=60)                :: dataname, rootname, filename

  filename = "QSE_table.jld"
  !print*, filename
  !stop
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
      character(len=60) :: path
      !real(kind=8)                :: yrange(75), rrange(64),trange(75),srange(75)
      !integer :: j
      !call range_params(yrange, "/yrange")
      !call range_params(srange, "/yrange")
      !call data_qse(y,path)
      !stop
      !print*, yrange
      !trange = linspace(3e9_dp,5e9_dp,75) !trange
      !rrange =  linspace(1e8_dp,5e9_dp,64) !rrange
      call interpolate4D(1.e8_dp,3.e9_dp,0.49_dp,-3._dp,2,x(1))
      !print*, x(1)
end subroutine qse

subroutine interpolate4D(rho,tem,ye,cl,index_part,x_inter)
      ! f11: rho(i_r),T(i_t)
      ! f12: rho(i_r), T(i_t+1)
      ! f21: rho(i_r+1),T(i_t)
      ! f22: rho(i_r+1),T(i_t+1)

      implicit none
      real(kind=8), intent(in)  :: rho,tem,ye,cl
      integer, intent(in)       ::  index_part
      real(kind=8), intent(out) :: x_inter
      real(kind=8)              :: yrange(75), rrange(64),trange(75),srange(75)
      real(kind=8)              :: fl(75,1,1,75,5371)
      ! real(kind=8)              :: f11(75,1,1,75,5371)
      ! real(kind=8)              :: f22(75,1,1,75,5371)
      ! real(kind=8)              :: f12(75,1,1,75,5371)
      ! real(kind=8)              :: f21(75,1,1,75,5371)
      integer                   :: i_r,i_t,i_ye,i_cl
      real(kind=8)              :: x1,x2,x3,x4,&
                                   x5,x6,x7,x8,&
                                   x9,x10,x11,x12,&
                                   x13,x14,x15,x16
      real(kind=8)              :: w10,w20,w30,w40,&
                                   w11,w21,w31,w41
      character(len=60)         :: path1,path2,path3,path4

      rrange = linspace(1e8_dp,1e9_dp,64)
      trange = linspace(5e9_dp,3e9_dp,75)
      yrange = linspace(0.5_dp,0.47_dp,75)
      print*, yrange
      srange = linspace(-1.5._dp,-4._dp,75)
      path1 = ""
      path2 = ""
      path3 = ""
      path4 = ""
      !call range_params(yrange, "/yrange")
      !call range_params(srange, "/srange")
      !print*,"test"
      call binary_search(rrange,rho,i_r)
      call binary_search(trange,tem,i_t)
      call binary_search(yrange,ye,i_t)
      call binary_search(srange,cl,i_t)
      !print*,"test"
      ! now convert i_r,i_t into dir-names !!!
      call data_qse(fl,path1)
      print*,i_ye,i_cl
      !i_cl = 1
      !print*,fl(:,1,1,:,1)
      x1 = fl(i_ye,1,1,i_cl,1)      ! 0000
      x5 = fl(i_ye+1,1,1,i_cl,1)    ! 0010
      x6 = fl(i_ye,1,1,i_cl+1,1)    ! 0001
      x8 = fl(i_ye+1,1,1,i_cl+1,1)  ! 0011
      call data_qse(fl,path2)
      print*,"test"
      x4  = fl(i_ye,1,1,i_cl,1)     ! 0100
      x10 = fl(i_ye,1,1,i_cl+1,1)   ! 0101
      x11 = fl(i_ye+1,1,1,i_cl,1)   ! 0110
      x16 = fl(i_ye+1,1,1,i_cl+1,1) ! 0111
      call data_qse(fl,path3)
      x3  = fl(i_ye,1,1,i_cl,1)     ! 1000
      x9  = fl(i_ye+1,1,1,i_cl,1)   ! 1010
      x12 = fl(i_ye,1,1,i_cl+1,1)   ! 1001
      x15 = fl(i_ye+1,1,1,i_cl+1,1) ! 1011
      call data_qse(fl,path4)
      x2  = fl(i_ye+1,1,1,i_cl+1,1)  ! 1111
      x7  = fl(i_ye,1,1,i_cl,1)      ! 1100
      x13 = fl(i_ye+1,1,1,i_cl,1)    ! 1110
      x14 = fl(i_ye,1,1,i_cl+1,1)    ! 1101

! =================================================
! TODO: i_cl, i_cl not working because reverse:
! TODO: path directories einpflegen
! TODO: Loop fuer alle elements adden
! =================================================
! testing_QSE(ARGS,
! collect(LinRange(0.5,0.47,75)),
! collect(LinRange(5e9, 3e9, 75)),
! collect(LinRange(1e8, 1e9, 64)),
! collect(LinRange(-1.5, -4, 75)));

      ! (rho,temp,ye,xcl)
      ! xi is a 1D value of element i:
      ! x1 = f11(i_ye,i_cl)      ! 0000
      ! x2 = f22(i_ye+1,i_cl+1)  ! 1111
      ! x3 = f21(i_ye,i_cl)      ! 1000
      ! x4 = f12(i_ye,i_cl)      ! 0100
      ! x5 = f11(i_ye+1,i_cl)    ! 0010
      ! x6 = f11(i_ye,i_cl+1)    ! 0001
      ! x7 = f22(i_ye,i_cl)      ! 1100
      ! x8 = f11(i_ye+1,i_cl+1)  ! 0011
      ! x9 = f21(i_ye+1,i_cl)    ! 1010
      ! x10 = f12(i_ye,i_cl+1)   ! 0101
      ! x11 = f12(i_ye+1,i_cl)   ! 0110
      ! x12 = f21(i_ye,i_cl+1)   ! 1001
      ! x13 = f22(i_ye+1,i_cl)   ! 1110
      ! x14 = f22(i_ye,i_cl+1)   ! 1101
      ! x15 = f21(i_ye+1,i_cl+1) ! 1011
      ! x16 = f12(i_ye+1,i_cl+1) ! 0111
      ! 12 weights:
      w10 = (rrange(i_r+1) - rho) / (rrange(i_r+1) - rrange(i_r))
      w11 = 1. - w10
      w20 = (trange(i_t+1) - tem) / (trange(i_t+1) - trange(i_t))
      w21 = 1. - w20
      w30 = (yrange(i_ye+1) - ye) / (yrange(i_ye+1) - yrange(i_ye))
      w31 = 1. - w30
      w40 = (srange(i_ye+1) - cl) / (srange(i_ye+1) - srange(i_ye))
      w41 = 1. - w40

      x_inter = x1 * (w10 * w20 * w30 * w40) + &
                x2 * (w11 * w21 * w31 * w41) + &
                x3 * (w11 * w20 * w30 * w40) + &
                x4 * (w10 * w21 * w30 * w40) + &
                x5 * (w10 * w20 * w31 * w40) + &
                x6 * (w10 * w20 * w30 * w41) + &
                x7 * (w11 * w21 * w30 * w40) + &
                x8 * (w10 * w20 * w31 * w41) + &
                x9 * (w11 * w20 * w31 * w40) + &
                x10 * (w10 * w21 * w30 * w41) + &
                x11 * (w10 * w21 * w31 * w40) + &
                x12 * (w11 * w20 * w30 * w41) + &
                x13 * (w11 * w21 * w31 * w40) + &
                x14 * (w11 * w21 * w30 * w41) + &
                x15 * (w11 * w20 * w31 * w41) + &
                x16 * (w10 * w21 * w31 * w41)
        print*, x_inter
      end subroutine interpolate4D

subroutine binary_search(x,x1,i0)
 implicit none
 real(kind=8), intent(in) :: x(:)
 real(kind=8), intent(in) :: x1
 integer                  :: i,i1
 integer, intent(inout)   :: i0
 i0 = 1
 i1 = size(x)
 ! test if x1 is outside of interval [x(1), x(end)]
 if (x1 < x(i0)) then
   i = i0
   print*, "Outside range in binary search"
   return
 end if
 if (x1 > x(i1)) then
   i = i1
   print*, "Outside range in binary search"
   return
 end if

 ! binary search
 do while (i1 > (i0 + 1))
   i = (i1 + i0) / 2
   if (x(i) < x1) then
     i0 = i
   else if (x(i) > x1) then
     i1 = i
   else
     return
   end if
 end do
 !i = merge(i0, i1, ((2 * x1) < (x(i0) + x(i1))))
end subroutine binary_search

end module read_qse
