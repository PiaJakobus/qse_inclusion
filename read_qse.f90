module read_qse
use hdf5
use hdf5_utils
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

subroutine range_params(yrange)
  ! This subroutine reads params.jld
  ! add path to filename
  real(kind=8), intent(inout) :: yrange(75)
  integer, parameter          :: ydims = 1
  integer(HID_T)              :: yshape(ydims)
  integer                     :: errorflag
  integer(HID_T)              :: file_id,group_id,root_id, d_id
  character(len=60)           :: dataname_y, filename, rootname

  yshape = shape(yrange)
  filename = "params.jld"
  rootname = "/"
  dataname_y = "/"//"yrange"
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
      real, intent(in)            :: rho, temp, ye, x_cl,enbyrst
      real(kind=8), intent(inout) :: x(5371)
      real(kind=8)                :: y(75,1,1,75,5371)
      real(kind=8)                :: yrange(75)
      call range_params(yrange)
      call data_qse(y)

end subroutine qse



end module read_qse
