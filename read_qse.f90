module read_qse
use hdf5 
use hdf5_utils
implicit none

contains 

subroutine qse(rho,temp,ye,x,x_cl,enbyrst)
      implicit none
      real, intent(in)    :: rho, temp, ye, x_cl 
      real, intent(inout) :: enbyrst
      
      real(kind=8), intent(inout) :: x(75)
      real(kind=8)           :: y(75,1,1,75,5371)
      INTEGER, PARAMETER :: ndims = 1
      INTEGER, PARAMETER :: ydims = 5
      INTEGER(HID_T)   :: xshape(ndims)
      INTEGER(HID_T)   :: yshape(ydims)
      
      integer            :: errorflag
      integer(HID_T)     :: file_id,group_id,root_id, d_id
      character(len=60)  :: dataname, filename, rootname 
       
      xshape = shape(x)
      filename = "params.jld"
      rootname = "/"
      dataname = "/"//"yrange"
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
      call H5DREAD_F(d_id, H5T_NATIVE_DOUBLE, x, xshape, errorflag)
      if (errorflag .lt. 0) then 
              print*, "Error data open"
              return
      end if 
 

      call h5dclose_f(d_id, errorflag)
      CALL h5gclose_f(root_id,errorflag)
      CALL h5fclose_f(file_id,errorflag)
      CALL h5close_f(errorflag)
      if (errorflag .lt. 0) then 
              print*, "Error closing files"
              return
      end if 
      yshape = shape(y)
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
      CALL h5gclose_f(root_id,errorflag)
      CALL h5fclose_f(file_id,errorflag)
      CALL h5close_f(errorflag)
      if (errorflag .lt. 0) then 
              print*, "Error closing files"
              return
      end if 
      print*, y

end subroutine qse

      

end module read_qse 
