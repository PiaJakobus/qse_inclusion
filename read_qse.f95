module read_qse
use hdf5
use hdf5_utils
use iso_fortran_env, only: dp => real64
use linspace_mod


implicit none

contains

subroutine data_qse(y,filename)
  ! This subroutine reads params.jld
  ! add path to filename
  implicit none

  real(kind=8), intent(inout)      :: y(75,1,1,75,5371)
  character(len=60), intent(in)    :: filename
  integer, parameter               :: ydims = 5
  integer(HID_T)                   :: yshape(ydims)
  integer                          :: errorflag
  integer(HID_T)                   :: file_id,group_id,root_id, d_id
  character(len=60)                :: dataname, rootname

  !filename = path!"QSE_table.jld"
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
  !print*, y(5,1,1,5,:)
  !stop
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
      real(kind=8)                :: y(5371,75,1,1,75)
      character(len=60) :: path
      call interpolate4D(1.3e8_dp,4.95e9_dp,0.471_dp,-1.6_dp,2,x)
end subroutine qse

! r2/t2: 1.2857142857142857e8 / 4.972972972972973e9
! r2/t3: 1.2857142857142857e8 / 4.945945945945947e9
! r3/t2: 1.4285714285714287e8 / 4.972972972972973e9
! r3/t3: 1.4285714285714287e8 / 4.945945945945947e9
subroutine interpolate4D(rho,tem,ye,cl,index_part,x_inter)
  ! f11: rho(i_r),T(i_t)
  ! f12: rho(i_r), T(i_t+1)
  ! f21: rho(i_r+1),T(i_t)
  ! f22: rho(i_r+1),T(i_t+1)
  !
  ! This subroutine interpolates a 4 dimensional table to a single value
  ! which is the abundance of a particle i. Subroutine loops over particles
  !stored in x_inter. Therefore the return value x_inter has
  ! length 5371.
      implicit none
      real(kind=8), intent(in)  :: rho,tem,ye,cl
      integer, intent(in)       ::  index_part
      real(kind=8), intent(out) :: x_inter(5371)
      real(kind=8)              :: yrange(75), rrange(64),trange(75),srange(75)
      real(kind=8)              :: fl11(5371,75,1,1,75)
      real(kind=8)              :: fl12(5371,75,1,1,75)
      real(kind=8)              :: fl21(5371,75,1,1,75)
      real(kind=8)              :: fl22(5371,75,1,1,75)
      real(kind=8)              :: mass = 0., mass1 = 0.
      integer                   :: i_r,r,i_t,i_ye,i_cl, j
      real(kind=8)              :: x0000,x1111,x1000,x0100,&
                                   x0010,x0001,x1100,x0011,&
                                   x1010,x0101,x0110,x1001,&
                                   x1110,x1101,x1011,x0111
      real(kind=8)              :: w10,w20,w30,w40,&
                                   w11,w21,w31,w41
      character(len=60)         :: path1,path2,path3,path4,tmp1,tmp2

      rrange = linspace(1e8_dp,1e9_dp,64)
      trange = linspace(5e9_dp,3e9_dp,75)
      yrange = linspace(0.5_dp,0.47_dp,75)
      srange = linspace(-1.5_dp,-4._dp,75)
      ! I could also call the parameter space direclty. But maybe not in this routine:
      !call range_params(yrange, "/yrange")
      !call range_params(srange, "/srange")
      i_r = MINLOC(rrange, dim=1,mask=(rho < rrange))
      print*,'---'
      print*, "rrange(i_r-1)",rrange(i_r-1)
      print*, "rrange(i_r)",rrange(i_r)
      print*, "rho",rho
      print*, "i_r", i_r
      print*,'---'
      i_cl = MINLOC(abs(srange), dim=1,mask=(abs(cl) < abs(srange)))
      print*, "srange(i_cl-1)",srange(i_cl-1)
      print*, "srange(i_cl)",srange(i_cl)
      print*, "cl",cl
      print*, "i_cl", i_cl
      print*,'---'
      i_t = MINLOC(trange, dim=1,mask=(tem < trange))
      print*, "trange(i_t+1)",trange(i_t+1)
      print*, "trange(i_t)",trange(i_t)
      print*, "tem",tem
      print*, "i_t", i_t
      print*,'---'
      i_ye = MINLOC(yrange, dim=1,mask=(ye < yrange))
      print*, "yrange(i_ye+1)",yrange(i_ye+1)
      print*, "yrange(i_ye)",yrange(i_ye)
      print*, "ye",ye
      print*, "i_ye", i_ye
      print*,'---'

      ! i_r = 1
      ! i_t = 1
      ! i_cl = 10
      ! i_ye = 10
      i_r = i_r - 1 ! directories for rho start at zero
      write(tmp1,'(a,i0,a)') 'r',i_r-1
      write(tmp2,'(a,i0,a)') '/t',i_t
      path1 = "qse_table/"//trim(tmp1) // trim(tmp2) // "/QSE_table.jld"
      write(tmp1,'(a,i0,a)') 'r',i_r-1
      write(tmp2,'(a,i0,a)') '/t',i_t+1
      path2 = "qse_table/"//trim(tmp1) // trim(tmp2) // "/QSE_table.jld"
      write(tmp1,'(a,i0,a)') 'r',i_r
      write(tmp2,'(a,i0,a)') '/t',i_t
      path3 = "qse_table/"//trim(tmp1) // trim(tmp2) // "/QSE_table.jld"
      write(tmp1,'(a,i0,a)') 'r',i_r
      write(tmp2,'(a,i0,a)') '/t',i_t+1
      path4 = "qse_table/"//trim(tmp1) // trim(tmp2) // "/QSE_table.jld"
      i_r = i_r + 1

      call data_qse(fl11,path1)
      call data_qse(fl12,path2)
      call data_qse(fl21,path3)
      call data_qse(fl22,path4)

      ! ***** calculating weights *****
      ! see example in  https://en.wikipedia.org/wiki/Linear_interpolation
      ! rrange,srange have ascening order: rrange(i_r-1) < rho < rrange(i_r)
      ! w10: dimension 1 at index i
      ! w11: dimension 1 at index i+1
      w10 = (rrange(i_r) - rho) / (rrange(i_r) - rrange(i_r-1))
      w11 = 1. - w10
      ! w40: dimension 4 at index i
      ! w41: dimension 4 at index i+1
      w40 = (srange(i_cl) - cl) / (srange(i_cl) - srange(i_cl-1))
      w41 = 1. - w40
      !
      ! trange,yrange have descending order: trange(i_t) > tem > trange(i_t+1)
      ! w21: dimension 2 at index i (the sign doesn't matter, only ratio)
      ! w20: dimension 2 at index i+1
      w21 = (trange(i_t) - tem) / (trange(i_t) - trange(i_t+1))
      w20 = 1. - w21
      ! w30: dimension 2 at index i (the sign doesn't matter, only ratio)
      ! w31: dimension 2 at index i+1
      w31 = (yrange(i_ye) - ye) / (yrange(i_ye) - yrange(i_ye+1))
      w30 = 1. - w31
      print*, w10,w20,w30,w40,w11,w21,w31,w41
      ! **********


      do j = 1, size(x_inter)
        !(size(yrange, 1), size(trange, 1), size(rrange, 1), size(srange, 1))


            x0000 = fl11(j,i_ye,1,1,i_cl)
            x0010 = fl11(j,i_ye+1,1,1,i_cl)
            x0001 = fl11(j,i_ye,1,1,i_cl+1)
            x0011 = fl11(j,i_ye+1,1,1,i_cl+1)

            x0100 = fl12(j,i_ye,1,1,i_cl)
            x0101 = fl12(j,i_ye,1,1,i_cl+1)
            x0110 = fl12(j,i_ye+1,1,1,i_cl)
            x0111 = fl12(j,i_ye+1,1,1,i_cl+1)

            x1000 = fl21(j,i_ye,1,1,i_cl)
            x1010 = fl21(j,i_ye+1,1,1,i_cl)
            x1001 = fl21(j,i_ye,1,1,i_cl+1)
            x1011 = fl21(j,i_ye+1,1,1,i_cl+1)

            x1111 = fl22(j,i_ye+1,1,1,i_cl+1)
            x1100 = fl22(j,i_ye,1,1,i_cl)
            x1110 = fl22(j,i_ye+1,1,1,i_cl)
            x1101 = fl22(j,i_ye,1,1,i_cl+1)


! =================================================
! TODO: i_cl, i_cl not working because reverse: DONE!
! TODO: path directories einpflegen: DONE!
! TODO: Loop fuer alle elements adden : DONE!
! TODO: fix binary search -> DONE
! TODO: run on harddrive -> DONE
! TODO: write tests -> DONE
! TODO: one_zone.f95
! =================================================


      ! each direction is multiplied with weights
      ! Example for 1D interpolation:
      !         x0000: value of x,y,z,h at (i_x,i_y,i_z,i_h) - 1
      !       * w10:   (x_i - x) / (x_i - x_i-1)
      !       + w11    (x - x_i-1) / (x_i - x_i-1)
      !       * x1000: value of x,y,z,h at i_x,(i_y,i_z,i_h) - 1)
      x_inter(j) = x0000 * (w10 * w20 * w30 * w40) + &
                   x0010 * (w10 * w20 * w31 * w40) + &
                   x0001 * (w10 * w20 * w30 * w41) + &
                   x0011 * (w10 * w20 * w31 * w41) + &
                   x0100 * (w10 * w21 * w30 * w40) + &
                   x0101 * (w10 * w21 * w30 * w41) + &
                   x0110 * (w10 * w21 * w31 * w40) + &
                   x0111 * (w10 * w21 * w31 * w41) + &
                   x1000 * (w11 * w20 * w30 * w40) + &
                   x1010 * (w11 * w20 * w31 * w40) + &
                   x1001 * (w11 * w20 * w30 * w41) + &
                   x1011 * (w11 * w20 * w31 * w41) + &
                   x1111 * (w11 * w21 * w31 * w41) + &
                   x1100 * (w11 * w21 * w30 * w40) + &
                   x1110 * (w11 * w21 * w31 * w40) + &
                   x1101 * (w11 * w21 * w30 * w41)
            print*,'x1111    ',x1111, "x_inter  ",x_inter(j)

        if (x_inter(j) .ne. x_inter(j)) then
          x_inter(j) = 0
        else
          mass = mass + x_inter(j)
          mass1 = mass1 + x1111
        endif
        !print*, x_inter(j), x0000, x0100, x1000, x1111
      end do
        print*,"--------------", mass, mass1, i_ye, i_cl,path4

      end subroutine interpolate4D



end module read_qse
