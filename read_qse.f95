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
      !call interpolate4D(1.3e8_dp,4.1e9_dp,0.471_dp,-1.6_dp,2,x)
      call interpolate4D(1.e8_dp,5.e9_dp,0.5_dp,-1.5_dp,2,x)
      !call interpolate4D(1.e9_dp,3.e9_dp,0.47_dp,-4._dp,2,x)
end subroutine qse



subroutine find_index(rho,rrange,i_r,tem,trange,i_t,ye,yrange,i_ye,cl,srange,i_cl)
  implicit none
  real(kind=8), intent(in) :: yrange(75), rrange(64),trange(75),srange(75)
  real(kind=8), intent(in) :: rho,tem,ye,cl
  real(kind=8)           :: drho,dcl,dtem,dye
  integer, intent(inout) :: i_r,i_t,i_ye,i_cl
  drho   = rrange(2) - rrange(1)
  dcl    = abs(srange(2)) - abs(srange(1))
  dtem   = trange(1) - trange(2)
  dye    = yrange(1) - yrange(2)

  i_r  = floor((rho - rrange(1)) / drho + 1.) + 1
  i_cl = floor((abs(cl) -  abs(srange(1))) / dcl + 1.) + 1
  i_t =  floor((trange(1) - tem) / dtem + 1.) + 1
  i_ye = floor((yrange(1) - ye) / dye + 1.) + 1

  if (rho == rrange(1)) then
      i_r = 1
  end if
  if (cl == srange(1)) then
      i_cl = 1
  end if
  if (tem == trange(1)) then
      i_t = 1
  end if
  if (ye == yrange(1)) then
      i_ye = 1
  end if
  if (rho == rrange(size(rrange))) then
      i_r = 64
  end if
  if (cl == srange(size(srange))) then
      i_cl = 75
  end if
  if (tem == trange(size(trange))) then
      i_t = 75
  end if
  if (ye == yrange(size(yrange))) then
      i_ye = 75
  end if
  print*,'---'
  print*, "rrange(i_r-1)",rrange(i_r-1)
  print*, "rrange(i_r)",rrange(i_r)
  print*, "rho",rho
  print*, "i_r", i_r
  print*,'---'

  print*, "srange(i_cl-1)",srange(i_cl-1)
  print*, "srange(i_cl)",srange(i_cl)
  print*, "cl",cl
  print*, "i_cl", i_cl
  print*,'---'

  print*, "trange(i_t-1)",trange(i_t-1)
  print*, "trange(i_t)",trange(i_t)
  print*, "tem",tem
  print*, "i_t", i_t
  print*,'---'

  print*, "yrange(i_ye-1)",yrange(i_ye-1)
  print*, "yrange(i_ye)",yrange(i_ye)
  print*, "ye",ye
  print*, "i_ye", i_ye
  print*,'---'
  stop
end subroutine find_index

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
      integer, intent(in)       :: index_part
      real(kind=8), intent(out) :: x_inter(5371)
      real(kind=8)              :: yrange(75), rrange(64),trange(75),srange(75)
      real(kind=8)              :: fl11(5371,75,1,1,75)
      real(kind=8)              :: fl12(5371,75,1,1,75)
      real(kind=8)              :: fl21(5371,75,1,1,75)
      real(kind=8)              :: fl22(5371,75,1,1,75)
      real(kind=8)              :: mass = 0., mass1 = 0.
      integer                   :: i_r,i_t,i_ye,i_cl, j
      real(kind=8)              :: x0000,x1111,x1000,x0100,&
                                   x0010,x0001,x1100,x0011,&
                                   x1010,x0101,x0110,x1001,&
                                   x1110,x1101,x1011,x0111
      real(kind=8)              :: w10,w20,w30,w40,&
                                   w11,w21,w31,w41,&
                                   drho,dcl,dtem,dye
      character(len=60)         :: path1,path2,path3,path4,tmp1,tmp2

      rrange = linspace(1e8_dp,1e9_dp,64)
      trange = linspace(5e9_dp,3e9_dp,75)
      yrange = linspace(0.5_dp,0.47_dp,75)
      srange = linspace(-1.5_dp,-4._dp,75)
      drho   = rrange(2) - rrange(1)
      dcl    = abs(srange(2)) - abs(srange(1))
      dtem   = trange(1) - trange(2)
      dye    = yrange(1) - yrange(2)
      ! I could also call the parameter space direclty. But maybe not in this routine:
      !call range_params(yrange, "/yrange")
      !call range_params(srange, "/srange")


      ! i_r = 1
      ! i_t = 1
      ! i_cl = 10
      ! i_ye = 10
      call find_index(rho,rrange,i_r,tem,trange,i_t,ye,yrange,i_ye,cl,srange,i_cl)

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
      ! https://en.wikipedia.org/wiki/Linear_interpolation
      w10 = (rrange(i_r) - rho) / drho
      w11 = 1. - w10
      w40 = (srange(i_cl) - cl) / dcl
      w41 = 1. - w40
      w20 = (- trange(i_t) + tem) / dtem
      w21 = 1. - w20
      w30 = (- yrange(i_ye) + ye) / dye
      w31 = 1. - w30
      print*, w10,w20,w30,w40,w11,w21,w31,w41

      ! ***** checking for corner cases *****
      if ((i_ye == 1) .and. (i_cl == 1)) then
        x0011 = fl11(j,i_ye,1,1,i_cl)
        x0111 = fl12(j,i_ye,1,1,i_cl)
        x1011 = fl21(j,i_ye,1,1,i_cl)
        x1111 = fl22(j,i_ye,1,1,i_cl)
        x_inter(j) = w10 * (x0011 * w20 + x0111 * w21) + &
                     w11 * (x1011 * w20 + x1111 * w21)
        return
      else if (i_ye == 1) then
        x0010 = fl11(j,i_ye,1,1,i_cl-1)
        x0011 = fl11(j,i_ye,1,1,i_cl)
        x0110 = fl12(j,i_ye,1,1,i_cl-1)
        x0111 = fl12(j,i_ye,1,1,i_cl)
        x1010 = fl21(j,i_ye,1,1,i_cl-1)
        x1011 = fl21(j,i_ye,1,1,i_cl)
        x1111 = fl22(j,i_ye,1,1,i_cl)
        x1110 = fl22(j,i_ye,1,1,i_cl-1)
        x_inter(j) = w10 * ( &
                     x0010 * (w20 * w40) + &
                     x0011 * (w20 * w41) + &
                     x0110 * (w21 * w40) + &
                     x0111 * (w21 * w41)) + &
                     w11 * ( &
                     x1010 * (w20 * w40) + &
                     x1011 * (w20 * w41) + &
                     x1111 * (w21 * w41) + &
                     x1110 * (w21 * w40))
        return
      else if (i_cl == 1) then
        x0001 = fl11(j,i_ye-1,1,1,i_cl)
        x0011 = fl11(j,i_ye,1,1,i_cl)
        x0101 = fl12(j,i_ye-1,1,1,i_cl)
        x0111 = fl12(j,i_ye,1,1,i_cl)
        x1001 = fl21(j,i_ye-1,1,1,i_cl)
        x1011 = fl21(j,i_ye,1,1,i_cl)
        x1111 = fl22(j,i_ye,1,1,i_cl)
        x1101 = fl22(j,i_ye-1,1,1,i_cl)
        x_inter(j) = w10 * ( &
                     x0001 * (w20 * w30) + &
                     x0011 * (w20 * w31) + &
                     x0101 * (w21 * w30) + &
                     x0111 * (w21 * w31)) + &
                     w11 * ( &
                     x1001 * (w20 * w30) + &
                     x1011 * (w20 * w31) + &
                     x1111 * (w21 * w31) + &
                     x1101 * (w21 * w30))
        return
      end if

      ! ***** looping over entire array *****
      do j = 1, size(x_inter)
            ! Julia: fl__ -> ye, t, rho, s
            ! x____ -> rho,tem,ye,cl
            x0000 = fl11(j,i_ye-1,1,1,i_cl-1)
            x0010 = fl11(j,i_ye,1,1,i_cl-1)
            x0001 = fl11(j,i_ye-1,1,1,i_cl)
            x0011 = fl11(j,i_ye,1,1,i_cl)
            x0100 = fl12(j,i_ye-1,1,1,i_cl-1)
            x0101 = fl12(j,i_ye-1,1,1,i_cl)
            x0110 = fl12(j,i_ye,1,1,i_cl-1)
            x0111 = fl12(j,i_ye,1,1,i_cl)
            x1000 = fl21(j,i_ye-1,1,1,i_cl-1)
            x1010 = fl21(j,i_ye,1,1,i_cl-1)
            x1001 = fl21(j,i_ye-1,1,1,i_cl)
            x1011 = fl21(j,i_ye,1,1,i_cl)
            x1111 = fl22(j,i_ye,1,1,i_cl)
            x1100 = fl22(j,i_ye-1,1,1,i_cl-1)
            x1110 = fl22(j,i_ye,1,1,i_cl-1)
            x1101 = fl22(j,i_ye-1,1,1,i_cl)


! =================================================
! TODO: i_cl, i_cl not working because reverse: DONE!
! TODO: path directories einpflegen: DONE!
! TODO: Loop fuer alle elements adden : DONE!
! TODO: fix binary search -> DONE
! TODO: run on harddrive -> DONE
! TODO: write tests -> DONE
! TODO: write if clauses for edges
! TODO: one_zone.f95
! =================================================


      x_inter(j) = w10 * ( &
                   x0000 * (w20 * w30 * w40) + &
                   x0010 * (w20 * w31 * w40) + &
                   x0001 * (w20 * w30 * w41) + &
                   x0011 * (w20 * w31 * w41) + &
                   x0100 * (w21 * w30 * w40) + &
                   x0101 * (w21 * w30 * w41) + &
                   x0110 * (w21 * w31 * w40) + &
                   x0111 * (w21 * w31 * w41)) + &
                   w11 * ( &
                   x1000 * (w20 * w30 * w40) + &
                   x1010 * (w20 * w31 * w40) + &
                   x1001 * (w20 * w30 * w41) + &
                   x1011 * (w20 * w31 * w41) + &
                   x1111 * (w21 * w31 * w41) + &
                   x1100 * (w21 * w30 * w40) + &
                   x1110 * (w21 * w31 * w40) + &
                   x1101 * (w21 * w30 * w41))
            if (x_inter(j) > 1.e-6_dp) then
              print*, 100._dp*abs(1.0-x1111/x_inter(j)),x1111,x_inter(j)
            end if
        if (x_inter(j) .ne. x_inter(j)) then
          x_inter(j) = 0
        else
          mass = mass + x_inter(j)
          mass1 = mass1 + x1111
        endif
      end do
        print*,"--------------", mass, mass1, abs(1._dp-mass/mass1)

      end subroutine interpolate4D



end module read_qse
