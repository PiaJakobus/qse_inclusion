module read_qse
use hdf5
use hdf5_utils
!use iso_fortran_env, only: dp => real64
use linspace_mod


implicit none

contains

subroutine data_qse(y,filename)
  ! This subroutine reads params.jld
  ! add path to filename
  implicit none

  real(dp), intent(inout)      :: y(75,1,1,75,5371)
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
  real(dp), intent(inout) :: yrange(75)
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
      real(dp), intent(in)    :: rho, temp, ye, x_cl,enbyrst
      real(dp), intent(inout) :: x(5371)
      real(dp)                :: y(5371,75,1,1,75)
      character(len=60) :: path
      ! r2/t2: 1.2857142857142857e8 / 4.972972972972973e9
      ! r2/t3: 1.2857142857142857e8 / 4.945945945945947e9
      ! r3/t2: 1.4285714285714287e8 / 4.972972972972973e9
      ! r3/t3: 1.4285714285714287e8 / 4.945945945945947e9

      !call interpolate4D(4.e8_dp,3.1e9_dp,0.49_dp,-0.6_dp,2,x)
      call interpolate4D(1.3e8_dp,4.95e9_dp,0.471_dp,-1.6_dp,2,x)
      !call interpolate4D(1.e8_dp,5.e9_dp,0.5_dp,-1.5_dp,2,x)
      !call interpolate4D(1.e9_dp,3.e9_dp,0.47_dp,-4._dp,2,x)
end subroutine qse



subroutine find_index(rho,rrange,i_r,tem,trange,i_t,ye,yrange,i_ye,cl,srange,i_cl)
  implicit none
  real(dp), intent(in) :: yrange(75), rrange(64),trange(75),srange(75)
  real(dp), intent(in) :: rho,tem,ye,cl
  real(dp)           :: drho,dcl,dtem,dye
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
!  print*,'---'
!  print*, "rrange(i_r-1)",rrange(i_r-1)
!  print*, "rrange(i_r)",rrange(i_r)
!  print*, "rho",rho
!  print*, "i_r", i_r
!  print*,'---'
!
!  print*, "srange(i_cl-1)",srange(i_cl-1)
!  print*, "srange(i_cl)",srange(i_cl)
!  print*, "cl",cl
!  print*, "i_cl", i_cl
!  print*,'---'
!
!  print*, "trange(i_t-1)",trange(i_t-1)
!  print*, "trange(i_t)",trange(i_t)
!  print*, "tem",tem
!  print*, "i_t", i_t
!  print*,'---'
!
!  print*, "yrange(i_ye-1)",yrange(i_ye-1)
!  print*, "yrange(i_ye)",yrange(i_ye)
!  print*, "ye",ye
!  print*, "i_ye", i_ye
!  print*,'---'
end subroutine find_index

subroutine inject_to_file(f_red,rrange,trange,yrange,srange)
  ! Julia: yy,tt,rr,ss
  real(dp), allocatable, intent(inout) :: f_red(:,:,:,:,:)
  real(dp) :: fl(5371,75,1,1,75)
  real(dp),intent(in):: yrange(75), rrange(64),trange(75),srange(75)
  integer :: yy,ss,rr,tt
  character(len=100)         :: path
  character(len=10)       :: tmp1,tmp2
  allocate (f_red(21,75,62,64,75))
  do rr = 1,64
    do tt = 1,62
      write(tmp1,'(A1,i0)') 'r',rr-1
      write(tmp2,'(A2,i0)') '/t',tt
      print*, tmp1, tmp2
      path = "/c/pia/output/"//trim(tmp1)//trim(tmp2)//"/QSE_table.jld"
      call data_qse(fl,path)
      f_red(1:20,1:75,tt,rr,1:75) = fl((/1,2,3,5,6,7,8,31,84,145,215,292,376,465,560,659,661,763,518,769/),1:75,1,1,1:75)
      f_red(21,1:75,tt,rr,1:75) = 0.0
        do yy = 1,75
          do ss = 1,75
            f_red(:,yy,tt,rr,ss) = f_red(:,yy,tt,rr,ss) / sum(f_red(:,yy,tt,rr,ss))
          end do
        end do
    end do
    print*, path
  end do
  print*, "saving table..."
  open(unit=90,file="/c/pia/qse_inject.dat",status='replace',form='unformatted')
  write(90) trange,rrange,yrange, srange
  write(90) f_red
  close(90)
  print*, "saving complete!"
end subroutine inject_to_file

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
      real(dp), intent(in)  :: rho,tem,ye,cl
      integer, intent(in)   :: index_part
      real(dp), intent(out) :: x_inter(21)
      real(dp)              :: yrange(75), rrange(64),trange(75),srange(75)
     real(dp), allocatable  :: f_red(:,:,:,:,:)
      real(dp)              :: mass = 0., mass1 = 0.
      integer               :: i_r,i_t,i_ye,i_cl, j
      real(dp)              :: x0000,x1111,x1000,x0100,&
                                   x0010,x0001,x1100,x0011,&
                                   x1010,x0101,x0110,x1001,&
                                   x1110,x1101,x1011,x0111
      real(dp)              :: w10,w20,w30,w40,&
                                   w11,w21,w31,w41,&
                                   drho,dcl,dtem,dye
      character(len=60)     :: path1,path2,path3,path4,tmp1,tmp2

      rrange = linspace(1e8_dp,1e9_dp,64)
      trange = linspace(5e9_dp,3e9_dp,75)
      yrange = linspace(0.5_dp,0.47_dp,75)
      srange = linspace(-1.5_dp,-4._dp,75)
      drho   = rrange(2) - rrange(1)
      dcl    = abs(srange(2)) - abs(srange(1))
      dtem   = trange(1) - trange(2)
      dye    = yrange(1) - yrange(2)
      !call range_params(yrange, "/yrange")
      !call range_params(srange, "/srange")


      call find_index(rho,rrange,i_r,tem,trange,i_t,ye,yrange,i_ye,cl,srange,i_cl)


      !call inject_to_file(f_red,rrange,trange,yrange,srange)

      allocate (f_red(21,75,62,64,75))
      open(unit=95,file="/c/pia/qse_inject.dat",form='unformatted')
      read(95) trange,rrange,yrange,srange
      read(95) f_red
      close(95) 
      print*, "Does it sum to one?", sum(f_red(:,1,1,1,1))

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


!      ! ***** looping over entire array *****
      do j = 1, size(x_inter)
            ! Julia: fl__ -> ye, t, rho, s
            ! x____ -> rho,tem,ye,cl
            x0000 = f_red(j,i_ye-1,i_t-1,i_r-1,i_cl-1)
            x0010 = f_red(j,i_ye,i_t-1,i_r-1,i_cl-1)
            x0001 = f_red(j,i_ye-1,i_t-1,i_r-1,i_cl)
            x0011 = f_red(j,i_ye,i_t-1,i_r-1,i_cl)
            x0100 = f_red(j,i_ye-1,i_t,i_r-1,i_cl-1)
            x0101 = f_red(j,i_ye-1,i_t,i_r-1,i_cl)
            x0110 = f_red(j,i_ye,i_t,i_r-1,i_cl-1)
            x0111 = f_red(j,i_ye,i_t,i_r-1,i_cl)
            x1000 = f_red(j,i_ye-1,i_t-1,i_r,i_cl-1)
            x1010 = f_red(j,i_ye,i_t-1,i_r,i_cl-1)
            x1001 = f_red(j,i_ye-1,i_t-1,i_r,i_cl)
            x1011 = f_red(j,i_ye,i_t-1,i_r,i_cl)
            x1111 = f_red(j,i_ye,i_t,i_r,i_cl)
            x1100 = f_red(j,i_ye-1,i_t,i_r,i_cl-1)
            x1110 = f_red(j,i_ye,i_t,i_r,i_cl-1)
            x1101 = f_red(j,i_ye-1,i_t,i_r,i_cl)


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
        if (x_inter(j) .ne. x_inter(j)) then
          x_inter(j) = 0
        else
          mass = mass + x_inter(j)
        endif
      end do
      deallocate(f_red)
      

      if ((mass - 1.0) > 1e-3) then 
        print*,"ERROR: mass not summing to 1", mass
        stop 
      end if 
      stop 
      end subroutine interpolate4D



end module read_qse
