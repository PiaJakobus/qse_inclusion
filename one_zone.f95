program one_zone

use read_qse

implicit none


integer, parameter :: n = 50 ! guesses for T
integer :: i,j
real(kind=8)   :: ye,den, p, temp, eden, eden_int, enbyrst, x_cl
real(kind=8), dimension(5371) ::  x
real(kind=8)   :: delta_eden, eps, eden_int_old, eden_0

eps = 0.1
enbyrst = 1.
temp = 10.

call hydro(den,eden,ye,x_cl)
call eos(den,eden-enbyrst,ye,p,temp)


! what if T falls behind T_qse in iteration?
! kepdata.py in python/source (kepler) approx21
! 1) mapping all nuclei -> 21 nuclei

! --------- check via eden --------
delta_eden = 100.
j = 1
eden_int = eden_0 - enbyrst
eden_int_old = eden_int

do while ((delta_eden > 0.01) .and. (j < 10))
      call eos(den,eden_int,ye,p,temp)
      call qse(den,temp,ye,x,x_cl,enbyrst)
      eden_int = eden_0 - enbyrst
      delta_eden = abs(eden_int - eden_int_old)
      eden_int_old = eden_int
      j = j + 1
enddo



! --------------------------------------------------------------
! --------------------------------------------------------------

contains
subroutine eos(rho,eden,ye,p,T)
      implicit none
      real(kind=8), intent(in) :: rho, eden, ye
      real(kind=8), intent(inout) :: T
      real(kind=8), intent(out) :: p
      T = 1.01 * T
      p = 100. + T

end subroutine eos

subroutine hydro(den,eden,ye,x_cl)
      implicit none
      real(kind=8), intent(inout) :: den, eden,ye,x_cl
      den     = 10.**i + 10.
      eden    = 10.**i + 10.
      ye      = i + 1
end subroutine hydro


end program one_zone
