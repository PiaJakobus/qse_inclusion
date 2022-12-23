program main
!use iso_fortran_env, only: dp => real64
      use read_qse
      include 'implno.dek'
      include 'vector_eos.dek'
      include 'const.dek'

      integer, parameter :: ionmax = 21, nromax = 1
      integer            :: i,j

      real(dp)   :: temp, den, ye, x_cl, p
      real(dp)   :: pres, ener, enullm 
      real(dp)   :: eden_0, eden_int, enbyrst,&
                        eden_int_old, delta_eden_int
      real(dp)   :: abar,zbar  
      real(dp), dimension(ionmax) :: xmass,aion,zion
      !real(dp), dimension(21) ::  x

      enbyrst        = 0.0_dp
      delta_eden_int = 100.0_dp
      enullm         = 8.8_dp * 1.60218e-6_dp 

      den  = 1.3e8_dp
      temp = 4.95e9_dp
      ye   = 0.471
      x_cl = -1.6_dp 
 


      aion(1)   = 1.0 ; zion(1)  = 0.0
      aion(2)   = 1.0 ; zion(2)  = 1.0
      aion(3)   = 2.0 ; zion(3)  = 1.0
      aion(4)   = 3.0 ; zion(4)  = 2.0
      aion(5)   = 4.0 ; zion(5)  = 2.0
      aion(6)   = 2.0 ; zion(6)  = 6.0
      aion(7)   = 16.0 ; zion(7)  = 8.0
      aion(8)   = 20.0 ; zion(8)  = 10.0
      aion(9)   = 24.0 ; zion(9)  = 12.0
      aion(10)  = 28.0 ; zion(10)  = 14.0
      aion(11)  = 32.0 ; zion(11)  = 16.0
      aion(12)  = 48.0 ; zion(12)  = 18.0
      aion(13)  = 40.0 ; zion(13)  = 20.0
      aion(14)  = 44.0 ; zion(14)  = 22.0
      aion(15)  = 48.0 ; zion(15)  = 24.0
      aion(16)  = 52.0 ; zion(16)  = 26.0
      aion(17)  = 54.0 ; zion(17)  = 26.0
      aion(18)  = 56.0 ; zion(18)  = 28.0
      aion(19)  = 52.0 ; zion(19)  = 23.0
      aion(20)  = 62.0 ; zion(20)  = 28.0
      aion(21)  = 14.0 ; zion(21)  = 7.0

     
      ! set the input vector. pipeline is only 1 element long
      temp_row(1) = temp ; den_row(1)  = den ; abar_row(1) = 6 ; zbar_row(1) = 6
      jlo_eos = 1 ; jhi_eos = 1 ; ye_row(1) = ye 
            
      !call pretty_eos_out('eosfxt:  ')
      
      eden_0 = 1e17_dp 
      etot_row(1) = 5e16_dp 



      ! TODO: recalculate x_cl
      ! --------- check via eden ------

      !call qse(xmass_row(:,1),den_row(1),temp_row(1),ye_row(1),x_cl)
      j = 0
      do while ((delta_eden_int > 0.01) .and. (j < 10))
            call eosfxt
            call qse(xmass_row(:,1),den_row(1),temp_row(1),ye_row(1),x_cl)
            abar   = 1.0/sum(xmass_row(1:ionmax,1)/aion(1:ionmax))
            zbar   = abar * sum(xmass_row(1:ionmax,1) * zion(1:ionmax)/aion(1:ionmax))
            abar_row(1) = abar ; zbar_row(1) = zbar
            do i = 1, ionmax
              enbyrst = enbyrst + ((zion(i) * mp + (aion(i) - zion(i)) * mn) &
                        / (amu * aion(i))) * den_row(1) * xmass_row(i,1) 
            end do 
            enbyrst = 1.66e-24_dp * enbyrst + den_row(1) * (-mn + enullm + ye_row(1) * me)
            eden_int = eden_0 - enbyrst
            delta_eden_int = abs(etot_row(1) - eden_int)
            etot_row(1) = eden_int 
            j = j + 1
!            print*, enbyrst,delta_eden_int, etot_row(1), den_row(1),&
!                    temp_row(1),ye_row(1)
            enbyrst = 0.0_dp
      enddo


      end   



      include 'eosfxt.f90'

