program main
!use iso_fortran_env, only: dp => real64
      use read_qse
      include 'implno.dek'
      include 'vector_eos.dek'

      integer, parameter :: ionmax = 20 
      integer            :: i,j

      real(dp)   :: temp, den, ye, x_cl, p
      real(dp)   :: eden_0, eden_int, enbyrst,&
                        eden_int_old, delta_eden_int
      real(dp)   :: abar,zbar  
      real(dp), dimension(ionmax) :: xmass,aion,zion
      real(dp), dimension(5371) ::  x

      enbyrst        = 1.0_dp
      temp           = 1e9_dp
      den            = 1e9_dp
      delta_eden_int = 100.0_dp
      eden_int       = eden_0 - enbyrst
      eden_int_old   = eden_int

      call qse(den,temp,ye,x,x_cl,enbyrst)
stop

      ! ionmax  = number of isotopes in the network
      ! xmass   = mass fraction of isotope i
      ! aion    = number of nucleons in isotope i
      ! zion    = number of protons in isotope i

      ! set the mass fractions, z's and a's of the composition
      ! hydrogen, heliu, and carbon
      xmass(1) = 0.75 ; aion(1)  = 1.0  ; zion(1)  = 1.0
      xmass(2) = 0.23 ; aion(2)  = 4.0  ; zion(2)  = 2.0
      xmass(3) = 0.02 ; aion(3)  = 12.0 ; zion(3)  = 6.0

      ! average atomic weight and charge
      abar   = 1.0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

      ! set the input vector. pipeline is only 1 element long
      temp_row(1) = 1.e8 ; den_row(1)  = 1.e6 ; abar_row(1) = abar ; zbar_row(1) = zbar
      jlo_eos = 1 ; jhi_eos = 1

      !call eosfxt

      !call pretty_eos_out('eosfxt:  ')

      ! what if T falls behind T_qse in iteration?
      ! kepdata.py in python/source (kepler) approx21
      ! 1) mapping all nuclei -> 21 nuclei

      ! --------- check via eden --------
      j = 1
      do while ((delta_eden_int > 0.01) .and. (j < 10))
            !call eos(den,eden_int,ye,p,temp)
            eden_int = eden_0 - enbyrst
            delta_eden_int = abs(eden_int - eden_int_old)
            eden_int_old = eden_int
            j = j + 1
      enddo


      end   



      include 'eosfxt.f90'

