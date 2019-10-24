      PROGRAM MAIN

      IMPLICIT NONE

      integer wavelength_tab_size
      PARAMETER(wavelength_tab_size=4)
      double precision wavelength_tab(wavelength_tab_size)

      wavelength_tab(1)=3.00D-1
      wavelength_tab(2)=4.00D-1
      wavelength_tab(3)=6.00D-1
      wavelength_tab(4)=9.99D-1

      call compute_tab(wavelength_tab,wavelength_tab_size)

      END
