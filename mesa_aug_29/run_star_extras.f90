! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine axion_other_energy(id, ierr)
         use star_def
         use auto_diff
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k

         integer :: i_na23
         real(dp) :: delta_E, g_P, g_N ! physical parameters
         real(dp) :: g ! intermediate values
         real(dp), pointer :: dm_na23(:), T(:), E(:) ! intermediate arrays
         real(dp) :: dm0, T0, N0, mu0 ! baseline values

         ierr = 0
         i_na23 = 31 ! index of species

         delta_E = 7.05e-7 ! ergs; energy of M1 transition; =440 keV
         g_P = 1.0e-9 ! unitless; proton coupling
         g_N = 0.0 ! unitless; neutron coupling

         dm0 = 9.942e31 ! grams; baseline; =0.05 solar masses
         T0 = 5.106e+9 ! Kelvins; baseline
         N0 = 1.21e48 ! axions/s; baseline
         mu0 = 1.5 ! unitless; chemical potential

         g = (g_P - 0.058*g_N) * 1e9
         dm_na23(:) = (s% dm(:)) * (s% xa_start(i_na23, :)) / dm0
         T(:) = T0 / (s% T(:))
         E(:) = delta_E * N0 * dm_na23(:) * g*g / ( exp(-T(:)) + 1.5 )

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extra_heat(:) = E(:)
         !s% extra_heat(:) = s% extra_power_source
         ! note that extra_heat is type(auto_diff_real_star_order1) so includes partials.
      end subroutine axion_other_energy

      subroutine extras_controls(id, ierr)
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! this is the place to set any procedure pointers you want to change
           ! e.g., other_wind, other_mixing, other_energy  (see star_data_procedures.inc)
           s% other_neu => axion_other_energy

        end subroutine extras_controls

      end module run_star_extras
      
