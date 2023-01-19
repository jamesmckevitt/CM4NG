      !**************************************************
      ! gravitational potential from the analytic solution
      ! used to calculate the boundary potential at the coarsest grid
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !**************************************************
      double precision function gravity_profile(i,k,l,m)
        use machine
        use mesh
        implicit none

        double precision, external :: gravity_profile_an

        ! input index of cell
        integer, intent(in) :: i, k, l, m

        gravity_profile=gravity_profile_an(i,k,l,m)

      end function gravity_profile
      !**************************************************