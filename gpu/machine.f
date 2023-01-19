      !****************************************************
      ! module for machine constants
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      module machine
        implicit none

        ! PI value
        double precision, parameter :: pi = 3.141592653589793238462d0

        ! EPSILON value
        double precision, parameter :: eps = 1.0d-8

      end module machine