      !****************************************************
      ! module for basic parameters and variables
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      module mesh
        implicit none

        !**************************************************
        ! mesh constants & variables
        !**************************************************
        double precision, parameter :: length = 4.5d0         ! length of computational domain
        double precision center                               ! center of computational domain
        integer, parameter :: N = 128                         ! number of cells on each nested grid
        integer, parameter :: Ndepth = 4                      ! depth of the nested mesh
        double precision h                                    ! size of the cell on the coarsest grid

        !**************************************************
        ! initial density setup
        !**************************************************

        double precision :: rho_sp=1.d0                       ! density of the sphere or ellipsoid

        ! homogeneous oblate ellipsoid (centered at the coordinate center)
        logical, parameter :: ellipsoid=.true.
        double precision :: rsp1=1.5d0                        ! major semi-axis of the ellipsoid
        double precision :: rsp3=0.75d0                       ! minor semi-axis of the ellipsoid

        ! homogeneous sphere (centered in the coordinate center)
        logical, parameter :: sphere=.false.
        double precision :: rsp=2.0d0                         ! radius of the sphere

        ! wide-separation binary
        logical, parameter :: binary=.false.
        double precision :: rsp_b1=0.2d0                      ! radius of the primary
        double precision :: rsp_b2=0.2d0                      ! radius of the cecondary
        double precision :: xb1=-0.5d0                        ! position of the primary center along the x-coordinate (y,z=0)
        double precision :: xb2=0.5d0                         ! position of the secondary center along the x-coordinate (y,z=0)
        double precision :: rho_sp1=2.d0                      ! density of the primary
        double precision :: rho_sp2=1.d0                      ! density of the secondary

        ! density distribution after Wang & Yen 2020, ApJS, 247, 2
        logical, parameter :: wang=.false.
        double precision :: rho_wang = 1.d0
        double precision :: r_wang = 0.25d0

        ! density distribution for a tight binary (dipole)
        logical, parameter :: dipole=.false.
        double precision :: rho_dp1=1.d0                      ! density of the secondary
        double precision :: rho_dp2=2.d0                      ! density of the primary

        !--------------------------------------------------
        ! physical constants 
        !--------------------------------------------------
        double precision, parameter :: GravConst=1.d0         ! Gravitational constant

        !**************************************************
        ! conservative and physics variables on mesh
        !**************************************************
        double precision density(N+2, N+2, N+2, Ndepth)       ! conservative density
        double precision densityCut(N+2, N+2, N+2, Ndepth)    !  density in the cut grid
        double precision gravity(N+2, N+2, N+2, Ndepth)       ! gravity
        double precision gravityCut(N+2, N+2, N+2, Ndepth)    ! gravity with the nested grid cut out
        double precision gravityCutInt(N+2, N+2, N+2, Ndepth) ! gravity with the nested grid cut out
        double precision gravityRed(N/2+2,N/2+2,N/2+2)        ! supplimentary file for averaged potential
        double precision densityDBL(N*2+2,N*2+2,N*2+2)        ! density on a doubled domain
        double precision gravityDBL(N*2+2,N*2+2,N*2+2)        ! gravity on a doubled domain
        !**************************************************

        !**************************************************
        ! OPENMP number of threads
        !**************************************************
        integer, parameter :: Nthreads=32

      end module mesh