      !****************************************************
      ! load physics problem subroutine
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      subroutine load
        use mesh
        implicit none
        double precision, external :: density_profile

      !*************************************************
      !  seting the initial density distribution. See density_init_profile.f
      !*************************************************

        ! variables
        double precision densval
        integer i, k, l, m
        logical is_empty

        ! calculate the size of the coarsest cell
        h = length/N

        ! calculate the center of the computational domain
        center = length/2.d0

        ! set initial profiles
        do m=1,Ndepth
         do l=2,N+1
          do k=2,N+1
           do i=2,N+1

            ! initialize density
            densval = density_profile(i,k,l,m)
            density(i,k,l,m) = densval

            is_empty = (m < Ndepth) .and. (i > N/4+1)
     =                  .and. (i < 3*N/4+2)
     =                  .and. (k > N/4+1)
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l > N/4+1)
     =                  .and. (l < 3*N/4+2)

            If (.not.is_empty) then
               densityCut(i,k,l,m)=densval
            else
               densityCut(i,k,l,m)=0.d0
            endif

            ! initialize gravitational potential
            gravity(i,k,l,m) = 0.d0

           enddo
          enddo
         enddo
        enddo

      end subroutine load
      !**************************************************