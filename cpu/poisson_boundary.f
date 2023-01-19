      !****************************************************
      ! boundary condition for poisson, linear interpolation, fiducial
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      subroutine poisson_boundary(m)
        use mesh
        implicit none
        double precision, external :: gravity_profile

        ! input mesh level
        integer, intent(in) :: m

        ! variables
        integer i, k, l, ib, kb, lb

        ! set boundary conditions
        do l=1,N+2
         do k=1,N+2
           i = 1        ! left boundary condition
           if(m == 1) then    ! basic mesh
             gravity(i,k,l,m) = gravity_profile(i,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1

             if (mod(l,2).ne.0.and.mod(k,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb+1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(k,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb+1,lb-1,m-1)
             endif
             if (mod(l,2).ne.0.and.mod(k,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb-1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(k,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb-1,lb-1,m-1)
             endif

           endif

           i = N+2      ! right boundary condition
           if(m == 1) then    ! basic mesh
             gravity(i,k,l,m) = gravity_profile(i,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1

             if (mod(l,2).ne.0.and.mod(k,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb+1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(k,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb+1,lb-1,m-1)
             endif
             if (mod(l,2).ne.0.and.mod(k,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb-1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(k,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb-1,lb-1,m-1)
             endif

           endif

         enddo
        enddo

        do l=1,N+2
         do i=1,N+2
           k = 1        ! top boundary condition
           if(m == 1) then    ! basic mesh
             gravity(i,k,l,m) = gravity_profile(i,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1

             if (mod(l,2).ne.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb+1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb+1,lb-1,m-1)
             endif
             if (mod(l,2).ne.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb+1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb+1,lb-1,m-1)
             endif

           endif

           k = N+2      ! bottom boundary condition
           if(m == 1) then    ! basic mesh
             gravity(i,k,l,m) = gravity_profile(i,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1

             if (mod(l,2).ne.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb-1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb-1,lb-1,m-1)
             endif
             if (mod(l,2).ne.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb-1,lb+1,m-1)
             endif
             if (mod(l,2).eq.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb-1,lb-1,m-1)
             endif

           endif

         enddo
        enddo

        do k=1,N+2
         do i=1,N+2
           l = 1        ! down boundary condition
           if(m == 1) then    ! basic mesh
             gravity(i,k,l,m) = gravity_profile(i,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1

             if (mod(k,2).ne.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb+1,lb+1,m-1)
             endif
             if (mod(k,2).eq.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb-1,lb+1,m-1)
             endif
             if (mod(k,2).ne.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb+1,lb+1,m-1)
             endif
             if (mod(k,2).eq.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb-1,lb+1,m-1)
             endif

           endif

           l = N+2      ! up boundary condition
           if(m == 1) then    ! basic mesh
             gravity(i,k,l,m) = gravity_profile(i,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1

             if (mod(k,2).ne.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb+1,lb-1,m-1)
             endif
             if (mod(k,2).eq.0.and.mod(i,2).ne.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib+1,kb-1,lb-1,m-1)
             endif
             if (mod(k,2).ne.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb+1,lb-1,m-1)
             endif
             if (mod(k,2).eq.0.and.mod(i,2).eq.0) then
                 gravity(i,k,l,m) = 0.75d0*gravity(ib,kb,lb,m-1) +
     =                            0.25d0*gravity(ib-1,kb-1,lb-1,m-1)
             endif

           endif

         enddo
        enddo

      end subroutine poisson_boundary
      !**************************************************