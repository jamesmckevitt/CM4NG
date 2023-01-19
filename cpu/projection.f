      !****************************************************
      ! projection from the m-grid to the (m-1)-grid
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      subroutine projection(m)
        use mesh
        implicit none 

        ! input mesh level
        integer, intent(in) :: m

        ! variables
        double precision dmdens, dmimpx, dmimpy, dmimpz
        integer i, k, l, ib, kb, lb

        do l=2,N+1,2
         do k=2,N+1,2
          do i=2,N+1,2

            ! calculate median density and momentum
            dmdens = ( density(i+1,k+1,l+1,m) +
     =                 density(i+1,k+1,l  ,m) +
     =                 density(i+1,k  ,l+1,m) +
     =                 density(i+1,k  ,l  ,m) +
     =                 density(i  ,k+1,l+1,m) +
     =                 density(i  ,k+1,l  ,m) +
     =                 density(i  ,k  ,l+1,m) +
     =                 density(i  ,k  ,l  ,m) ) / 8.d0

            ib = N/4 + i/2 + 1
            kb = N/4 + k/2 + 1
            lb = N/4 + l/2 + 1

            density(ib,kb,lb,m-1) = dmdens

          enddo
         enddo
        enddo

      end subroutine projection
      !**************************************************