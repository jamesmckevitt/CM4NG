      !****************************************************
      ! saving the results
      ! only XY and XZ cuts are provided
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      subroutine save(timer)
        use machine
        use mesh 
        implicit none

        ! input timer
        double precision, intent(in) :: timer
        double precision, external :: gravity_profile_an
        double precision, external :: gradient_profile_an

        ! variables
        integer i, k, l, m
        logical is_empty
        double precision xcord, ycord, zcord, radius
        double precision hloc, centerloc, gravityAN,gradPhi
        double precision gradPhiAn,gravity_xz_abs,theta
        double precision gravity_x,gravity_y,gravity_z
        character*32 filenamexy, filenamexz

201     format(15(1x,f21.14))

        ! creating output files
        write(filenamexy, 101) timer
101     format(f6.3,"_xy.dat")
        open(501, file = filenamexy)

        write(filenamexz, 102) timer
102     format(f6.3,"_xz.dat")
        open(502, file = filenamexz)

        ! writing the results to ascii files
        do m=1,Ndepth

         hloc      = h     /(2.d0**(m-1))
         centerloc = length/(2.d0**m)

         ! saving the XY cut (at z=0)
         l = N/2 + 1
         do k=2,N+1
          do i=2,N+1

           ! current cell is empty
           is_empty = (m < Ndepth)
     =                .and. (i > N/4+1)
     =                .and. (i < 3*N/4+2)
     =                .and. (k > N/4+1)
     =                .and. (k < 3*N/4+2)
     =                .and. (l > N/4+1)
     =                .and. (l < 3*N/4+2)

           ! current cell is not empty
           if( .not. is_empty  ) then
            xcord = (i-1)*hloc - centerloc - hloc/2.d0
            ycord = (k-1)*hloc - centerloc - hloc/2.d0
            zcord = (l-1)*hloc - centerloc - hloc/2.d0

            radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)   ! radial distance from the coordinate center

            gravityAN=gravity_profile_an(i,k,l,m)                     ! analytical potential

            gradPhi=dsqrt(                                            ! module of the acceleration 
     =        ((gravity(i+1,k,l,m)-gravity(i-1,k,l,m))/(hloc*2))**2
     =      + ((gravity(i,k+1,l,m)-gravity(i,k-1,l,m))/(hloc*2))**2
     =      + ((gravity(i,k,l+1,m)-gravity(i,k,l-1,m))/(hloc*2))**2 )

            gradPhiAn=gradient_profile_an(i,k,l,m)                    ! analytic acceleration

            if (m.eq.1) then 
              If (k.eq.2.or.k.eq.N+1.or.i.eq.2.or.i.eq.N+1) then      ! skipping the outer boundary of the coarsest grid

              else
                 write(501,201) xcord,ycord,density(i,k,l,m),
     =           gravity(i,k,l,m),gravityAN,
     =           100*abs(gravity(i,k,l,m)-gravityAN)/abs(gravityAN),
     =           gradPhi,gradPhiAn, 
     =           100*abs(gradPhi-gradPhiAn)/abs(gradPhiAn)
             endif

            Else
                write(501,201) xcord,ycord,density(i,k,l,m),
     =          gravity(i,k,l,m),gravityAN,
     =          100*abs(gravity(i,k,l,m)-gravityAN)/abs(gravityAN),
     =          gradPhi,gradPhiAn,
     =          100*abs(gradPhi-gradPhiAn)/abs(gradPhiAn)

           endif

           endif
          enddo
         enddo

         ! saving the XZ-cut (Y=0)
         k = N/2 + 1
         do l=2,N+1
          do i=2,N+1

           ! current cell is empty
           is_empty = (m < Ndepth)
     =                .and. (i > N/4+1)
     =                .and. (i < 3*N/4+2)
     =                .and. (k > N/4+1)
     =                .and. (k < 3*N/4+2)
     =                .and. (l > N/4+1)
     =                .and. (l < 3*N/4+2)

           ! current cell is not empty
           if ( .not. is_empty  ) then
           xcord = (i-1)*hloc - centerloc - hloc/2.d0
           ycord = (k-1)*hloc - centerloc - hloc/2.d0
           zcord = (l-1)*hloc - centerloc - hloc/2.d0

           radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)

           gravityAN=gravity_profile_an(i,k,l,m)

           gravity_x=-(gravity(i+1,k,l,m)-gravity(i-1,k,l,m))/(hloc*2)  ! acceleration components along x, y, and z
           gravity_y=-(gravity(i,k+1,l,m)-gravity(i,k-1,l,m))/(hloc*2)
           gravity_z=-(gravity(i,k,l+1,m)-gravity(i,k,l-1,m))/(hloc*2)

           gradPhi=dsqrt(gravity_x**2 + gravity_y**2 + gravity_z**2)
           gradPhiAn=gradient_profile_an(i,k,l,m)

           gravity_xz_abs=dsqrt(                                        ! gravitational acceleration in the xz-plane
     =        ((gravity(i+1,k,l,m)-gravity(i-1,k,l,m))/(hloc*2))**2
     =      + ((gravity(i,k,l+1,m)-gravity(i,k,l-1,m))/(hloc*2))**2 )

            If (gravity_z.ge.0) then                                    ! calculating the angle between the acceleration vector in the xz-plane and the x-coordinate
                If (gravity_xz_abs.eq.0) then
                   theta=0.d0
                Else
                   theta=dacos(gravity_x/gravity_xz_abs)
                Endif
            Else
                If (gravity_xz_abs.eq.0) then
                   theta=0.d0
                else
                   theta=2.d0*pi-dacos(gravity_x/gravity_xz_abs)
                endif
            Endif 

           if (m.eq.1) then
              If (l.eq.2.or.l.eq.N+1.or.i.eq.2.or.i.eq.N+1) then

              else
              write(502,201) xcord,zcord,density(i,k,l,m),
     =                     gravity(i,k,l,m),gravityAN,
     =          100*abs(gravity(i,k,l,m)-gravityAN)/abs(gravityAN),
     =          gradPhi,gradPhiAn,
     =          100*abs(gradPhi-gradPhiAn)/abs(gradPhiAn),
     =          gravity_xz_abs,theta*180/pi                           ! these two values can be used to visualize the gravitational acceleration field
              endif

           Else
              write(502,201) xcord,zcord,density(i,k,l,m),
     =                     gravity(i,k,l,m),gravityAN,
     =          100*abs(gravity(i,k,l,m)-gravityAN)/abs(gravityAN),
     =          gradPhi,gradPhiAn,
     =          100*abs(gradPhi-gradPhiAn)/abs(gradPhiAn),
     =          gravity_xz_abs,theta*180/pi

          endif
 
          Endif
          enddo
         enddo

        enddo

        close(501)
        close(502)

      end subroutine save
      !**************************************************