      !**************************************************
      ! initial density profiles
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !**************************************************
      double precision function density_profile(i,k,l,m)
        use machine
        use mesh
        implicit none

        ! input index of cell
        integer, intent(in) :: i, k, l, m

        ! variables
        double precision xcord, ycord, zcord, radius
        double precision hloc, centerloc
        double precision radius1,radius2
        double precision density_profile_1,density_profile_2

        ! computing the cell size and the center of a particular nested mesh
        hloc      = h     /(2.d0**(m-1))
        centerloc = length/(2.d0**m)

        ! computing coordinates on a particular mesh
        xcord = (i-1)*hloc - centerloc - hloc/2.d0
        ycord = (k-1)*hloc - centerloc - hloc/2.d0
        zcord = (l-1)*hloc - centerloc - hloc/2.d0

        If (dipole) then
          density_profile=0.d0
          if (m.eq.Ndepth) then
              if (i.eq.N/2+2.and.k.eq.N/2+2.and.l.eq.N/2+2) then
                    density_profile=rho_dp1
                endif
              if (i.eq.N/2+2.and.k.eq.N/2+2.and.l.eq.N/2+3) then
                    density_profile=rho_dp2
                endif
          endif
        Endif

        If (sphere) then
           radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)
           if(radius <= rsp) then
             density_profile = rho_sp
           else
             density_profile = eps
           endif
        Endif

        If (ellipsoid) then
           radius = dsqrt(xcord*xcord + ycord*ycord)
           if((radius**2/rsp1**2+zcord**2/rsp3**2.le.1.d0)) then
             density_profile = rho_sp
           else
             density_profile = eps
           endif
        Endif

        If (binary) then
           radius1 = dsqrt((xcord-xb1)*(xcord-xb1) +
     =               ycord*ycord + zcord*zcord)

           radius2 = dsqrt((xcord-xb2)*(xcord-xb2) + ycord*ycord +
     =                      zcord*zcord)

           if(radius1 <= rsp_b1) then
             density_profile_1 = rho_sp1
           else
             density_profile_1 = eps
           endif

           if(radius2 <= rsp_b2) then
             density_profile_2 = rho_sp2
           else
             density_profile_2 = eps
           endif

           density_profile=density_profile_1+density_profile_2

        Endif

        If (wang) then

           radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)

           if (radius.le.r_wang) then
               density_profile = rho_wang*(1-(radius/r_wang)**2)**2
           else
               density_profile = 0.d0
           endif

        endif

      end function density_profile
      !**************************************************