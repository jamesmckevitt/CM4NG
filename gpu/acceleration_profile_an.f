      !**************************************************
      ! analytical solutions for the gravitational accelerations
      ! note that an analytic solution for a tigth binary is not provided
      ! an analytic solution for a uniform sphere can be found in Stone & Norman 1992
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !**************************************************
      double precision function gradient_profile_an(i,k,l,m)
        use machine
        use mesh
        implicit none

        ! input index of cell
        integer, intent(in) :: i, k, l, m

        ! variables
        double precision xcord, ycord, zcord, radius
        double precision hloc, centerloc, I0

        double precision prom1,prom2,lambda,I1,sqroot
        double precision AA1,AA3,E

        double precision gravity_profile_an_b1,gravity_profile_an_b2
        double precision gradient_profile_an_b1_x,
     =  gradient_profile_an_b1_y,gradient_profile_an_b1_z
        double precision gradient_profile_an_b2_x,
     =  gradient_profile_an_b2_y,gradient_profile_an_b2_z
        double precision radius1, radius2,prom3,prom4,prom5

        ! compute actual h and center for local mesh
        hloc      = h     /(2.d0**(m-1))
        centerloc = length/(2.d0**m)

        ! compute coordinates
        xcord = (i-1)*hloc - centerloc - hloc/2.d0
        ycord = (k-1)*hloc - centerloc - hloc/2.d0
        zcord = (l-1)*hloc - centerloc - hloc/2.d0

! ============ uniform oblate ellipsoid ==================
        If (ellipsoid) then

         radius = dsqrt(xcord*xcord+ycord*ycord)

         if (radius**2/rsp1**2+zcord**2/rsp3**2.le.1.d0) then

            E=dsqrt(1-(rsp3/rsp1)**2)
            AA1=( dsqrt(1-E**2)/E**3 )*dasin(E)-(1-E**2)/E**2
            AA3=2./E**2 - 2.*( dsqrt(1-E**2)/E**3 )*dasin(E)

            gradient_profile_an = 2*pi*GravConst*rho_sp *
     =          dsqrt(AA1**2*radius**2 + AA3**2*zcord**2)

         else

          prom1=radius**2+zcord**2 - rsp1**2-rsp3**2
          prom2=4.d0*(rsp1**2*rsp3**2-radius**2*rsp3**2 -
     =          zcord**2*rsp1**2)
          lambda=(prom1+dsqrt(prom1**2-prom2))/2

          I1=pi/Dsqrt(rsp1**2-rsp3**2) -
     =       2.d0/Dsqrt(rsp1**2-rsp3**2)*
     =       datan(Dsqrt( (rsp3**2+lambda)/(rsp1**2-rsp3**2) ))

          prom1 = radius*I1/(2.d0*(rsp3**2-rsp1**2))
          prom2 = radius*dsqrt(rsp3**2+lambda)/
     =           ((rsp3**2-rsp1**2)*(rsp1**2+lambda))
          prom3 = zcord*I1/(rsp3**2-rsp1**2)
          prom4 = 2.d0/((rsp1**2+lambda)*dsqrt(rsp3**2+lambda))
          prom5 = 2.d0*dsqrt(rsp3**2+lambda)/
     =            ((rsp3**2-rsp1**2)*(rsp1**2+lambda))

         gradient_profile_an=2*pi*GravConst*rho_sp*rsp1**2*rsp3*
     =   dsqrt( (prom1-prom2)**2 + (prom3+zcord*(prom4-prom5))**2 )

          endif

        Endif

! ============ a wide-separation binary system =====================
        If (binary) then

           radius1 = dsqrt((xcord-xb1)*(xcord-xb1) +
     =               ycord*ycord + zcord*zcord)
           
           radius2 = dsqrt((xcord-xb2)*(xcord-xb2) + ycord*ycord +
     =                      zcord*zcord)

           if(radius1 <= rsp_b1) then
              gradient_profile_an_b1_x=-4.d0/3.d0*pi*GravConst*
     =                               rho_sp1*(xcord-xb1)            !radius1
              gradient_profile_an_b1_y=-4.d0/3.d0*pi*GravConst*
     =                               rho_sp1*ycord                  !radius1
              gradient_profile_an_b1_z=-4.d0/3.d0*pi*GravConst*
     =                               rho_sp1*zcord                  !radius1
           else
              gradient_profile_an_b1_x=-4.d0*pi*GravConst*rho_sp1*
     =                rsp_b1**3/(3.d0*radius1**3)*(xcord-xb1)
              gradient_profile_an_b1_y=-4.d0*pi*GravConst*rho_sp1*
     =                rsp_b1**3/(3.d0*radius1**3)*ycord
              gradient_profile_an_b1_z=-4.d0*pi*GravConst*rho_sp1*
     =                rsp_b1**3/(3.d0*radius1**3)*zcord
           endif

           if(radius2 <= rsp_b2) then
              gradient_profile_an_b2_x=-4.d0/3.d0*pi*GravConst*
     =                               rho_sp2*(xcord-xb2)         !radius2
              gradient_profile_an_b2_y=-4.d0/3.d0*pi*GravConst*
     =                               rho_sp2*ycord               !radius2
              gradient_profile_an_b2_z=-4.d0/3.d0*pi*GravConst*
     =                               rho_sp2*zcord               !radius2
           else
              gradient_profile_an_b2_x=-4.d0*pi*GravConst*rho_sp2*
     =                rsp_b2**3/(3.d0*radius2**3)*(xcord-xb2)
              gradient_profile_an_b2_y=-4.d0*pi*GravConst*rho_sp2*
     =                rsp_b2**3/(3.d0*radius2**3)*ycord
              gradient_profile_an_b2_z=-4.d0*pi*GravConst*rho_sp2*
     =                rsp_b2**3/(3.d0*radius2**3)*zcord
           endif

           gradient_profile_an = dsqrt(
     =      (gradient_profile_an_b1_x + gradient_profile_an_b2_x)**2 +
     =      (gradient_profile_an_b1_y + gradient_profile_an_b2_y)**2 +
     =      (gradient_profile_an_b1_z + gradient_profile_an_b2_z)**2)

        Endif

! =========density distribution after Wang & Yen 2020 =======================
        If (wang) then

            radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)

           if (radius.le.r_wang) then
               gradient_profile_an = 4*pi*rho_wang*(radius/3.d0
     =         - 2.d0/5.d0*radius**3/r_wang**2
     =         + 1.d0/7.d0*radius**5/r_wang**4)
           else
               gradient_profile_an =  
     =         32*pi*rho_wang*r_wang**3/105.d0/radius**2
           endif

        Endif

      end function gradient_profile_an
      !**************************************************