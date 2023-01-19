      !**************************************************
      ! analytic solutions for the gravitational potential
      ! note that an analytic solution for a tigth binary is not provided
      ! part of Poisson solver, the convolution method (Vorobyov, McKevitt, et al. 2023)
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !**************************************************
      double precision function gravity_profile_an(i,k,l,m)
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
        double precision radius1, radius2

        ! compute actual h and center for local mesh
        hloc      = h     /(2.d0**(m-1))
        centerloc = length/(2.d0**m)

        ! compute coordinates
        xcord = (i-1)*hloc - centerloc - hloc/2.d0
        ycord = (k-1)*hloc - centerloc - hloc/2.d0
        zcord = (l-1)*hloc - centerloc - hloc/2.d0

! ============ uniform sphere =====================
        If (sphere) then
        ! return exact gravity
        radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)
           if(radius <= rsp) then
              gravity_profile_an=-2*pi*GravConst*rho_sp*(rsp**2 -
     =                radius**2/3.d0)
           else
              gravity_profile_an= -4.d0*pi*GravConst*rho_sp*
     =                rsp**3/(3.d0*radius)
           endif
        Endif

! ============ uniform oblate ellipsoid ==================
        If (ellipsoid) then

         radius = dsqrt(xcord*xcord+ycord*ycord)

         if (radius**2/rsp1**2+zcord**2/rsp3**2.le.1.d0) then

            E=dsqrt(1-(rsp3/rsp1)**2)
            AA1=( dsqrt(1-E**2)/E**3 )*dasin(E)-(1-E**2)/E**2
            AA3=2./E**2 - 2.*( dsqrt(1-E**2)/E**3 )*dasin(E)

            gravity_profile_an = - pi*GravConst*rho_sp *
     =              (rsp1**2*AA1+rsp1**2*AA1+rsp3**2*AA3 -
     =               AA1*radius**2 - AA3*zcord**2 )

         else

          prom1=radius**2+zcord**2 - rsp1**2-rsp3**2
          prom2=4.d0*(rsp1**2*rsp3**2-radius**2*rsp3**2 -
     =          zcord**2*rsp1**2)
          lambda=(prom1+dsqrt(prom1**2-prom2))/2
          I1=pi/Dsqrt(rsp1**2-rsp3**2) -
     =       2.d0/Dsqrt(rsp1**2-rsp3**2)*
     =       datan(Dsqrt( (rsp3**2+lambda)/(rsp1**2-rsp3**2) ))

          sqroot=Dsqrt(rsp3**2+lambda)

          gravity_profile_an=-pi*GravConst*rho_sp*rsp1**2*rsp3*
     =    ( (1.0d0 + radius**2/(2.0d0*(rsp3**2-rsp1**2)) -
     =    zcord**2/(rsp3**2-rsp1**2) )*I1 -
     =    radius**2*sqroot/((rsp3**2-rsp1**2)*(rsp1**2+lambda)) -
     =    zcord**2*( 2.d0/((rsp1**2+lambda)*sqroot) -
     =       2.d0*sqroot/((rsp3**2-rsp1**2)*(rsp1**2+lambda)) ) )

          endif

        Endif

! ============ a wide-separation binary system =====================
        If (binary) then

           radius1 = dsqrt((xcord-xb1)*(xcord-xb1) +
     =               ycord*ycord + zcord*zcord)

           radius2 = dsqrt((xcord-xb2)*(xcord-xb2) + ycord*ycord +
     =                      zcord*zcord)

           if(radius1 <= rsp_b1) then
          gravity_profile_an_b1=-2*pi*GravConst*rho_sp1*(rsp_b1**2-
     =                radius1**2/3.d0)
           else
              gravity_profile_an_b1=-4.d0*pi*GravConst*rho_sp1*
     =                rsp_b1**3/(3.d0*radius1)
           endif

           if(radius2 <= rsp_b2) then
          gravity_profile_an_b2=-2*pi*GravConst*rho_sp2*(rsp_b2**2-
     =                radius2**2/3.d0)
           else
              gravity_profile_an_b2=-4.d0*pi*GravConst*rho_sp2*
     =                rsp_b2**3/(3.d0*radius2)
           endif

           gravity_profile_an = gravity_profile_an_b1 +
     =                          gravity_profile_an_b2

        Endif

! =========density distribution after Wang & Yen 2020 =======================
        If (wang) then

            radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)

           if (radius.le.r_wang) then
               gravity_profile_an = -2.d0/3.d0*pi*rho_wang*r_wang**2
     =         + 4*pi*rho_wang*(radius**2/6.d0
     =         - 0.1d0*radius**4/r_wang**2
     =         + 1.d0/42.d0*radius**6/r_wang**4)
           else
               gravity_profile_an =
     =          - 32*pi*rho_wang*r_wang**3/105.d0 / radius
           endif

        Endif

      end function gravity_profile_an
      !**************************************************