      !****************************************************
      ! Poisson solver, the convolution method, after Vorobyov, McKevitt, et al. 2023.
      ! Copyright © 2023 Eduard Vorobyov and James McKevitt. All rights reserved.
      !****************************************************
      program Convolution
        ! using modules
        use machine
        use mesh
        use PoissonIntegral
        implicit none

        double precision timer, tend
        integer i, m, nn
        !****************************************************
        ! choose the version: basic or advanced !
        ! basic - see sections 3-6              !
        ! advanced - see section 7              !
        !****************************************************
        logical, parameter :: advanced=.false.

        !****************************************************
        ! choose the summation method: fft or direct summation
        ! see Appendix B
        ! hint: choose FFT if N>=32
        logical, parameter :: direct=.false.

        !****************************************************
        ! Main code
        !****************************************************

        timer = 0.d0                   ! set the timer
        call load                      ! load initial density configuration
        Call OMP_SET_NUM_Threads(Nthreads)

        do m = Ndepth, 2, -1           ! projection of density from m to (m-1) grid level
            call projection(m)
        enddo

        call Integration               ! calculation of the integral (13), sect. 3. If the grid is static, has to be done only once

        Do m=1,Ndepth
           call SetGravity(m)          ! calculating the FFT transform of the inverse distance (see eq. 10)
           If (advanced) then
                if (m.gt.1.and..not.direct)
     =              call setGravityDBL(m)  ! FFT transform of the inverse distance on a doubled mesh (see sect. 7)
           endif
        Enddo

        IF (advanced) then

          If (Ndepth.eq.1) then
             call SolveGravity(Ndepth)          ! calculating the gravitational potential  on the single mesh (see eq.10)
          Else   
          do m = 2, Ndepth
               call SolveGravityCut(m-1)        ! calculating the gravitational potential on the mesh with a central hole (see Fig. 5)
               If (direct) then 
                   call DirectSummation(m-1)    ! calculating the input from coarser meshes using direct summation (see eq. B.7)
               Else    
                   call DoubleMeshNew(m)        ! initializing a doubled mesh for a particular mesh level m (Omega(m))
                   call SolveGravityDBL(m)      ! solving the potential on a doubled mesh (see Fig. 14)
                   call ProjectionGravityDBL(m) ! adding potentials from Omega_dbl(m) and Omega_hole(m-1)
               Endif    
               call SolveGravity(m)             ! calculating the gravitational potential on a particular mesh level m , Omega(m) (see eq.10, Fig. 14)
               Call gravity_reduction(m-1)      ! projection of the potential from a finer grid onto a coarser grid. Part of the interpolation procedure in deriving \tilde\Phi (e.g. eq. 14)
               Call AddPotential(m)             ! adding the potentials from the decomposed grid together (see Fig. 5 and eqs. 14)
          enddo
          Endif

        ELSE
          m=1
          call SolveGravity(m)                  ! calculating the gravitational potential  on the single mesh (see eq.10)
          do m = 2, Ndepth
                call SolveGravityCut(m-1)       ! calculating the gravitational potential on the mesh with a central hole (see Fig. 5)
                call SolveGravity(m)            ! calculating the gravitational potential on a particular mesh level m , Omega(m) (see eq.10, Fig. 14)
                Call gravity_reduction(m-1)     ! projection of the potential from a finer grid onto a coarser grid. Part of the interpolation procedure in deriving \tilde\Phi (e.g. eq. 14)
                Call AddPotential(m)            ! adding the potentials from the decomposed grid together (see Fig. 14 and eqs. 15, and also B1-B2)
          enddo

        ENDIF

        do m = Ndepth, 2, -1                    ! projection of the potential from the m-grid to (m-1)-grid
           call projection_gravity(m)
        enddo  
        
        Do m=1,Ndepth
          call poisson_boundary(m)              ! applying boundary conditions at the nested grid interfaces. Required when calculating the gravitational acceleration
        enddo

        call save(timer)                        ! final output

      end program Convolution