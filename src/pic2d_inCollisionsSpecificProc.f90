
!---------------------------------------------------------------------------------------------------------------------
!
subroutine calculate_thermal_cx_probab

  USE CurrentProblemValues, ONLY: delta_t_s, N_subcycles, e_Cl, kB_JK, m_e_kg
  USE MCCollisions
  USE IonParticles, ONLY : N_spec, Ms

  implicit none

  integer s, n
  real(8) sigma_m2_1eV, alpha, Tgas_eV, ngas_m3, sigma_m2_therm, Vmean_ms
!function
  real(8) sigma_rcx_m2

  if (no_rcx_collisions) return

  DO s = 1, N_spec
     if (.not.collision_rcx(s)%rcx_on) cycle
     n = collision_rcx(s)%neutral_species_index
     sigma_m2_1eV = neutral(n)%sigma_rcx_m2_1eV
     alpha = neutral(n)%alpha_rcx
     Tgas_eV = neutral(n)%T_K * kB_JK / e_Cl
     ngas_m3 = neutral(n)%N_m3
     sigma_m2_therm = sigma_rcx_m2(1.5_8 * Tgas_eV, sigma_m2_1eV, alpha)
     Vmean_ms = sqrt(2.54647908947033_8 * Tgas_eV * e_Cl / (Ms(s) * m_e_kg))   ! 2.5464... = 8/pi, note that Ms(s) is (amu_kg * M_i_amu(s)) / m_e_kg
     collision_rcx(s)%probab_thermal = ngas_m3 * sigma_m2_therm * Vmean_ms * delta_t_s * N_subcycles
  END DO

end subroutine calculate_thermal_cx_probab

!------------------------------------------------------------------------- 
!
real(8) function sigma_rcx_m2(energy_eV, sigma_m2_1eV, alpha)

  implicit none

  real(8) energy_eV, sigma_m2_1eV, alpha, Tgas_eV

  sigma_rcx_m2 = sigma_m2_1eV * (1.0_8 - alpha * log(MAX(energy_eV,0.001_8)))**2   ! the limiter avoids singularity
                                                                                   ! for comparison, 300 K = 0.026 eV

end function sigma_rcx_m2

!---------------------------------------------------------------------------------------------------------------------
!
subroutine PERFORM_RESONANT_CHARGE_EXCHANGE

  USE MCCollisions
  USE IonParticles
  USE CurrentProblemValues, ONLY : energy_factor_eV, delta_t_s, N_subcycles, V_scale_ms, T_cntr
  USE rng_wrapper

  IMPLICIT NONE

  INTEGER s, n, i
  real(8) ngas_m3, sigma_m2_1eV, alpha, probab_rcx_therm_2
  real(8) factor_eV, vfactor, prob_factor

  real(8) vx, vy, vz, energy_eV, vr_ms
  real(8) probab_rcx
  real(8) vxn, vyn, vzn, vr

! functions
  real(8) neutral_density_normalized, sigma_rcx_m2
  
  if (no_rcx_collisions) return

! clear collision counters
  DO s = 1, N_spec
     collision_rcx(s)%counter = 0
  END DO

  DO s = 1, N_spec

    if (.not.collision_rcx(s)%rcx_on) cycle

    n = collision_rcx(s)%neutral_species_index

    ngas_m3 = neutral(n)%N_m3
    sigma_m2_1eV = neutral(n)%sigma_rcx_m2_1eV
    alpha =        neutral(n)%alpha_rcx
    probab_rcx_therm_2  = (collision_rcx(s)%probab_thermal)**2

    factor_eV = Ms(s) * energy_factor_eV         ! instead of collision_rcx(s)%factor_eV
    vfactor = collision_rcx(s)%vfactor           ! to convert Maxwellian sample
    prob_factor = ngas_m3 * delta_t_s * N_subcycles

    DO i = 1, N_ions(s)
      ! create a virtual neutral particle
      call GetMaxwellVelocity(vxn)
      call GetMaxwellVelocity(vyn)
      call GetMaxwellVelocity(vzn)
      vxn = vxn * vfactor
      vyn = vyn * vfactor
      vzn = vzn * vfactor

      vx = ion(s)%part(i)%VX
      vy = ion(s)%part(i)%VY
      vz = ion(s)%part(i)%VZ

      ! relative velocity between colliding particles
      vr = sqrt((vx-vxn)**2 + (vy-vyn)**2 + (vz-vzn)**2)
      energy_eV = vr * factor_eV
      vr_ms = vr * V_scale_ms 

      probab_rcx = prob_factor * vr_ms * sigma_rcx_m2(energy_eV, sigma_m2_1eV, alpha)
      probab_rcx = neutral_density_normalized(n, ion(s)%part(i)%x, ion(s)%part(i)%y) * probab_rcx  ! account for the nonuniform density and the low-energy correction
      
      if (ngas_m3.le.1E21) then ! less than about 4 Pa? ... 
        if ((delta_t_s*T_cntr).le.10E-6) then ! ... then first have it running at higher pressure for some time to remove waves.
          probab_rcx = ((delta_t_s*T_cntr)/10E-6) * probab_rcx + probab_rcx/ngas_m3 * 1E21 * (1 - (delta_t_s*T_cntr)/10E-6) ! linear decrease to final pressure over time
        end if
      end if

      if (ngas_m3.ge.3E21) then ! more than about 10 Pa? ... 
        if ((delta_t_s*T_cntr).le.6E-6) then ! ... then first have it running at lower pressure for some time have it converge faster.
          probab_rcx = ((delta_t_s*T_cntr)/6E-6) * probab_rcx + probab_rcx/ngas_m3 * 2E21 * (1 - (delta_t_s*T_cntr)/10E-6) ! linear decrease to final pressure over time
        end if
      end if

      if (well_random_number().le.probab_rcx) then
        ion(s)%part(i)%VX = vxn
        ion(s)%part(i)%VY = vyn
        ion(s)%part(i)%VZ = vzn
        collision_rcx(s)%counter = collision_rcx(s)%counter + 1
      end if
    END DO
  END DO   

END SUBROUTINE PERFORM_RESONANT_CHARGE_EXCHANGE


subroutine PERFORM_ION_NEUTRAL_COLLISION
! Julian Held, jheld@umn.edu, 2024
! After: J Trieschmann, PHD Thesis, 2015:
! https://hss-opus.ub.ruhr-uni-bochum.de/opus4/frontdoor/deliver/index/docId/5307/file/diss.pdf
! partly based on: Nanbu and Kitatani: J. Phys. D: Appl. Phys. 28 (1995) 324-330

  USE MCCollisions
  USE IonParticles
  USE CurrentProblemValues, ONLY : energy_factor_eV, delta_t_s, N_subcycles, V_scale_ms, T_cntr, pi, e_Cl, true_eps_0_Fm, amu_kg, kB_JK
  USE rng_wrapper
  !use stdlib_specialfunctions_gamma, only: gamma

  IMPLICIT NONE
  INCLUDE 'mpif.h'


  INTEGER s, n, i
  INTEGER ierr
  LOGICAL CX ! did charge exchange occur?
  INTEGER case ! Trieschmann's three cases, 1 -> A, 2 -> B, 3 -> C
  real(8) ngas_m3, sigma_m2_1eV, probab_rcx_therm_2
  real(8) factor_eV, vfactor, prob_factor

  real(8) vx, vy, vz, vx_, vy_, vz_, Ekin, vr_ms
  real(8) probab_rcx
  real(8) vxn, vyn, vzn, vxn_, vyn_, vzn_ ! neutral velocities before and after
  real(8) Rx, Ry, Rz
  real(8) gx, gy, gz, g, g_perp, hx, hy, hz
  real(8) Mi, Mn
  real(8) E_ratio, E_ratio_ion

  real(8) neutral_density_normalized, sigma_rcx_m2 ! functions

  real(8) beta_inf, dref, omega, Tref, alpha0, k1, k2, vth, q, mr, t_func, bcx
  real(8) sigmaL, sigmaP, sigmaT, sigma_cx, d0
  real(8) p_col, p_cx, p, bmax_col, bmax_cx, bmax, b
  real(8) xi0, xi1, xi, chi, beta, beta0, theta, theta0, Fel, dtheta
  real(8) cos_chi, sin_chi, phi

  beta_inf = 4 ! impact parameter cutoff

  ! Constants (for Ar) from Trieschmann, PHD thesis, p 19
  dref = 4.614D-10    
  omega = 0.721
  Tref = 273.0    
  alpha0 = 1.6411D-30

  ! Constants for charge exchange (Ar) Trieschmann, PHD thesis, p 67
  k1 = 0.525D-11
  k2 = 1.007D-9

  vth = 15000.0 ! threshold velocity above which collision acts like neutral-neutral col

  if (no_rcx_collisions) return

  ! clear collision counters
  DO s = 1, N_spec
      collision_rcx(s)%counter = 0
  END DO
  
  DO s = 1, N_spec
    if (.not.collision_rcx(s)%rcx_on) cycle
    n = collision_rcx(s)%neutral_species_index
    ngas_m3 = neutral(n)%N_m3
    Mn = neutral(n)%M_amu * amu_kg
    vfactor = collision_rcx(s)%vfactor

    DO i = 1, N_ions(s)
      CX = .False.
      ! create a virtual neutral particle
      call GetMaxwellVelocity(vxn)
      call GetMaxwellVelocity(vyn)
      call GetMaxwellVelocity(vzn)
      vxn = vxn * vfactor * V_scale_ms
      vyn = vyn * vfactor * V_scale_ms
      vzn = vzn * vfactor * V_scale_ms

      vx = ion(s)%part(i)%VX * V_scale_ms
      vy = ion(s)%part(i)%VY * V_scale_ms
      vz = ion(s)%part(i)%VZ * V_scale_ms

      q = Qs(s) * e_Cl 
      Mi = M_i_amu(s) * amu_kg

      ! relative velocity between colliding particles
      gx = vxn-vx
      gy = vyn-vy
      gz = vzn-vz
      g = sqrt((gx)**2 + (gy)**2 + (gz)**2)
      g_perp = sqrt(gy**2 + gz**2)

      mr = Mi*Mn/(Mi+Mn) 
      Ekin = 0.5 * mr * g**2

      t_func = 3*((g - 0.2*vth)/(0.8*vth))**2 - 2*((g - 0.2*vth)/(0.8*vth))**3

      sigmaL = sqrt(pi*alpha0*(q**2)/(true_eps_0_Fm * mr))/g
      d0 = dref * sqrt((((kB_JK * Tref)/Ekin)**(omega-0.5)) * 1/gamma(5.0/2.0 - omega)) 
      sigmaT = pi * d0**2

      ! Trieschmann's A,B,C cases, splitting by velocity:
      ! A -> slow, Nanbus model
      ! B -> half and half
      ! C -> fast, M1 elleastic scattering (like neutral-neutral col)
      if (g.LT.(0.2*vth)) then ! A  
        case = 1
        sigmaP = sigmaL * beta_inf**2
      end if

      if ((g.GE.(0.2*vth)).AND.(g.LT.vth)) then ! B
        case = 2
        sigmaP = sigmaL * beta_inf**2 + sigmaT * t_func
      end if

      if (g.GE.vth) then ! C
        case = 3
        sigmaP = sigmaT
      end if

      ! charge exchange 
      bcx = -k1 * log(g) + k2
      sigma_cx = 0.5 * pi * bcx**2
   
      ! Larger cross section gives total probalility for joined evaluation of CX
      ! and ellastic collision
      p_col = g * sigmaP * ngas_m3 * delta_t_s * N_subcycles
      p_cx = g * sigma_cx * ngas_m3 * delta_t_s * N_subcycles
      p = max((2.0*p_cx), p_col)

      ! <-------------------- probe project specific adjustments ------------------------>
      if (ngas_m3.le.1E21) then ! less than about 4 Pa? ... 
        if ((delta_t_s*T_cntr).le.10E-6) then ! ... then first have it running at higher collision probability for some time to remove waves.
          p = ((delta_t_s*T_cntr)/10E-6) * p + p/ngas_m3 * 1E21 * (1 - (delta_t_s*T_cntr)/10E-6) ! linear decrease to final collision probability over time
        end if
      end if

      if (ngas_m3.ge.3E21) then ! more than about 10 Pa? ... 
        if ((delta_t_s*T_cntr).le.6E-6) then ! ... then first have it running at lower collision probability for some time for faster convergence.
          p = ((delta_t_s*T_cntr)/6E-6) * p + p/ngas_m3 * 2E21 * (1 - (delta_t_s*T_cntr)/10E-6)  ! linear increase to final collision probability over time
        end if
      end if
      ! <------------------- end project specific adjustments --------------------------->

      if (well_random_number().le.p) then ! perform collision
        collision_rcx(s)%counter = collision_rcx(s)%counter + 1
        bmax_col = sqrt(sigmaP/pi)
        bmax_cx = sqrt(sigma_cx/pi)
        bmax = max(bmax_col,bmax_cx)
        b = bmax * sqrt(well_random_number())

        phi = 2.0*pi*well_random_number() ! phi is always random
        hx = g_perp * cos(phi)
        hy = -(gy*gx*cos(phi) + g*gz*sin(phi)) / g_perp                                                
        hz = -(gz*gx*cos(phi) - g*gy*sin(phi)) / g_perp

        IF (case.EQ.2) THEN ! Split case B randomly in A and C buckets
          if (well_random_number().le.(sigmaT * t_func / sigmaP)) then ! -> case A
            case = 1
          else ! -> case C
            case = 3
          end if
        END IF

        IF (case.EQ.1) THEN ! collision for case A
          beta = b*(2*pi*true_eps_0_Fm*Ekin/(alpha0*q**2))**0.25
          beta0 = 1.001 
          if (beta.GE.beta0) then ! beta > 1  ->  polarization scattering
            xi0 = sqrt(beta**2 - sqrt(beta**4 - 1))
            xi1 = sqrt(beta**2 + sqrt(beta**4 - 1))
            xi = xi0/xi1

            Fel = 0.0
            dtheta = 0.0002 ! unsure about stepsize TODO check
            do theta = 0, pi/2.0, dtheta ! incomplete elliptic integral of first kind (note missprint in Trieschmann's thesis)
              Fel = Fel + 1.0/sqrt((1 - (xi**2)*sin(theta)**2))*dtheta
            end do

            theta0 = sqrt(2.0)*beta/xi1 * Fel
            chi = pi - 2.0 * theta0
            cos_chi = cos(chi)
            sin_chi = abs(sin(chi))
            
            ! post collision velocities
            vx_ =  vx +  Mn/(Mi + Mn) * (gx*(1-cos_chi) + hx*sin_chi)
            vy_ =  vy +  Mn/(Mi + Mn) * (gy*(1-cos_chi) + hy*sin_chi)
            vz_ =  vz +  Mn/(Mi + Mn) * (gz*(1-cos_chi) + hz*sin_chi)
            vxn_ = vxn - Mi/(Mi + Mn) * (gx*(1-cos_chi) + hx*sin_chi)
            vyn_ = vyn - Mi/(Mi + Mn) * (gy*(1-cos_chi) + hy*sin_chi)
            vzn_ = vzn - Mi/(Mi + Mn) * (gz*(1-cos_chi) + hz*sin_chi)

            if ((well_random_number().le.0.5).AND.(b.LE.bcx)) then  ! CX: identity switch (50% chance if b<bcx)
              CX = .True.
            end if

          end if  

          if (beta.LT.beta0) then ! spiraling motion, VHS (variable hard sphere) model
            ! VHS model -> random direction
            ! https://stackoverflow.com/questions/5408276/
            ! https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
            theta = acos(1.0 - 2.0*well_random_number())
            Rx = sin(theta) * cos(phi)
            Ry = sin(theta) * sin(phi)
            Rz = cos(theta)
            
            ! post collision velocities
            vx_ =  (1/(Mi + Mn)) * (Mi*vx + Mn*vxn - Mn*g*Rx)
            vy_ =  (1/(Mi + Mn)) * (Mi*vy + Mn*vyn - Mn*g*Ry)
            vz_ =  (1/(Mi + Mn)) * (Mi*vz + Mn*vzn - Mn*g*Rz)
            vxn_ = (1/(Mi + Mn)) * (Mi*vx + Mn*vxn + Mi*g*Rx)
            vyn_ = (1/(Mi + Mn)) * (Mi*vy + Mn*vyn + Mi*g*Ry)
            vzn_ = (1/(Mi + Mn)) * (Mi*vz + Mn*vzn + Mi*g*Rz)

            if (well_random_number().le.0.5) then ! CX: (50% chance when spiraling)
              CX = .TRUE.
            end if 

          end if
        END IF

        IF (case.EQ.3) THEN ! collision for case C, M1 model
          chi = pi*(1.0 - b/d0)         
          cos_chi = cos(chi)
          sin_chi = abs(sin(chi))

          ! post-collision velocities
          vx_ =  vx +  Mn/(Mi + Mn) * (gx*(1-cos_chi) + hx*sin_chi)
          vy_ =  vy +  Mn/(Mi + Mn) * (gy*(1-cos_chi) + hy*sin_chi)
          vz_ =  vz +  Mn/(Mi + Mn) * (gz*(1-cos_chi) + hz*sin_chi)
          vxn_ = vxn - Mi/(Mi + Mn) * (gx*(1-cos_chi) + hx*sin_chi)
          vyn_ = vyn - Mi/(Mi + Mn) * (gy*(1-cos_chi) + hy*sin_chi)
          vzn_ = vzn - Mi/(Mi + Mn) * (gz*(1-cos_chi) + hz*sin_chi)

          if ((well_random_number().le.0.5).AND.(b.LE.bcx)) then  ! CX: identity switch (50% chance if b<bcx)
            CX = .TRUE.
          end if

        END IF 

        if (CX) then ! assign neutral post-collision velcoity to the ion (col + identity switch)
          ion(s)%part(i)%VX = vxn_/V_scale_ms
          ion(s)%part(i)%VY = vyn_/V_scale_ms
          ion(s)%part(i)%VZ = vzn_/V_scale_ms
        else ! assign ion post-collision velocity (no identity switch)
          ion(s)%part(i)%VX = vx_/V_scale_ms
          ion(s)%part(i)%VY = vy_/V_scale_ms
          ion(s)%part(i)%VZ = vz_/V_scale_ms
        end if

        ! Check for energy conservation. 
        E_ratio = (0.5*Mi*(vx_**2+vy_**2+vz_**2) + 0.5*Mn*(vxn_**2 + vyn_**2 + vzn_**2))/(0.5*Mi*(vx**2+vy**2+vz**2) + 0.5*Mn*(vxn**2 + vyn**2 + vzn**2))
        E_ratio_ion = (0.5*Mi*(vx_**2+vy_**2+vz_**2))/(0.5*Mi*(vx**2+vy**2+vz**2))
        if ((E_ratio-1.0).GT.1E-6) then
          PRINT *, "Error in PERFORM_ION_NEUTRAL_COLLISION. Energy not conserved. This should never happen. Energy lost (1-E_before/E_after): ", (E_ratio-1)
          PRINT *, "Collision case: ", case  
          PRINT *, "Charge exchange occured?: ", CX  
          CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

      end if ! if col
    END DO ! ion loop
  END DO ! ion species loop

END SUBROUTINE PERFORM_ION_NEUTRAL_COLLISION
