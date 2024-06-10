
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

      if (well_random_number().le.probab_rcx) then
        ion(s)%part(i)%VX = vxn
        ion(s)%part(i)%VY = vyn
        ion(s)%part(i)%VZ = vzn
        collision_rcx(s)%counter = collision_rcx(s)%counter + 1
      end if
    END DO
  END DO   

END SUBROUTINE PERFORM_RESONANT_CHARGE_EXCHANGE
