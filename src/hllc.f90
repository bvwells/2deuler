module solver_constants

   implicit none
!==============================================================================
   double precision, parameter :: TOL = 0.00000000001d0
   double precision, parameter :: gamma = 1.4d0
   double precision, parameter :: g1 = (gamma - 1.0d0)/(2.0d0*gamma)
   double precision, parameter :: g3 = 2.0d0*gamma/(gamma - 1.0d0)
   double precision, parameter :: g5 = 2.0d0/(gamma + 1.0d0)
   double precision, parameter :: g6 = (gamma - 1.0d0)/(gamma + 1.0d0)
   double precision, parameter :: g7 = 0.5d0*(gamma - 1.0d0)
!==============================================================================

end module solver_constants

subroutine hllc_ale(flux, cons_l, cons_r, x_dot)

   use solver_constants

   implicit none
!==============================================================================
   double precision, intent(INOUT), dimension(1:4) :: flux
   double precision, intent(IN), dimension(1:4) :: cons_l, cons_r
   double precision, intent(IN) :: x_dot
!==============================================================================
   double precision, dimension(1:4) :: U_star_left, U_star_right
   double precision :: rho_l, u_l, v_l, a_l, H_l, p_l, q_l, s_l, E_l
   double precision :: rho_r, u_r, v_r, a_r, H_r, p_r, q_r, s_r, E_r
   double precision :: p_star, s_star
!==============================================================================
   U_star_left = 0.0d0; U_star_right = 0.0d0
   rho_l = 0.0d0; u_l = 0.0d0; v_l = 0.0d0; a_l = 0.0d0; H_l = 0.0d0; p_l = 0.0d0; q_l = 0.0d0; s_l = 0.0d0; E_l = 0.0d0
   rho_r = 0.0d0; u_r = 0.0d0; v_r = 0.0d0; a_r = 0.0d0; H_r = 0.0d0; p_r = 0.0d0; q_r = 0.0d0; s_r = 0.0d0; E_r = 0.0d0
   p_star = 0.0d0; s_star = 0.0d0

! left state (primitive variables)
   rho_l = cons_l(1)
   u_l = cons_l(2)/cons_l(1)
   v_l = cons_l(3)/cons_l(1)
   p_l = (gamma - 1.0d0)*(cons_l(4) - 0.5d0*rho_l*(u_l**2 + v_l**2))

! right state (primitive variables)
   rho_r = cons_r(1)
   u_r = cons_r(2)/cons_r(1)
   v_r = cons_r(3)/cons_r(1)
   p_r = (gamma - 1.0d0)*(cons_r(4) - 0.5d0*rho_r*(u_r**2 + v_r**2))

! calculate local sound speeds
   a_l = sqrt(gamma*p_l/rho_l)
   a_r = sqrt(gamma*p_r/rho_r)

! Compute left middle and right signal velocities for the left and right discontinuities

   call initial_guess(p_star, S_star, u_l, u_r, rho_l, rho_r, p_l, p_r, a_l, a_r)

   H_l = p_star/p_l
   H_r = p_star/p_r

   if (H_l <= 1.0d0) then
      q_l = 1.0d0
   else
      q_l = sqrt(1.0d0 + ((gamma + 1.0d0)/(2.0d0*gamma))*(H_l - 1.0d0))
   endif

   if (H_r <= 1.0d0) then
      q_r = 1.0d0
   else
      q_r = sqrt(1.0d0 + ((gamma + 1.0d0)/(2.0d0*gamma))*(H_r - 1.0d0))
   endif

   S_l = u_l - a_l*q_l
   S_r = u_r + a_r*q_r

! Compute Numerical Flux

   E_l = p_l/(gamma - 1.0d0) + 0.5d0*rho_l*(u_l**2 + v_l**2)
   E_r = p_r/(gamma - 1.0d0) + 0.5d0*rho_r*(u_r**2 + v_r**2)

   if ((S_l - x_dot) >= 0.0d0) then

      flux(1) = (u_l - x_dot)*rho_l
      flux(2) = (u_l - x_dot)*rho_l*u_l + p_l
      flux(3) = (u_l - x_dot)*rho_l*v_l
      flux(4) = (u_l - x_dot)*E_l + u_l*p_l

   elseif (((S_l - x_dot) < 0.0d0) .and. ((S_star - x_dot) > 0.0d0)) then

      U_star_left(1) = rho_l*((S_l - u_l)/(S_l - S_star))
      U_star_left(2) = rho_l*((S_l - u_l)/(S_l - S_star))*S_star
      U_star_left(3) = rho_l*((S_l - u_l)/(S_l - S_star))*v_l
      U_star_left(4) = rho_l*((S_l - u_l)/(S_l - S_star))* &
                       ((E_l/rho_l) + (S_star - u_l)*(S_star + p_l/(rho_l*(S_l - u_l))))

      flux(1) = ((u_l - x_dot)*rho_l) + (S_l - x_dot)*(U_star_left(1) - rho_l)
      flux(2) = ((u_l - x_dot)*rho_l*u_l + p_l) + (S_l - x_dot)*(U_star_left(2) - rho_l*u_l)
      flux(3) = ((u_l - x_dot)*rho_l*v_l) + (S_l - x_dot)*(U_star_left(3) - rho_l*v_l)
      flux(4) = ((u_l - x_dot)*E_l + u_l*p_l) + (S_l - x_dot)*(U_star_left(4) - E_l)

   elseif (abs(S_star - x_dot) .lt. 1d-10) then

      U_star_left(1) = rho_l*((S_l - u_l)/(S_l - S_star))
      U_star_left(2) = rho_l*((S_l - u_l)/(S_l - S_star))*S_star
      U_star_left(3) = rho_l*((S_l - u_l)/(S_l - S_star))*v_l
      U_star_left(4) = rho_l*((S_l - u_l)/(S_l - S_star))* &
                       ((E_l/rho_l) + (S_star - u_l)*(S_star + p_l/(rho_l*(S_l - u_l))))

      flux(1) = ((u_l - x_dot)*rho_l) + (S_l - x_dot)*(U_star_left(1) - rho_l)
      flux(2) = ((u_l - x_dot)*rho_l*u_l + p_l) + (S_l - x_dot)*(U_star_left(2) - rho_l*u_l)
      flux(3) = ((u_l - x_dot)*rho_l*v_l) + (S_l - x_dot)*(U_star_left(3) - rho_l*v_l)
      flux(4) = ((u_l - x_dot)*E_l + u_l*p_l) + (S_l - x_dot)*(U_star_left(4) - E_l)

      U_star_right(1) = rho_r*((S_r - u_r)/(S_r - S_star))
      U_star_right(2) = rho_r*((S_r - u_r)/(S_r - S_star))*S_star
      U_star_right(3) = rho_r*((S_r - u_r)/(S_r - S_star))*v_r
      U_star_right(4) = rho_r*((S_r - u_r)/(S_r - S_star))* &
                        ((E_r/rho_r) + (S_star - u_r)*(S_star + p_r/(rho_r*(S_r - u_r))))

      flux(1) = 0.5d0*(flux(1) + ((u_r - x_dot)*rho_r) + (S_r - x_dot)*(U_star_right(1) - rho_r))
      flux(2) = 0.5d0*(flux(2) + ((u_r - x_dot)*rho_r*u_r + p_r) + (S_r - x_dot)*(U_star_right(2) - rho_r*u_r))
      flux(3) = 0.5d0*(flux(3) + ((u_r - x_dot)*rho_r*v_r) + (S_r - x_dot)*(U_star_right(3) - rho_r*v_r))
      flux(4) = 0.5d0*(flux(4) + ((u_r - x_dot)*E_r + u_r*p_r) + (S_r - x_dot)*(U_star_right(4) - E_r))

   elseif (((S_star - x_dot) < 0.0d0) .and. ((S_r - x_dot) > 0.0d0)) then

      U_star_right(1) = rho_r*((S_r - u_r)/(S_r - S_star))
      U_star_right(2) = rho_r*((S_r - u_r)/(S_r - S_star))*S_star
      U_star_right(3) = rho_r*((S_r - u_r)/(S_r - S_star))*v_r
      U_star_right(4) = rho_r*((S_r - u_r)/(S_r - S_star))* &
                        ((E_r/rho_r) + (S_star - u_r)*(S_star + p_r/(rho_r*(S_r - u_r))))

      flux(1) = ((u_r - x_dot)*rho_r) + (S_r - x_dot)*(U_star_right(1) - rho_r)
      flux(2) = ((u_r - x_dot)*rho_r*u_r + p_r) + (S_r - x_dot)*(U_star_right(2) - rho_r*u_r)
      flux(3) = ((u_r - x_dot)*rho_r*v_r) + (S_r - x_dot)*(U_star_right(3) - rho_r*v_r)
      flux(4) = ((u_r - x_dot)*E_r + u_r*p_r) + (S_r - x_dot)*(U_star_right(4) - E_r)

   else

      flux(1) = (u_r - x_dot)*rho_r
      flux(2) = (u_r - x_dot)*rho_r*u_r + p_r
      flux(3) = (u_r - x_dot)*rho_r*v_r
      flux(4) = (u_r - x_dot)*E_r + u_r*p_r

   endif

   return

end subroutine hllc_ale

subroutine initial_guess(p_star, u_star, u_l, u_r, rho_l, rho_r, p_l, p_r, a_l, a_r)

! Hybrid initial guess using primitive variables, two-rarefraction and two-shock approximations.

   use solver_constants

   implicit none
!==============================================================================
   double precision, intent(IN) :: u_l, u_r, rho_l, rho_r, p_l, p_r, a_l, a_r
   double precision, intent(INOUT) :: p_star, u_star
!==============================================================================
   double precision :: ltol, p_pv, p_min, p_max, g_l, g_r, P_LR, z, a_bar, rho_bar
!==============================================================================

   ltol = 0.000001

   p_pv = 0.5d0*(p_l + p_r) - 0.125d0*(u_r - u_l)*(rho_l + rho_r)*(a_l + a_r)
   p_min = MIN(p_l, p_r)
   p_max = MAX(p_l, p_r)

   if (((p_max/p_min) <= 2.0d0) .and. ((p_min <= p_pv) .and. (p_pv <= p_max))) then

!  Use primitive variables as guess

      a_bar = (a_l + a_r)/2.0d0
      rho_bar = (rho_l + rho_r)/2.0d0
      p_star = MAX(tol, p_pv)
      u_star = (u_l + u_r)/2.0d0 - (p_r - p_l)/(2.0d0*a_bar*rho_bar)

   else
      if (p_pv < p_min) then

!  Use two-rarefraction solution

         p_star = ((a_l + a_r - g7*(u_r - u_l))/(a_l/(p_l**g1) + a_r/(p_r**g1)))**g3
         z = (gamma - 1.0d0)/(2.0d0*gamma)
         P_LR = (p_l/p_r)**z
         u_star = (P_LR*u_l/a_l + u_r/a_r + 2.0d0*(P_LR - 1.0d0)/(gamma - 1.0d0))/(P_LR/a_l + 1.0d0/a_r)

      else

!  Use two-shock approximation with p_pv as estimate

         g_l = sqrt((g5/rho_l)/(g6*p_l + MAX(ltol, p_pv)))
         g_r = sqrt((g5/rho_r)/(g6*p_r + MAX(ltol, p_pv)))
         p_star = (g_l*p_l + g_r*p_r - (u_r - u_l))/(g_l + g_r)
         p_star = MAX(ltol, p_star)
         u_star = (u_l + u_r)/2.0d0 + ((p_star - p_r)*g_r - (p_star - p_l)*g_l)/2.0d0

      endif
   endif

   return

end subroutine initial_guess
