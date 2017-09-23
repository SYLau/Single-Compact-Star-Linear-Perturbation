module static_profile_2layer

!Description:
!    This module contains the subroutines for calculating the static profile of 2-layer models
!    The subroutine 'st2layer_body' contains the main flow of the TOV equation integration
!    Metric here defined as e^(2nu) + ...
!
!    'st2layer_body' input are contained in the module 'parameters':
!        -rho_0 central density
!        -delr_0 first integration step for integration
!        -delr_adj scale for variable step size for integration
!        -delr_ref refined variable step size for integration near interface
!        -P_min cut off pressure
!    'st2layer_body' output:
!        -P,rho,m,nu,r,dr background profile
!        -st_N1,st_N2 position of interface and stellar surface in the array
!
!    Physical quantities of the array elements in numerical integration
!    x(1): pressure, x(2): mass, x(3): nu, x(4): rho
!

implicit none

public :: st2layer_body
private :: st_tov_eqts, st_refine, st_dr, st_display

contains

subroutine st2layer_body(P_,rho_,m_,nu_,r_,dr_,st_N1_,st_N2_, drhodP_, mu_)
use parameters !requires module 'parameters' (variables: )
use EOS
use shear_modulus
use tov_equation
use RK4_Set
implicit none
real(8), intent(out) :: P_(0:),rho_(0:),m_(0:),nu_(0:),r_(0:),dr_(0:), drhodP_(0:), mu_(0:)
integer, intent(out) :: st_N1_,st_N2_ !st_N1: last point of core; st_N1+1: first point of crust
real(8) :: x(1:4), y(1:4), P_0, dr_scale
integer :: i
    P_0 = eos_fcn('p(rho)',rho_0)
    if (P_0 <= P_t) pause,'err: P_0 <= P_t'
    call tov_ic(P_0, rho_0, x(1:4)) !found in module 'tov_equation': initial condition for TOV equations
	dr_(0) = delr_0
	dr_scale = delr_adj
    i = 0
    do while (P(i) > P_min) !numerical integration
        r_(i+1) = r_(i) + dr_(i)
		call RK4(st_tov_eqts, 3, r_(i), r_(i+1), x(1:3), y(1:3)) !found in module 'RK4_Set'
		y(4) = eos_fcn('rho(p)',y(1))
		call st_refine(x(1:4), y(1:4), i, delr_adj, delr_ref, dr_scale, st_N1_) !found below
		x = y
		i = i + 1
		P_(i) = x(1)
		m_(i) = x(2)
		nu_(i) = x(3)
		rho_(i) = x(4)
		dr_(i) = st_dr(x(1:4),dr_scale, r(i)) !found below
		drhodP_(i) = eos_fcn('drhodp(p)', P_(i))
	enddo

	call tov_nu_BC(m_(st_N2), r_(st_N2), nu_(st_N2), nu_(0:st_N2)) !found in module 'tov_equation': boundary condition for TOV equations
	
	call shear_calculate(P_(0:st_N2),rho_(0:st_N2),m_(0:st_N2),nu_(0:st_N2),r_(0:st_N2), st_N1_,st_N2_, drhodP_(0:st_N2), mu_(0:st_N2)) !found in module 'shear_modulus'
	
	call st_display(P_(0:st_N2),rho_(0:st_N2),m_(0:st_N2),nu_(0:st_N2),r_(0:st_N2), st_N1_,st_N2_, drhodP_(0:st_N2), mu_(0:st_N2)) !found below
   
endsubroutine st2layer_body

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

subroutine st_tov_eqts(n, t, xi, fcn) !The set of TOV equations
use EOS
use tov_equation
implicit none
integer :: n
real(8), dimension(n) :: xi, fcn
real(8) :: t, x4
    x4 = eos_fcn('rho(p)',xi(1))
    if (t == 0.d0) fcn(1) = 0.d0
    if (t /= 0.d0) fcn(1) = tov_dPdr(xi(1), x4, xi(2), t) !found in module 'tov_equation'
    fcn(2) = tov_dmdr(xi(1), x4, xi(2), t) !found in module 'tov_equation'
    if (t == 0.d0) fcn(3) = 0.d0
    if (t /= 0.d0) fcn(3) = tov_g(xi(1), x4, xi(2), t) !found in module 'tov_equation'
	endif
endsubroutine st_tov_eqts

subroutine st_refine(x, y, i, delr_adj, delr_ref, dr_, st_N1_) !Refine the step sizes, determine the position of the interface
implicit none
real(8), intent(in) :: x(1:4), y(1:4), delr_adj, delr_ref
real(8), intent(out) :: dr_
integer :: i, st_N1_
integer, save :: count_ = 1, refined_step
	if (P_t < x(1) .and. P_t >= y(1) .and. count_ == 1) then
		i = i - 1
		refined_step = i
		count_ = count_ + 1
		dr_ = delr_ref
		y = x
	elseif (P_t < x(1) .and. P_t >= y(1) .and. count_ == 2) then 
		st_N1_ = i
		refined_step = i - refined_step
		count_ = count_ + 1
	elseif (count_ > 2 .and. count_ < refined_step +2) then
		count_ = count_ + 1
	elseif (count_ == refined_step +2) then
		dr_ = delr_adj
		count_ = 1
	else
		return
	endif
endsubroutine st_refine

function st_dr(x, dr_, r_) !variable grid suitable for integrating TOV equations (Baym, Pethick, Sutherland 1971)
use tov_equation
implicit none
real(8), intent(in) :: x(1:4), dr_, r_
real(8) :: st_dr
    st_dr = dabs(dr_/(tov_dmdr(x(1), x(4), x(2), r_)/x(2) - tov_dPdr(x(1), x(4), x(2), r_)/x(1)))
endfunction st_dr

subroutine st_display(P_,rho_,m_,nu_,r_, st_N1_,st_N2_, drhodP_, mu_) !display the results on screen
use constants
implicit none
real(8), intent(in) :: P_(0:), m_(0:), nu_(0:), rho_(0:), r_(0:), drhodP_(0:), mu_(0:)
integer, intent(in) :: st_N1_, st_N2_
    write(*,*) "Initial Conditions:"
    write(*,*) "rho_0(/g cm^-3)", rho_(0)
    write(*,*) "P_0(/dym cm^-2)", P_(0)
    write(*,*) "=========================================================="
    write(*,*) "crust core transition:"
    write(*,*) "rho_-(/g cm^-3):", rho_(st_N1)
    write(*,*) "rho_+(/g cm^-3):", rho_(st_N1+1)
    write(*,*) "density jump(/rho_+):", (rho_(st_N1) - rho_(st_N1+1))/rho_(st_N1+1)
    write(*,*) "P_t(dyn cm^-2)", P_(st_N1)
    write(*,*) "mu_t(dyn cm^-2)", mu_(st_N1)
    write(*,*) "shear modulus jump(/P_t):", mu_(st_N1)/P(st_N1)
    write(*,*) "=========================================================="
    write(*,*) "Summary:"
    write(*,*) "Mmid(/solar mass) = ", m_(st_N1)/M_Solar
    write(*,*) "Rmid(/km) = ", r_(st_N1_)/1.d5
    write(*,*) "M(/solar mass) = ", m_(st_N2_)/M_Solar
    write(*,*) "R(/km) = ", r_(st_N2_)/1.d5
    write(*,*) "dMc(/solar mass) = ", (m_(st_N2_) - m_(st_N1_))/M_Solar
    write(*,*) "dRc(/km) = ", (r_(st_N2_) - r_(st_N1_))/ 1.d5
    write(*,*) "Compactness = ", Grav_const * m_(st_N2_)/c**2/r_(st_N2_)
    write(*,*) "P_cutoff(/dym cm^-2) = ", P_min
    write(*,*) "rho_cutoff(/g cm^-3) = ", rho_(st_N2_)
    write(*,*) "shear modulus (R) = ", mu_(st_N2)
    write(*,*) "shear speed^2 (R) = ", mu_(st_N2)/rho(st_N2)
    write(*,*) "Sqrt(M/R^3) = ", dsqrt(Grav_Const * m_(st_N2_)/r_(st_N2_)**3)
    write(*,*) "=========================================================="
    write(*,*) "LD formalism singularity freq(/Hz) = ", dsqrt( l_0*(l_0+1.d0)/2.d0 * Grav_Const *m_(st_N2_)/r_(st_N2_)**3 )/2.d0/pi
    write(*,*) "=========================================================="
endsubroutine st_display

endmodule static_profile_2layer