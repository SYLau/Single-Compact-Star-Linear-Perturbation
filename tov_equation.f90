module tov_equation

!Description:
!    This module contains the TOV equations and the boundary conditions for integration
!    The function 'tov_dpdr' contains the TOV equation dP/dr
!    The function 'tov_g' contains the TOV equation dnu/dr
!    The input/ output are in cgs units
!    Metric here is defined as e^(2nu) + ...
!
!    The module can be used separately in different unit systems by modifying 
!    the interface subroutine 'tov_constants' which defines the constants for the TOV equations
!

implicit none

public :: tov_dpdr, tov_g, tov_ic, tov_nu_BC
private :: tov_constants

contains

function tov_dpdr(P, rho, m, r)
REAL(8)	:: P, rho, m, r, G, c, pi, tov_dpdr
    CALL tov_constants(G, c, pi) !found below
    tov_dpdr = - (rho + P/c**2)* tov_g(P, rho, m, r)
endfunction tov_dpdr

function tov_g(P, rho, m, r)
REAL(8)	:: P, rho, m, r, G, c, pi, tov_g
    CALL tov_constants(G, c, pi) !found below
    tov_g = G * (m + 4.D0*pi*r**3*P/c**2)/ r**2/ (1.D0 - 2.D0*G*m/c**2/r)
endfunction tov_g

function tov_dmdr(P, rho, m, r)
REAL(8)	:: P, rho, m, r, G, c, pi, tov_dmdr
    CALL tov_constants(G, c, pi) !found below
    tov_dmdr = 4.d0*pi*r**2 * rho
endfunction tov_dmdr

subroutine tov_ic(P_0, rho_0, x)
implicit none
real(8) :: P_0, rho_0, x(1:4)
	x(1) = P_0
	x(2) = 0.d0
	x(3) = 0.d0
	x(4) = rho_0
endsubroutine tov_ic

subroutine tov_nu_BC(M0, R0, nu0, nu)
implicit none
real(8), intent(in) :: M0, R0, nu0
real(8), intent(inout) :: nu(0:)
real(8) :: Const, nuR, G, c, pi
	CALL tov_constants(G, c, pi) !found below
	nuR = 0.5d0 *dlog(1.d0 - 2.d0*G*M0/R0/c**2)
	Const = nuR - nu0
	nu = nu + Const
endsubroutine tov_nu_BC

! --------------------------------------------------------------------
!Interface
! --------------------------------------------------------------------

subroutine tov_constants(G_im, c_im, pi_im) !interface: import the values of constants from global variables
use constants
implicit none
REAL(8), intent(out) :: G_im, c_im, pi_im
    G_im = Grav_Const
    c_im = c
    pi_im = pi
endsubroutine tov_constants

endmodule tov_equation