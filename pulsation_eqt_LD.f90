module pulsation_eqt_LD
!Description:
!    This module contains the subroutines for calculating quasi-normal modes
!

implicit none

public :: peLD_fluid, peLD_solid, pe_region
private :: peLD_E_Fluid, peLD_E_Solid

contains

subroutine peLD_fluid(n, t, x, fcn)
use constants
use parameters
use pert_eqt_var
implicit none
integer :: n
real(8), dimension(1:n) :: x, fcn
real(8) :: t
complex(8), dimension(1:n/2,1:n/2) :: B
complex(8), dimension(1:2,1:n/2) :: E
complex(8), dimension(1:n/2) :: xc, fcnc, H0_C, V_C			!	H0_C, V_C	-	Coefficients in the basis {H1, K, W, X}
complex(8) :: OmeC_sq
real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, r_
real(8) :: pi4, l_1, V1, ddNu
integer :: i, j
character (len=30) :: region
    if (mod(n,2) /= 0) pause 'err: puls eqt set n not even'
    xc(1:n/2) = dcmplx(x(1 :n/2),x(n/2+1 :n))
    
    call pe_region('retrieve', region)
    call pe_background(region, t, P_, rho_, m_, nu_, r_, drhodP_, OmeC_sq)
    gamma_1 = (1.d0 + rho_/P_)/ drhodP_
	
    call peLD_E_Fluid(P_, rho_, m_, nu_, r_, OmeC_sq, E)

    pi4 = 4.d0*pi
    l_1 = l_0 * (l_0 + 1.d0)
    eLamb = pe_eLamb(rho_, m_, r_) ** 2
    eNu = dexp(2.d0*nu_)
    eLamb_hf = eLamb**(0.5d0)
    eNu_hf = eNu**(0.5d0)
    dLamb = pe_dLamb(rho_, m_, r_) * 2.d0
    dNu = pe_dnu(P_, rho_, m_, r_) * 2.d0
    ddNu = pe_ddnu(P_, rho_, m_, r_) * 2.d0
    V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0
    H0_C(1:4) = E(1,1:4)
    V_C(1:4) = E(2,1:4)

    B(1,1) = -(l_0 + 1.d0 + 2.d0*m_*eLamb/r_ + pi4 * r_**2 * eLamb* (P_ - rho_))
    B(1,2) = eLamb
    B(1,3) = 0.d0
    B(1,4) = 0.d0
    B(1, 1:4) = B(1, 1:4) + eLamb*(H0_C(1:4) - 4.d0*pi4 * (P_ + rho_) * V_C(1:4) )	! Adding the algebraic solution of H0 & V

    B(2,1) = l_1/2.d0
    B(2,2) = -((l_0 + 1.d0) - dNu/2.d0*r_)
    B(2,3) = -2.d0*pi4 * eLamb_hf * (P_ + rho_) 
    B(2,4) = 0.d0
    B(2, 1:4) = B(2, 1:4) + H0_C(1:4)

    B(3,1) = 0.d0
    B(3,2) = r_**2 * eLamb_hf
    B(3,3) = - (l_0 + 1.d0)
    B(3,4) = r_**2 * eLamb_hf *(eNu_hf**(-1.d0)/gamma_1/P_)
    B(3, 1:4) = B(3, 1:4) + r_**2 * eLamb_hf *( 0.5d0 *H0_C(1:4) +  (-l_1/r_**2) *V_C(1:4))

    B(4,1) = (P_ + rho_)*eNu_hf * r_* (0.5d0*((r_ *OmeC_sq * eNu**(-1.d0)) + 0.5d0 * l_1/r_))
    B(4,2) = (P_ + rho_)*eNu_hf * r_* (0.5d0*(3.d0/2.d0 * dNu - 1.d0/r_) )
    B(4,3) = (P_ + rho_)*eNu_hf * r_* (-1.d0/r_*(pi4*(P_ + rho_)* eLamb_hf + OmeC_sq * eLamb_hf/eNu - r_**2/2.d0* eLamb_hf**(-1.d0)/r_**2 *(ddNu - dLamb*dNu/2.d0 - 2.d0/r_*dNu) ) )
    B(4,4) = -l_0
    B(4, 1:4) = B(4, 1:4) + (P_ + rho_)*eNu_hf *r_* ( 0.5d0 * (1.d0/r_ - dNu/2.d0)*H0_C(1:4) - l_1/2.d0/r_**2 * dNu *V_C(1:4))

    do i = 1,n/2
        fcnc(i) = 0.d0
        do j = 1,n/2
            fcnc(i) = fcnc(i) + B(i,j)*xc(j)/(1.d0 +V1)
        enddo
    enddo
    fcn(1:n/2) = dreal(fcnc(1:n/2))
    fcn(n/2+1 :n) = dimag(fcnc(1:n/2))
endsubroutine peLD_fluid

subroutine peLD_solid(n, t, x, fcn)
use constants
use parameters
use pert_eqt_var
implicit none
integer :: n
real(8), dimension(1:n) :: x, fcn
real(8) :: t
complex(8), dimension(1:n/2,1:n/2) :: B
complex(8), dimension(1:3,1:n/2) :: E
complex(8), dimension(1:n/2) :: xc, fcnc, H2_C, X_C, T1_C
complex(8) :: OmeC_sq
real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, mu_, r_
real(8) :: pi4, l_1, l_2_hf, V1, ddNu, dP
integer :: i, j
    if (mod(n,2) /= 0) pause 'err: puls eqt set n not even'
    xc(1:n/2) = dcmplx(x(1 :n/2),x(n/2+1 :n))

    call pe_region('retrieve', region)
    call pe_background(region, t, P_, rho_, m_, nu_, r_, drhodP_, OmeC_sq, mu_)
    gamma_1 = (1.d0 + rho_/P_)/ drhodP_
    call peLD_E_Solid(P_, rho_, m_, nu_, mu_, gamma_1, r_, OmeC_sq, E)

    pi4 = 4.d0*pi
    l_1 = l_0 * (l_0 + 1.d0)
    l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0
    eLamb = pe_eLamb(rho_, m_, r_) ** 2
    eNu = dexp(2.d0*nu_)
    eLamb_hf = eLamb**(0.5d0)
    eNu_hf = eNu**(0.5d0)
    dLamb = pe_dLamb(rho_, m_, r_) * 2.d0
    dNu = pe_dnu(P_, rho_, m_, r_) * 2.d0
    ddNu = pe_ddnu(P_, rho_, m_, r_) * 2.d0
    dP = -dNu/2.d0 * (P_ + rho_)
    V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0
    H2_C(1:6) = E(1,1:6)
    X_C(1:6) = E(2,1:6)
    T1_C(1:6) = E(3,1:6)

    B = 0.d0

    B(1,1) = 0.5d0*(dLamb - dNu)*r_ - (l_0 + 1.d0)
    B(1,2) = eLamb
    B(1,5) = -eLamb*4.d0*pi4*(P_ + rho_)
    B(1, 1:6) = B(1, 1:6) + eLamb*H2_C(1:6)	! Adding the algebraic solution of H2 & X & T1

    B(2,1) = l_2_hf + 1.d0
    B(2,2) = 0.5d0* dNu*r_ - (l_0 + 1.d0)
    B(2,4) = -2.d0*pi4*(P_ + rho_)*eLamb_hf
    B(2, 1:6) = B(2, 1:6) + H2_C(1:6)

    B(3,1) = -r_**2 * eNu**(-1.d0) * OmeC_sq
    B(3,2) = l_0
    B(3,3) = -(0.5d0* dNu*r_ + (l_0 - 1.d0))
    B(3,6) = -4.d0*pi4
    B(3, 1:6) = B(3, 1:6) + B(2, 1:6) - (0.5d0* dNu*r_ + 1.d0)*H2_C(1:6)

    B(4,2) = r_**2 * eLamb_hf
    B(4,4) = -(l_0 + 1.d0)
    B(4,5) = -eLamb_hf * l_1
    B(4, 1:6) = B(4, 1:6) + 0.5d0 * r_**2 * eLamb_hf * H2_C(1:6) +  r_**2 * eLamb_hf * eNu_hf**(-1.d0)/ gamma_1 / P_ * X_C(1:6)

    B(5,4) = eLamb_hf
    B(5,5) = (2.d0 - l_0)
    B(5,6) = 0.5d0/mu_

    B(6,3) = -0.5d0* r_**2 *eLamb * (P_ + rho_)
    B(6,4) = eLamb_hf * dP * r_
    B(6,5) = 4.d0 * l_2_hf * eLamb * mu_ - eLamb/eNu * r_**2 * OmeC_sq * (P_ + rho_)
    B(6,6) = 0.5d0 * (dLamb - dNu) * r_ - (l_0 + 1.d0)
    B(6, 1:6) = B(6, 1:6) + r_**2 * eLamb/eNu_hf * (X_C(1:6) - 0.5d0/r_**2 * eNu_hf * T1_C(1:6))

    do i = 1,n/2
        fcnc(i) = 0.d0
        do j = 1,n/2
            fcnc(i) = fcnc(i) + B(i,j)*xc(j)/(1.d0 +V1)
        enddo
    enddo
    fcn(1:n/2) = dreal(fcnc(1:n/2))
    fcn(n/2+1 :n) = dimag(fcnc(1:n/2))
endsubroutine peLD_solid

! --------------------------------------------------------------------
!Interface
! --------------------------------------------------------------------

subroutine pe_region(mode, region)
implicit none
character (len=*), intent(in) :: mode
character (len=30), save :: region_store
character (len=*), intent(inout) :: region
    if (mode = 'read') then
	    region_store = region
	elseif (mode ='retrieve') then
	    region = region_store
	else
	    pause 'err: pe_region'
	endif
endsubroutine pe_region

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

subroutine peLD_E_Fluid(P_, rho_, m_, nu_, r_, OmeC_sq, E_)
use constants
use parameters
use pert_eqt_var
implicit none
real(8) :: P_, rho_, m_, nu_, r_, l_1, l_2, pi4, dP
real(8) :: eLamb, eNu
complex(8) :: Alpha_(1:3), E_(1:2, 1:4), OmeC_sq
    pi4 = 4.d0 * pi
    l_1 = l_0 * (l_0 + 1.d0)
    l_2 = (l_0 + 2.d0) * (l_0 - 1.d0)
    eLamb = pe_eLamb(rho_, m_, r_)**2
    eNu = dexp(nu_)**2
    dP = -(rho_ + P_)*pe_dnu(P_, rho_, m_, r_)
    Alpha_(1) = (3.d0 * m_ + 0.5d0 * l_2 * r_ + 4.d0 * pi * r_**3 * P_)	! Coefficient of H0 in eqt (1)
    Alpha_(2) = OmeC_sq * (P_ + rho_) * eNu **  (-0.5d0)	! Coefficient of V in eqt (2)
    Alpha_(3) = 0.5d0 * (P_ + rho_) * eNu **  (0.5d0)		! Coefficient of H0 in eqt (2)

    E_(1,1) = -1.d0/Alpha_(1) * ( 0.5d0 * l_1 * (m_ +pi4*r_**3 * P_) - OmeC_sq*r_**3 *(eLamb*eNu)**(-1.d0))
    E_(1,2) = 1.d0/Alpha_(1) * (0.5d0 * l_2 * r_ - OmeC_sq*r_**3 *(eNu)**(-1.d0) - eLamb/r_ * (m_ + pi4*r_**3 * P_) * (3.d0*m_ - r_ + pi4*r_**3 * P_))
    E_(1,3) = (0.d0, 0.d0)
    E_(1,4) = 1.d0/Alpha_(1) * 2.d0*pi4*r_**3 *(eNu)**(-0.5d0)
    if (OmeC_sq /= (0.d0, 0.d0)) then
        E_(2,1) = -Alpha_(3)/Alpha_(2) * E_(1,1)
        E_(2,2) = -Alpha_(3)/Alpha_(2) * E_(1,2)
        E_(2,3) = 1.d0/Alpha_(2) * dP/r_ * (eNu/eLamb) ** (0.5d0)
        E_(2,4) = 1.d0/Alpha_(2) - Alpha_(3)/Alpha_(2)*E_(1,4)
    else
        pause 'err: peLD_E_Fluid, frequency 0'
    endif
endsubroutine peLD_E_Fluid

subroutine peLD_E_Solid(P_, rho_, m_, Nu_, mu_, gamma_1, r_, OmeC_sq, E_)
use constants
use parameters
use pert_eqt_var
implicit none
real(8) :: P_, rho_, m_, Nu_, r_, l_1, l_2_hf, pi4, dP, Q, dNu, mu_, gamma_1
real(8) :: eLamb, eNu
complex(8) :: A_(1:4), B_(1:2, 1:6), E_(1:3, 1:6), OmeC_sq
integer :: i
		pi4 = 4.d0 * pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = 0.5d0 * (l_0 + 2.d0) * (l_0 - 1.d0)
		eLamb = pe_eLamb(rho_, m_, r_)**2
		eNu = dexp(Nu_)**2
		dNu = pe_dnu(P_, rho_, m_, r_)*2.d0
		dP = -(rho_ + P_)*pe_dnu(P_, rho_, m_, r_)
		Q = 0.5d0 * r_**2 * eLamb**(-1.d0) * dNu

		A_(1) = 2.d0*pi4*r_**3*eNu**(-0.5d0)			! Coefficient of X in eqt (2)
		A_(2) = 2.d0*pi4*r_								! Coefficient of T1 in eqt (2)
		A_(3) = 2.d0/3.d0*eNu**(-0.5d0)*mu_ *r_**2		! Coefficient of X in eqt (3)
		A_(4) = -0.25d0*gamma_1*P_						! Coefficient of T1 in eqt (3)

		B_ = 0.d0
		B_(1,1) = (l_2_hf + 1.d0)*Q - OmeC_sq*r_**3 * (eLamb*eNu)**(-1.d0)						! Coefficient of H1 in eqt (2)
		B_(1,2) = -(l_2_hf*r_ - OmeC_sq*r_**3 *eNu**(-1.d0) - eLamb/r_*Q *(2.d0*m_ + Q - r_))	! Coefficient of K in eqt (2)
		B_(1,3) = 2.d0*m_ + Q + l_2_hf*r_														! Coefficient of H0 in eqt (2)
		B_(1,6) = 4.d0*pi4*r_ * eLamb**(-1.d0)													! Coefficient of T2 in eqt (2)

		B_(2,2) = - mu_ * gamma_1 * P_ * r_**2													! Coefficient of K in eqt (3)
		B_(2,4) = 2.d0 * mu_ * gamma_1 * P_ * eLamb**(-0.5d0)									! Coefficient of W in eqt (3)
		B_(2,5) = l_1 * mu_ * gamma_1 * P_														! Coefficient of V in eqt (3)

		E_ = 0.d0																				! H2
		E_(1,3) = 1.d0
		E_(1,5) = 16.d0*pi4*mu_
		
		E_(2,1:6) = (B_(1,1:6)/A_(2) - B_(2,1:6)/A_(4))/(A_(1)/A_(2) - A_(3)/A_(4))				! X		! Crammer's Rule
		E_(3,1:6) = (B_(1,1:6)/A_(1) - B_(2,1:6)/A_(3))/(A_(2)/A_(1) - A_(4)/A_(3))				! T1	! Crammer's Rule
endsubroutine peLD_E_Solid

endmodule pulsation_eqt_LD