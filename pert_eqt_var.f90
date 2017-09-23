module pert_eqt_var
!Description:
!    This module is the interface containing the variables used in perturbation equations
!    It obtains the public variables from the module 'global_var' and transfer them to the perturbation equations

implicit none

public :: pevar

contains

subroutine pe_background(region, t, P_, rho_, m_, nu_, r_, drhodP_, OmeC_sq_, mu_) !retrieve background variables from global_var
use global_var
implicit none
character (len=*), intent(in) :: region
real(8), intent(in) :: t
real(8), intent(out) :: P_, rho_, m_, nu_, r_, drhodP_
complex(8), intent(out) :: OmeC_sq_
real(8), intent(out), optional :: mu_
integer :: ii
    if (region == 'core') then
		ii = nint((t - XCo_i)/dx_Co)*2
        P_ = P_Co(ii)
        rho_ = rho_Co(ii)
        m_ = m_Co(ii)
        nu_ = nu_Co(ii)
        r_ = r_Co(ii)
        drhodP_ = drhodP_Co(ii)
		if (present(mu_)) mu_ = mu_Co(ii)/2.d0 !Andersson 2015 mistake?
	elseif (region == 'crust') then
		ii = nint((t - XCr_i)/dx_Cr)*2
        P_ = P_Cr(ii)
        rho_ = rho_Cr(ii)
        m_ = m_Cr(ii)
        nu_ = nu_Cr(ii)
        r_ = r_Cr(ii)
		drhodP_ = drhodP_Cr(ii)
		if (present(mu_)) mu_ = mu_Cr(ii)/2.d0 !Andersson 2015 mistake?
	else
	    pause 'err: pe_background'
	endif
	OmeC_sq_ = OmeC_sq
endsubroutine pe_background

function pe_eLamb(rho_, m_, r_)
use parameters
implicit none
real(8) :: pe_eLamb, rho_, m_, r_
    pe_eLamb = (1.d0 - rel * 2.d0*m_/r_)**(-0.5d0)
endfunction pe_eLamb

function pe_dLamb(rho_, m_, r_)
use constants
use parameters
implicit none
real(8) :: pe_dLamb, rho_, m_, r_, c1
    c1 = pe_dLamb(rho_, m_, r_) **2
    pe_dLamb = c1 / r_ * rel * (4.d0*pi*r_**2*rho_ - m_/r_)
endfunction pes03_dLamb

function pe_dnu(P_, rho_, m_, r_)
use constants
use parameters
implicit none
real(8) :: pe_dnu, P_, rho_, m_, r_, c1
    c1 = pe_eLamb(rho_, m_, r_) **2
    pe_dnu = c1 / r_ * rel * ( rel * 4.d0*pi*r_**(2)*P_ + m_/r_)
endfunction pe_dnu

function pe_ddnu(P_, rho_, m_, r_)
use constants
implicit none
real(8) :: pe_ddnu, P_, rho_, m_, r_, c1
real(8) :: dLamb, dP, dnu
real(8) :: B(1:3)
    c1 = pe_eLamb(rho_, m_, r_) **2
    dLamb = pe_dLamb(rho_, m_, r_)
    dnu = pe_dnu(P_, rho_, m_, r_)
    dP = -dnu *(rho_ + P_)
		
    B(1) = 2.d0*dLamb*(4.d0*pi*r_*P_ + m_/r_**2)
    B(2) = 4.d0*pi*(P_+r_*dP)
    B(3) = 4.d0*pi*rho_ - 2.d0*m_/r_**3
    pe_ddnu = c1 * (B(1) +B(2) +B(3))
endfunction pe_ddnu
	
endmodule pert_eqt_var