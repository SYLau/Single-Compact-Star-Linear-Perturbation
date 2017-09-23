module shear_modulus

!Description:
!    This module contains the formula for calculating the shear modulus
!

implicit none

public :: shear_calculate
private :: shear_dh01_table, shear_nuclear_crust, shear_cs_table, shear_constant_cs, shear_Penner2011, shear_CCS, shear_constant_mu

contains

subroutine shear_calculate(P_,rho_,m_,nu_,r_, st_N1_,st_N2_,drhodP_,mu_)
use parameters
implicit none
real(8), intent(in) :: P_(0:),rho_(0:),m_(0:),nu_(0:),r_(0:), drhodP_(0:)
integer, intent(in) :: st_N1_,st_N2_
real(8), intent(inout) :: mu_(0:)
    if (shear_opt == 'neutron_star_dh01') then
	    call shear_dh01_table(rho_(0:st_N2_), st_N1_+1,st_N2_, mu_(0:st_N2_)) !found below
	elseif (shear_opt == 'neutron_star_cs_table') then
	    call shear_dh01_table(rho_(0:st_N2_), st_N1_+1,st_N2_, mu_(0:st_N2_)) !found below
    elseif (shear_opt == 'shear_constant_cs') then 
	    call shear_dh01_table(rho_(0:st_N2_), 0, st_N2_, mu_(0:st_N2_)) !found below
	elseif (shear_opt == 'neutron_star_Penner2011') then 
	    call shear_Penner2011(P_(0:st_N2_), 0, st_N2_, mu_(0:st_N2_)) !found below
	elseif (shear_opt == 'CCS_quark_star') then
	    call shear_CCS(rho_(0:st_N2_), 0,st_N2_, mu_(0:st_N2_)) !found below
	elseif (shear_opt == 'CCS_hybrid_star') then
	    call shear_CCS(rho_(0:st_N2_), 0,st_N1_, mu_(0:st_N2_)) !found below
	elseif (shear_opt == 'shear_constant_mu') then
	    call shear_constant_mu(0,st_N2_, mu_(0:st_N2_)) !found below
	else
	    pause 'shear_opt'
	endif
endsubroutine shear_calculate

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

subroutine shear_dh01_table(rho_, N1_,N2_, mu_)
use parameters
use interpolation
use FRead
implicit none
real(8), intent(in) :: rho_(0:)
integer, intent(in) :: N1_,N2_
real(8), intent(inout) :: mu_(0:)
integer :: j,comp_n
integer, parameter :: NMax_tab = 8000 !max number of rows in an input table
real(8), dimension(1:NMax_tab) :: n_in, nb_in, rho_in, nN_in, Z_in, A_in
real(8) :: nN_j, Z_j
    call RFile4R8(comp_tab, comp_n, n_in, rho_in, nN_in, Z_in)
	if (rho_(N1_) > max(rho_in(1), rho_in(comp_n))) pause 'err: shear_dh01_table out of range >'
	if (rho_(N2_) < min(rho_in(1), rho_in(comp_n))) pause 'err: shear_dh01_table out of range <'
    do j = N1_ , N2_
	    nN_j = Int_api('linear', 1, comp_n, rho_in, nN_in, rho_(j))
        Z_j = Int_api('linear', 1, comp_n, rho_in, Z_in, rho_(j))
        mu_(j) = shear_nuclear_crust(nN_j, Z_j) !found below
    enddo
endsubroutine shear_dh01_table

function shear_nuclear_crust(nN, Z)
implicit none
real(8) :: nN, Z, RCell, Ga, T
real(8) :: shear_nuclear_crust
    !Shear modulus formula used in McDermott et al. 1988:
    !shear_nuclear_crust = 0.3711d0*(Z*e)**2 * nN**(4.d0/3.d0)/2.d0**(1.d0/3.d0)
    
	!Shear modulus formula used in ~1990 onwards:
	T = 1.d8
    RCell = (3.d0/4.d0/pi/nN)**(1.d0/3.d0)
    Ga = (Z*e)**2/RCell/k_b/T
    shear_nuclear_crust = 0.1194d0*nN * (Z*e)**2.d0/RCell / (1.d0 + 0.595d0*(173.d0/Ga)**2)
endfunction shear_nuclear_crust

subroutine shear_cs_table(rho_, N1_,N2_, mu_)
use parameters
use interpolation
use FRead
implicit none
integer, parameter :: NMax_tab = 8000 !max number of rows in an input table
real(8), dimension(1:NMax_tab) :: r_, lgrho_, lgmu_
integer :: i, shear_n
real(8), intent(in) :: rho_(0:)
integer, intent(in) :: N1_,N2_
real(8), intent(inout) :: mu_(0:)
    call RFile2R8(shear_tab, shear_n, lgrho_, lgmu_)
	if (rho_(N1_) > max(10.d0**lgrho_(1), 10.d0**lgrho_(shear_n))) pause 'err: shear_dh01_table out of range >'
	if (rho_(N2_) < min(10.d0**lgrho_(1), 10.d0**lgrho_(shear_n))) pause 'err: shear_dh01_table out of range <'
    do i = N1_ , N2_
        mu_(i) = Int_api('linear', 1, shear_n, lgrho_, lgmu_, dlog10(rho_(i)))
        mu_(i) = 10.d0**(mu_(i))
    enddo
endsubroutine shear_cs_table

subroutine shear_constant_cs(rho_, N1_,N2_, mu_)
use parameters
implicit none
integer :: i
real(8), intent(in) :: rho_(0:)
integer, intent(in) :: N1_,N2_
real(8), intent(inout) :: mu_(0:)
    do i = N1_ , N2_
        mu_(i) = rho_(i)* poly_cs**2
    enddo
endsubroutine shear_constant_cs

subroutine shear_Penner2011(P_, N1_,N2_, mu_)
use parameters
implicit none
integer :: i
real(8), intent(in) :: P_(0:)
integer, intent(in) :: N1_,N2_
real(8), intent(inout) :: mu_(0:)
    do i = N1_ , N2_
        mu_(i) = P_(i)* Ki_poly + mu_poly* (1.d-10/(Grav_Const/c**4))
    enddo
endsubroutine shear_Penner2011

subroutine shear_CCS(rho_, N1_,N2_, mu_)
use constants
use parameters
implicit none
integer :: i
real(8) :: E_, C1_, A_, B_, C_, muq_2, muq_
real(8), intent(in) :: rho_(0:)
integer, intent(in) :: N1_,N2_
real(8), intent(inout) :: mu_(0:)
    !¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
    C1_ = 3.d0/4.d0/pi**(2.d0)
    A_ = C1_ * QS_a4
    B_ = - C1_ * QS_a2
    C_ = -B_eff
    do i = N1_ , N2_
        E_ = (h_bar*c)**3 * rho_(i)/((1.d13)*e_mks)**4 * c**2	! Energy density E_
        muq_2 = -B_/6.d0/A_ + dsqrt(B_**2 + 12.d0*A_*(C_+E_))/6.d0/A_
        mu_(i) = QS_mu0 *(QS_Gap/10.d0)**2 * muq_2/(400.d0)**2
        mu_(i) = mu_(i) * ((1.d13)*e_mks) * (1.d39)	! convert unit from MeV/fm^3 to cgs
    enddo
endsubroutine shear_CCS

subroutine shear_constant_mu(N1_,N2_, mu_)
use parameters
implicit none
integer, intent(in) :: N1_,N2_
real(8), intent(inout) :: mu_(0:)
integer :: i
    do i = N1_ , N2_
        mu_(i) = Constant_mu
    enddo
endsubroutine shear_constant_mu

endmodule shear_modulus