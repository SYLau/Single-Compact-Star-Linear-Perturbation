module global_var

!Description:
!    This module contains the public variables shared by different modules

implicit none

    logical :: eos_first_call

	integer, parameter :: st_Nmax = 200000  !upper limit of the number of grids in solving TOV eqt
	
	real(8), dimension(0:sp_Nmax) :: rho, P, m, mu, r, dr, nu, drhodP
	integer :: eos_n, comp_n, shear_n
	integer :: st_N1, st_N2


	character (len=30), dimension(1:NIface-1) :: pef_ziC, pef_ziD
	character (len=30), dimension(1:NIface-1) :: pef_LD_Cri1, pef_LD_Cri2, pef_LD_Cri3, pef_LD_Cri4, pef_LD_Cri5
	real(8), dimension(1:NIface-1) :: P_i
	integer :: sp_Ni(1:NIface-1)
	real(8), dimension(1:NIface-1, 0:pg_Ni*2) :: r_Cri, x_Cri, P_Cri, rho_Cri, m_Cri, mu_Cri, nu_Cri
	real(8), dimension(1:NIface-1) :: dx_Cri, XCri_i, XCri_f

	real(8), allocatable, dimension(:) :: r_Co, x_Co, P_Co, rho_Co, m_Co, mu_Co, nu_Co, drhodP_Co
	real(8), allocatable, dimension(:) :: r_Cr, x_Cr, P_Cr, rho_Cr, m_Cr, mu_Cr, nu_Cr, drhodP_Cr
	real(8), allocatable, dimension(:) :: r_Oc, x_Oc, P_Oc, rho_Oc, m_Oc, nu_Oc, drhodP_Oc
	real(8) :: P_Ou, rho_Ou, m_Ou, nu_Ou
	real(8) :: dx_Co, dx_Cr, dx_Oc, XCo_i, XCo_f, XCr_i, XCr_f, XOc_i, XOc_f
	real(8) :: dr_Ou, ROu_i, ROu_f
	
	real(8) :: afreqH, afreqL, afreq, Omega_sq
	complex(8) :: OmeC_sq
	logical :: pe_fcall

	real(8), allocatable :: Coe(:)
	complex(8), allocatable :: CoeC(:)
	integer :: gr_opt

	complex(8) :: k2
	logical :: zero_freq
	
endmodule global_var