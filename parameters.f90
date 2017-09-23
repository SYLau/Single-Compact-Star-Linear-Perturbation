module parameters

!Description:
!    This module contains the input parameters
!    The input parameters has the following options
	
	!shear_opt: 
	!   neutron_star_dh01
	!   neutron_star_cs_table
	!   shear_constant_cs
	!   neutron_star_Penner2011
	!   CCS_quark_star
	!   CCS_hybrid_star
	!   shear_constant_mu

implicit none
! 2a: grids
    real(8) :: delr_0 = 1.d0                                       !first step size for the variable grid in solving TOV eqt (cgs unit)
    real(8) :: delr_adj = 0.0115d0                                 !set a scale for the variable grid, i.e. variable grid dr \prop delr_adj
    real(8) :: delr_ref = delr_adj/1.d2                            !mesh refinement near the interface
	
	integer :: NMat = 4                                            !number of variables to match at the boundary for the perturbation problems (BVP/ eigen-value problem)
	
	character (len=30) :: eos_tab="input/eos_APR_ILQ.txt"
	character (len=30) :: comp_tab="input/comp_DH01_HP.txt"
	character (len=30) :: shear_tab="input/Rshear_Steiner_SLy4.txt"
	
	integer :: pg_N1 = 8000
	integer :: pg_N2 = 3000
	
! 2b: physical parameters/ options
	real(8) :: l_0 = 2.d0                                          !spherical harmonics degree
	real(8) :: rho_0 = 1.d15                                       !central density (cgs unit)
	real(8) :: P_t = 1.d30                                         !transition pressure (cgs unit)
	real(8) :: P_min = 1.d22                                       !cutoff pressure (cgs unit)

	character (len=30) :: sp_opt = '2layer'
	character (len=30) :: pg_opt = '2layer'
	character (len=30) :: pe_opt = 'fullgr'
	character (len=50) :: pef_opt = 'fullgr_solidcore_fluidcrust'
	character (len=50) :: eos_opt = 'eos_table'
	character (len=50) :: shear_opt = 'hybridstar'

	real(8) :: poly_K = 180.d0  ! K (geometrized unit, km)
	real(8) :: poly_n = 1.d0    ! polytropic index
	real(8) :: drho = 0.d0      ! density jump (fractional) in polytropic model
	real(8) :: rho_t = 9.d14    ! transition density in polytropic model (cgs unit)
	real(8) :: poly_cs = 2.d8   ! shear speed for the crust of polytropic star to mimic neutron star crust (cgs unit)
	real(8) :: QS_a4 = 0.7d0	! parameters of quark star models; (MeV unit); ref: PRD 89 103014 (2014)
	real(8) :: QS_a2 = 0.d0
	real(8) :: B_eff = 145.d0 ** 4
	real(8) :: QS_Gap = 5.d0
	real(8) :: QS_mu0 = 2.47d0	! unit MeV/fm^3; ref: PRD 89 103014 (2014)
	real(8) :: Ki_poly = 0.015d0	! shear modulus proportional contant Penner 2011
	real(8) :: mu_poly = 1.d-14	! shear modulus contant Penner 2011 (in km^-2)
	real(8) :: ic_Comp = 1.d-1	! one layer incompressible star
	real(8) :: ic_Comp_Co = 1.d-2	! two layers incompressible star core compactness
	real(8) :: Constant_mu = 1.d31 ! Constant Shear Modulus in cgs unit

	real(8) :: mu_red = 1.d0	! artificially reduce shear modulus
	real(8) :: rel = 1.d0	    ! adjust relativistic correction for relat cowling, range 0 < rel <= 1

	real(8) :: R_inf_f = 25.d0		! far-field limit of Zerilli eqt R_inf_f = r \omega (theoretically should be infinity)

endmodule parameters