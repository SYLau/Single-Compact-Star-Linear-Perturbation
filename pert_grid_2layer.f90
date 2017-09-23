module pert_grid_2layer

!Description:
!    This module contains the subroutines for reducing the no. of grids of 2-layer models
!    The subroutine 'pg2layer_body' contains the main flow of the mesh reduction
!    The grids are uniform in a new variable x = log(r/P) (McDermott et al. 1988)
!    Direct influence on the global_var
!

implicit none

public :: pg2layer_body
private :: pg2layer_xj, pg2layer_xi_Co,pg2layer_xi_Cr, pg2layer_ri_Ou

contains

subroutine pg2layer_body
implicit none
real(8) :: x_j(0:st_N2)
    call pg2layer_xj(x_j(0:st_N2)) !found below
    call pg2layer_xi_Co(x_j(0:st_N2)) !found below
    call pg2layer_xi_Cr(x_j(0:st_N2)) !found below
    call pg2layer_ri_Ou !found below
endsubroutine pg2layer_body

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

subroutine pg2layer_xj(x_j) !define a new variable x = log(r/P) (McDermott et al. 1988)
use global_var
implicit none
real(8) :: x_j(0:st_N2)
integer :: j
    do j = 0, st_N1
        x_j(j) = dlog(r(j)/P(j))
    enddo

    XCo_i = x_j(0)
    XCo_f = x_j(st_N1)

    do j = st_N1+1, st_N2
        x_j(j) = dlog(r(j)/P(j))
    enddo
    XCr_i = x_j(st_N1+1)
    XCr_f = x_j(st_N2)
endsubroutine pg2layer_xj

subroutine pg2layer_xi_Co(x_j)
use global_var
use interpolation
implicit none
real(8), dimension(0:st_N2) :: x_j
integer :: ii, iiN1
character(len=50) :: mode
    dx_Co = (XCo_f - XCo_i)/pg_N1
    r_Co(0) = r(0)
    x_Co(0) = XCo_i
    P_Co(0) = P(0)
    rho_Co(0) = rho(0)
    m_Co(0) = m(0)
    mu_Co(0) = mu(0)
    nu_Co(0) = nu(0)
	drhodP_Co(0) = drhodP(0)
	
	mode = '4pt_poly'
	iiN1 = pg_N1 *2
	do ii = 1, iiN1-1
		x_Co(ii) = dx_Co/2.d0 + x_Co(ii-1)
		r_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), r(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core radius
		P_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), P(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core pressure
		rho_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), rho(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core density
		m_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), m(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core mass
		mu_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), mu(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core shear modulus
		nu_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), nu(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core nu
	    drhodP_Co(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), drhodP(0:st_N2), x_Co(ii)) !found in module 'interpolation': interpolation; core drhodP
	enddo
	
    r_Co(pg_N1 *2) = r(st_N1)
    x_Co(pg_N1 *2) = XCo_f
    P_Co(pg_N1 *2) = P(st_N1)
    rho_Co(pg_N1 *2) = rho(st_N1)
    m_Co(pg_N1 *2) = m(st_N1)
    mu_Co(pg_N1 *2) = mu(st_N1)
    nu_Co(pg_N1 *2) = nu(st_N1)
	drhodP_Co(pg_N1 *2) = drhodP_Co(st_N1)
endsubroutine pg2layer_xi_Co
 
subroutine pg2layer_xi_Cr(x_j)
use global_var
use interpolation
implicit none
real(8), dimension(0:st_N2) :: x_j
integer :: ii, iiN2
character(len=50) :: mode
    dx_Cr = (XCr_f - XCr_i)/pg_N2
    r_Cr(0) = r(st_N1+1)
    x_Cr(0) = XCr_i
    P_Cr(0) = P(st_N1+1)
    rho_Cr(0) = rho(st_N1+1)
    m_Cr(0) = m(st_N1+1)
    mu_Cr(0) = mu(st_N1+1)
    nu_Cr(0) = nu(st_N1+1)
	drhodP_Cr(0) = drhodP(st_N1+1)

    mode = '4pt_poly'
    iiN2 = pg_N2*2
    do ii = 1, iiN2-1
        x_Cr(ii) = dx_Cr/2.d0 + x_Cr(ii-1)
        r_Cr(ii) = Int_api(mode, st_N1+1, st_N2, x_j(0:st_N2), r(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; crust radius
        P_Cr(ii) = Int_api(mode, st_N1+1, st_N2, x_j(0:st_N2), P(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; crust pressure
        rho_Cr(ii) = Int_api(mode, st_N1+1, st_N2, x_j(0:st_N2), rho(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; crust density
        m_Cr(ii) = Int_api(mode, st_N1+1, st_N2, x_j(0:st_N2), m(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; crust mass
        mu_Cr(ii) = Int_api(mode, st_N1+1, st_N2, x_j(0:st_N2), mu(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; crust shear modulus
        nu_Cr(ii) = Int_api(mode, st_N1+1, st_N2, x_j(0:st_N2), nu(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; crust nu
		drhodP_Cr(ii) = Int_api(mode, 0, st_N1, x_j(0:st_N2), drhodP(0:st_N2), x_Cr(ii)) !found in module 'interpolation': interpolation; core drhodP
    enddo
	
    r_Cr(pg_N2*2) = r(st_N2)
    x_Cr(pg_N2*2) = XCr_f
    P_Cr(pg_N2*2) = P(st_N2)
    rho_Cr(pg_N2*2) = rho(st_N2)
    m_Cr(pg_N2*2) = m(st_N2)
    mu_Cr(pg_N2*2) = mu(st_N2)
    nu_Cr(pg_N2*2) = nu(st_N2)
	drhodP_Cr(pg_N2*2) = drhodP(st_N2)
endsubroutine pg2layer_xi_Cr


subroutine pg2layer_ri_Ou
use global_var
implicit none
    P_Ou = P(st_N2)
    rho_Ou = rho(st_N2)
    m_Ou = m(st_N2)
    nu_Ou = nu(st_N2)
    ROu_i = r(st_N2)
endsubroutine pg2layer_ri_Ou

endmodule pert_grid_2layer