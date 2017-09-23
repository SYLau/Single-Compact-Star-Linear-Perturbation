module main

!Description:
!    This module contains the high-level subroutines of the program flow.
!    It branches out to call different subroutines with specific purposes
!    such as calculating the background stellar structures.
!
!Modifications:
!    Loops can be made in the subroutine 'main_body' to modify the program 
!    for scanning over certain parameters
!    e.g. scan along compactness by changing the central density rho_0

implicit none

public :: main_body
private :: main_initialization, main_stat_profile, main_pert_grid, main_pert_eqt

contains

subroutine main_body !main flow of the program body
implicit none
    call main_menu
	call main_initialization
    call main_stat_profile
    call main_pert_grid
    call main_pert_eqt
endsubroutine main_body

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

subroutine main_menu !load the settings and menu
use menu_01
implicit none
    call me01_body !found in module 'menu_01'
endsubroutine main_menu

subroutine main_initialization !initialize/ re-initialize parameters
use global_var
use parameters
use EOS
implicit none
    if (allocated(Coe)) deallocate(Coe)
    if (allocated(CoeC)) deallocate(CoeC)
	if (allocated(r_Co)) deallocate(r_Co,x_Co,P_Co,rho_Co,m_Co,mu_Co,nu_Co,drhodP_Co)
	if (allocated(r_Cr)) deallocate(r_Cr,x_Cr,P_Cr,rho_Cr,m_Cr,mu_Cr,nu_Cr,drhodP_Cr)
	
	if (eos_opt == 'eos_table') call eos_read_table('read') !read the EOS table
	
	allocate(r_Co(0:pg_N1*2), x_Co(0:pg_N1*2), P_Co(0:pg_N1*2), rho_Co(0:pg_N1*2), m_Co(0:pg_N1*2), mu_Co(0:pg_N1*2), nu_Co(0:pg_N1*2),drhodP_Co(0:pg_N1*2))
	allocate(r_Cr(0:pg_N2*2), x_Cr(0:pg_N2*2), P_Cr(0:pg_N2*2), rho_Cr(0:pg_N2*2), m_Cr(0:pg_N2*2), mu_Cr(0:pg_N2*2), nu_Cr(0:pg_N2*2),drhodP_Cr(0:pg_N2*2))
endsubroutine main_initialization

subroutine main_stat_profile !load the subroutines for calculating background stellar structures
use global_var !write the values of 'P,rho,m,nu,r,dr,st_N1,st_N2,mu' into the global variables
use parameters
use static_profile_2layer
implicit none
    if (sp_opt == '2layer') then
        call st2layer_body(P,rho,m,nu,r,dr,st_N1,st_N2,drhodP,mu)  !found in module 'static_profile_2layer'
    else
	    pause 'err: sp_opt'
	endif
endsubroutine main_stat_profile

subroutine main_pert_grid !load the subroutines for constructing the mesh for solving perturbation equations
use parameters
use pert_grid_2layer
implicit none
    if (pg_opt == '2layer') then
	    call pg2layer_body !found in module 'pert_grid_2layer'
	else
	    pause 'err: pg_opt'
	endif
endsubroutine main_pert_grid

subroutine main_pert_eqt !load the subroutines for solving perturbation equations
implicit none
endsubroutine main_pert_eqt

endmodule main