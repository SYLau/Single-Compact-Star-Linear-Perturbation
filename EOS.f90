module EOS

!Description:
!    This module contains the equation of state of the stellar models
!    Most of the subroutines require using the variables in the module 'global_var'/ 'parameters'/ 'constants'
!
!Input:
!f: the options for eos_fcn
!		-p(rho): output P as a function of rho
!		-rho(p): output rho as a function of P
!		-drhodp(p): output drho/ dP as a function of P
!x: the independent variable for the EOS function,
!can either be rho or P depending on the option f

implicit none

public :: eos_fcn, eos_read_table
private :: eos_eos_table, eos_polyFerrari2003, eos_quark_star_Alford2005, eos_incompressible_star

contains

function eos_fcn(f, x_) ! in rho or P in cgs unit, return P or rho in cgs unit
use parameters !requires module 'parameters' (variables: eos_opt)
implicit none
real(8) :: eos_fcn, x_
character(len=*), optional :: f
    if (eos_opt == 'eos_table') then
        eos_fcn = eos_eos_table(f, x_) !found below
	elseif (eos_opt == 'polyFerrari2003') then
	    eos_fcn = eos_polyFerrari2003(f, x_) !found below
    elseif (eos_opt == 'quark_star') then 
        eos_fcn = eos_quark_matter_Alford2005(f, x_) !found below
	elseif (eos_opt == 'incompressible_star') then
	    eos_fcn = eos_incompressible_star(f, x_) !found below
	elseif (eos_opt == '2layer_incompressible_star') then
	    eos_fcn = eos_2layer_inc_star(f, x_) !found below	
	else
		pause 'err: eos function'
	endif
endfunction eos_fcn

subroutine eos_read_table(mode, eos_n_in, n_in, nb_in, rho_in, P_in)
use parameters
use FRead
implicit none
integer, parameter :: NMax_tab = 1000 !max number of rows in an input table
real(8), dimension(1:NMax_tab), optional, intent(out) :: n_, nb_, rho_, P_
integer, optional, intent(out) :: eos_n
real(8), dimension(1:NMax_tab), save :: n_save, nb_save, rho_save, P_save
integer, save :: eos_n_save
character (len=30) :: mode
    if (mode == 'read') then
        if (present(eos_n_in) == .true.) pause 'err: eos_read_table 2'
		eos_n_save = 0
		n_save = 0.d0
		nb_save = 0.d0
		rho_save = 0.d0
		P_save = 0.d0
	    call RFile4R8(eos_tab,eos_n_save, n_save, nb_save, rho_save, P_save)
	elseif (mode == 'retrieve') then
	    if (present(eos_n_in) /= .true. .or. eos_n_save == 0) pause 'err: eos_read_table 3'
        eos_n_in = eos_n_save
		n_in = n_save
		nb_in = nb_save
		rho_in = rho_save
		P_in = P_save
	else
	    pause 'err: eos_read_table 1'
    endif
endsubroutine eos_read_table

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

function eos_eos_table(f, x_) !obtain EOS from EOS table
use parameters
use Interpolation
implicit none
integer, parameter :: NMax_tab = 1000 !max number of rows in an input table
real(8), dimension(1:NMax_tab) :: n_, nb_, rho_, P_
integer :: eos_n
real(8) :: eos_eos_table, x_
character (len=*) :: f
    call eos_read_table('retrieve', eos_n_in=eos_n, n_in=n_, nb_in=nb_, rho_in=rho_, P_in=P_)
write(*,*) 'EOS, eos_eos_table', P_(1), P_(2)
pause
    if (f == 'p(rho)') then
        eos_eos_table = Int_api('linear_Ln', 1, eos_n, rho_, P_, x_) !interpolate p(rho)
    elseif (f == 'rho(p)') then
        eos_eos_table = Int_api('linear_Ln', 1, eos_n, P_, rho_, x_) !interpolate rho(p)
    elseif (f == 'drhodp(p)') then
        eos_eos_table = Int_api('linear_Ln_slope', 1, eos_n, P_, rho_, x_) !interpolate drhodp(p)
    else
        pause 'err: static profile eos table'
    endif
endfunction eos_eos_table

function eos_polyFerrari2003(f, x_) !polytropic EOS from the polytropic model in Ferrari 2003
use parameters !requires module 'parameters' (variables: )
implicit none
real(8) :: eos_polyFerrari2003, x_
character (len=*) :: f
real(8) :: P_, rho_, P_t_, rho_t_, K_, gamma_
	rho_t_ = Grav_Const/c**2 * rho_t
	K_ = poly_K*(1.d5)**(2.d0/poly_n)/(1.d0 + drho)**((poly_n + 1.d0)/poly_n)
    P_t_ = K_*(rho_t*(1.d0 + drho))**((poly_n + 1.d0)/poly_n)
	P_t = P_t_/(Grav_Const/c**4) !modify the transition pressure

    if (drho == 0.d0) then !drho = 0
		if (f == 'p(rho)') then
            rho_ = Grav_Const/c**2 * x_
        	P_ = K_ * rho_**((poly_n + 1.d0)/poly_n)
	        eos_polyFerrari2003 = P_ / (Grav_Const/c**4)
        elseif (f == 'rho(p)') then
            P_ = Grav_Const/c**4 * x_
            rho_ = (P_ / K_ )**(poly_n/(poly_n + 1.d0))
	        eos_polyFerrari2003 = rho_ / (Grav_Const/c**2)
        elseif (f == 'drhodp(p)') then
            P_ = Grav_Const/c**4 * x_
            rho_ = (P_ / K_ )**(poly_n/(poly_n + 1.d0))
			eos_polyFerrari2003 = 1.d0/((poly_n + 1.d0)/(poly_n) * (P_/rho_))
        else
            pause 'err: static profile eos table'
        endif
    else !drho /= 0
	    if (f == 'p(rho)') then
            rho_ = Grav_Const/c**2 * x_
		    if (rho_ >= rho_t*(1.d0+drho)) P_ = K_ * rho_**((poly_n + 1.d0)/poly_n)
			if (rho_ < rho_t) P_ = K_ * (rho_*(1.d0 + drho))**((poly_n + 1.d0)/poly_n)
			if (rho_ < rho_t*(1.d0+drho) .and. rho_ >= rho_t) pause 'err: eos_polyFerrari2003'
	        eos_polyFerrari2003 = P_ / (Grav_Const/c**4)
        elseif (f == 'rho(p)') then
            P_ = Grav_Const/c**4 * x_
			if (P_ >= P_t) rho_ = (P_ / K_ )**(poly_n/(poly_n + 1.d0))
			if (P_ < P_t) rho_ = (P_ / K_ )**(poly_n/(poly_n + 1.d0))/(1.d0 + drho)
	        eos_polyFerrari2003 = rho_ / (Grav_Const/c**2)
        elseif (f == 'drhodp(p)') then
            P_ = Grav_Const/c**4 * x_
            if (P_ >= P_t) rho_ = (P_ / K_ )**(poly_n/(poly_n + 1.d0))
			if (P_ < P_t) rho_ = (P_ / K_ )**(poly_n/(poly_n + 1.d0))/(1.d0 + drho)
			eos_polyFerrari2003 = 1.d0/((poly_n + 1.d0)/(poly_n) * (P_/rho_))
			endif
        else
            pause 'err: static profile eos table'
        endif
    endif
endfunction eos_polyFerrari2003

function eos_quark_matter_Alford2005(f, x_)
use parameters !requires module 'parameters' (variables: )
use constants !requires module 'constants' (variables: )
implicit none
real(8) :: eos_quark_matter_Alford2005, x_
character (len=*) :: f
real(8) :: C1_, A_, B_, C_, mu_2, P_, E_
    C1_ = 3.d0/4.d0/pi**(2.d0)
    A_ = C1_ * QS_a4
    B_ = - C1_ * QS_a2
    C_ = -B_eff
    
	if (f == 'p(rho)') then
	        E_ = (h_bar*c)**3 * x_/((1.d13)*e_mks)**4 * c**2  !Energy density E_ !Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
			mu_2 = -B_/6.d0/A_ + dsqrt(B_**2 + 12.d0*A_*(C_+E_))/6.d0/A_
	        eos_quark_matter_Alford2005 = (A_*mu_2**2 + B_*mu_2 + C_) * ((1.d13)*e_mks)**4 /(h_bar*c)**3
        elseif (f == 'rho(p)') then
            P_ = (h_bar*c)**3 * x_/((1.d13)*e_mks)**4
			mu_2 = -B_/2.d0/A_ + dsqrt(B_**2 - 4.d0*A_*(C_-P_))/2.d0/A_
	        eos_quark_matter_Alford2005 = (3.d0*A_*mu_2**2 + B_*mu_2 - C_) * ((1.d13)*e_mks)**4/ (h_bar*c)**3/ c**2
        elseif (f == 'drhodp(p)') then
            P_ = (h_bar*c)**3 * x_/((1.d13)*e_mks)**4
			mu_2 = -B_/2.d0/A_ + dsqrt(B_**2 - 4.d0*A_*(C_-P_))/2.d0/A_
			eos_quark_matter_Alford2005 = (6.d0*A_*mu_2**2 + B_*mu_2)/(2.d0*A_*mu_2**2 + B_*mu_2)/ c**2 
        else
            pause 'err: static profile eos table'
    endif
endfunction eos_quark_matter_Alford2005

function eos_incompressible_star(f, x_) !incompressible model with homogeneous density
use parameters !requires module 'parameters' (variables: )
use constants !requires module 'constants' (variables: )
implicit none
real(8) :: eos_incompressible_star, x_
character (len=*) :: f
real(8) :: F1, rho_cm, P_cm
    F1 = (1.d0 - 2.d0*ic_Comp)**(0.5d0)
        if (f == 'p(rho)') then ! only applicable for finding Central Pressure
            rho_cm = x_ * Grav_Const/c**2
            P_cm = (F1 - 1.d0) / (1.d0 - 3.d0*F1) * rho_cm
            eos_incompressible_star = P_cm / (Grav_Const/c**4)
        elseif (f == 'rho(p)') then
            rho_cm = rho_0 * Grav_Const/c**2
            eos_incompressible_star = rho_cm / (Grav_Const/c**2)
        elseif (f == 'drhodp(p)') then
            eos_incompressible_star = 0.d0
        else
            pause 'err: static profile eos table'
        endif
endfunction eos_incompressible_star

function eos_2layer_inc_star(f, x_) !2layer incompressible model with a density gap
use parameters !requires module 'parameters' (variables: )
use constants !requires module 'constants' (variables: )
implicit none
real(8) :: eos_2layer_inc_star, x_
character (len=*) :: f
real(8) :: F1, rho_cm, P_cm, P_t_cm
    P_t_cm  = P_t *(Grav_Const/c**4)
    if (f == 'p(rho)') then 
        rho_cm = x_ * Grav_Const/c**2
        F1 = (1.d0 - 2.d0*ic_Comp_Co)**(0.5d0) * (P_t_cm + rho_cm)/(3.d0*P_t_cm + rho_cm)
        P_cm = (F1 - 1.d0) / (1.d0 - 3.d0*F1) * rho_cm
        eos_2layer_inc_star = P_cm / (Grav_Const/c**4)
    elseif (f == 'rho(p)') then
        if ( x_ < P_t) then
            rho_cm = rho_0 * (1.d0 - drho) * Grav_Const/c**2
        elseif (x_ >= P_t) then
            rho_cm = rho_0 * Grav_Const/c**2
        endif
        eos_2layer_inc_star = rho_cm / (Grav_Const/c**2)
    else
        pause 'err: static profile eos table'
    endif
endfunction eos_2layer_inc_star

endmodule EOS