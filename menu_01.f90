module menu_01

!Description:
!    This module contains the setting menu menu_01
!    The subroutine 'me01_body' contains the main flow

implicit none

public :: me01_body
private :: me01_modify_settings, me01_modify_parameters

contains

subroutine me01_body !main flow of the menu
implicit none
logical,intent(in) :: start
integer,intent(in) :: ans
    start = .false.
    do while (start /= .true.)
	    write(*,*) "=========================================================="
		write(*,*) "Options:"
		write(*,*) "1. Begin"
		write(*,*) "2. Modify settings"
		write(*,*) "3. Modify parameters"
		write(*,*) "=========================================================="
		read(*,*) ans
		write(*,*) "=========================================================="
		if (ans == 1) then
			start = .true.
			write(*,*) "=========================================================="
			write(*,*) "		Program Begins"
			write(*,*) "=========================================================="
		elseif (ans == 2) then
			call me01_modify_settings !found below
		elseif (ans == 3) then
			call me01_modify_parameters  !found below
		else
			write(*,*) 'invalid input'
		endif
	enddo
	return
endsubroutine me01_body

! --------------------------------------------------------------------
!Private subroutines/ functions
! --------------------------------------------------------------------

subroutine me01_modify_settings !modifies the settings/ options, e.g. options for integrating the perturbation equations
use parameters
implicit none
integer,intent(in) :: ans
    ans = 999
    do while (ans /= 0)
	    write(*,*) "Inputs:"
        write(*,*) "	(1) Static Profile Options"
        write(*,*) "	(2) Perturbation Grid Options"
        write(*,*) "	(3) Perturbation Problem Solver Options"
        write(*,*) "	(4) Perturbation Equation Iterations Options"
        write(*,*) "	(5) Perturbation Equation Combination Options"
        write(*,*) "	(6) EOS Options"
        write(*,*) "	(7) Shear Modulus Input Options"
        write(*,*) "	(8) Shear Modulus Formula Options"
3       write(*,*) "	(0) Back"
        write(*,*) "=========================================================="
        read(*,*) ans
        write(*,*) "=========================================================="
		    if (ans == 1) then
			elseif (ans == 2) then
			elseif (ans == 0) then
				return
		    else
			    write(*,*) 'invalid input'
			endif
    enddo
endsubroutine me01_modify_settings

subroutine me01_modify_parameters !modifies the parameters, e.g. central density rho_0
use parameters
implicit none
integer,intent(in) :: ans
    ans = 999
    do while (ans /= 0)
            write(*,*) "Inputs:"
            write(*,*) "	(1) Angular momentum degree"
            write(*,*) "	(2) Central Density"
            write(*,*) "	(3) Crust Core Transition Pressure"
            write(*,*) "	(4) Ocean Crust Transition Pressure"
            write(*,*) "	(5) Cutoff Pressure"
            write(*,*) "	(6) Polytrope: Density Jump Fraction"
            write(*,*) "	(7) Polytrope: Crust Core Transition Density"
            write(*,*) "	(8) Polytrope: Shear Speed throughout the crust"
            write(*,*) "	(9) Quark Star: a4"
            write(*,*) "	(10) Quark Star: a2"
            write(*,*) "	(11) Quark Star: B_eff"
            write(*,*) "	(12) Quark Star: Gap Parameter"
            write(*,*) "	(13) Polytrope: Polytropic index"
            write(*,*) "	(0) Back"
            write(*,*) "=========================================================="
            read(*,*) ans
            write(*,*) "=========================================================="
		    	if (ans == 1) then
		    	elseif (ans == 2) then
		    	elseif (ans == 0) then
			    	return
				else
			        write(*,*) 'invalid input'
			    endif
	enddo
endsubroutine me01_modify_parameters

endmodule menu_01