program cs_polar_deform
use main
implicit none
character (len=30) :: ans

01  continue
	write(*,*) " Main Menu:"
    write(*,*) "=========================================================="
! main body begin
    call main_body
! main body end
02  write(*,*) "Redo? (y/n)"
    read(*,*) ans
    if (ans == "y") then
        write(100,*) "END. Redo."
        write(100,*) "=========================================================="
        goto 01
    elseif (ans == "n") then
        write(*,*) "Bye"
    else 
        pause'err:wrong output!'
        write(*,*) "=========================================================="
        goto 02
    endif
endprogram cs_polar_deform