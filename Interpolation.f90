MODULE Interpolation
!*******************************************************************
!	Interpolation:
!	#	API for interpolation
!	Int_api		- function name, returns Int_api = y(xi) in a set of discrete data points (x,y)
!	
!	mode		- choose option for interpolation (e.g. linear interpolation, polynomial interpolation, cubic spline interpolation, etc)
!	i1, i2		- data range, must enclose the output range
!	x, y		- array of discrete data points (x,y)
!	xi			- desired x value
!
!*******************************************************************
CONTAINS

!************************************************************************************
!	Interpolation control unit
!************************************************************************************
function Int_api(mode, i1, i2, x, y, xi)
implicit none
!	Instructions: 
!				- input x, y as discrete points in j-space
!				- input i1, i2 as range of discrete points x, y in j-space
!				- input desired value xi
!				- the code returns Int_api as y(xi) value
!				- input parameter 'mode' according to the following descriptions
!
!     Input parameter 'mode' is a string. Put in the following codes inside '' to select 
!   the mode of interpolation. Modify this API to integrate more interpolation choices
!   to the code.
!
!	modes:
!		linear_Ln				- linear ln interpolation
!		linear					- linear interpolation
!		4pt_poly				- 4 pts polynomial interpolation
!		linear_Ln_slope			- linear ln interpolation finding slope
!
!------------------------------------------------------------------------------------
real(8) :: x(:), y(:), xi
real(8) :: Int_api, Polint_out
integer :: i1, i2, j, w
character(len=50) :: mode
	do j = i1, i2 - 1
		if ( (xi - x(j)) * (xi - x(j+1)) <= 0.d0) then
			if (mode == 'linear_Ln') then
				Int_api = linear_Ln_int(x(j), x(j+1), y(j), y(j+1), xi)
			elseif (mode == 'linear') then
				Int_api = linear_int(x(j), x(j+1), y(j), y(j+1), xi)
			elseif (mode == '4pt_poly') then
				w = j
				if (w == i1 .or. w == i1 +1) w = i1 + 2
				if (w == i2) w = i2 -1
				call polint(x(w-2:w+1),y(w-2:w+1),4,xi,Polint_out)
				Int_api = Polint_out
			elseif (mode == 'linear_Ln_slope') then
				Int_api = linear_Ln_Slope(x(j), x(j+1), y(j), y(j+1), xi)
			else
				pause 'err: interpolation mode'
			endif
			exit
		endif
		if (j == i2 - 1) then
			pause 'err: interpolation'
		endif
	enddo
endfunction Int_api

!************************************************************************************
!	Linear Interpolation
!************************************************************************************
FUNCTION linear_int(x1,x2,y1,y2,x)
REAL(8)	:: x1, x2, y1, y2, x
REAL(8)	:: linear_int
	linear_int = (y2 - y1)/(x2 - x1) * (x - x1) + y1
ENDFUNCTION linear_int

!************************************************************************************
!	Linear Ln Interpolation
!************************************************************************************
FUNCTION linear_Ln_int(x1,x2,y1,y2,x)
REAL(8)	:: x1, x2, y1, y2, x
REAL(8)	:: linear_Ln_int
	if (x1*x2*y1*y2*x == 0.d0) then
		write(*,*) "err: linear Ln interpolation; zero arguements"
		pause
	endif
	linear_Ln_int = dexp((dlog(y2) - dlog(y1))/(dlog(x2) - dlog(x1)) * (dlog(x) - dlog(x1)) + dlog(y1))
ENDFUNCTION linear_Ln_int

!************************************************************************************
!	Linear Ln Interpolation
!************************************************************************************
FUNCTION linear_Ln_Slope(x1,x2,y1,y2,x)
REAL(8), INTENT(IN) :: x1, x2, y1, y2, x
REAL(8) :: linear_Ln_Slope, y
	if (x1*x2*y1*y2*x == 0.d0) then
		write(*,*) "err: linear Ln interpolation; zero arguements"
		pause
	endif
	y = dexp((dlog(y2) - dlog(y1))/(dlog(x2) - dlog(x1)) * (dlog(x) - dlog(x1)) + dlog(y1))
	linear_Ln_Slope = (dlog(y2) - dlog(y1))/(dlog(x2) - dlog(x1)) * y/x
ENDFUNCTION linear_Ln_Slope

!************************************************************************************
!	Polynomial Interpolation (Numerical Recipes)
!************************************************************************************
SUBROUTINE polint(xa,ya,n,x,y)
INTEGER n,NMAX
REAL(8) dy,x,y,xa(n),ya(n)
PARAMETER (NMAX=10) !	Largest anticipated value of n.
!Given arrays xa and ya, each of length n, and given a value x, this routine returns a
!value y, and an error estimate dy. If P(x) is the polynomial of degree N ? 1 such that
!P(xai ) = yai, i = 1,..., n, then the returned value y = P(x).
INTEGER i,m,ns
REAL(8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

	ns=1
	dif=abs(x-xa(1))
	do i=1,n 
		dift=abs(x-xa(i))
	if (dift.lt.dif) then
		ns=i
		dif=dift
	endif
	c(i)=ya(i)
	d(i)=ya(i)
	enddo
	y=ya(ns) 
	ns=ns-1
	do m=1,n-1 
	do i=1,n-m 
		ho=xa(i)-x
		hp=xa(i+m)-x
		w=c(i+1)-d(i)
		den=ho-hp
	if(den.eq.0.)pause
		den=w/den
		d(i)=hp*den
		c(i)=ho*den
	enddo
	if (2*ns.lt.n-m)then
		dy=c(ns+1)
	else
		dy=d(ns)
		ns=ns-1
	endif
	y=y+dy
	enddo
	return
ENDSUBROUTINE polint

!!	Reference
SUBROUTINE POLINT1(XA,YA,N,X,Y,DY)
!*****************************************************
!*     Polynomial Interpolation or Extrapolation     *
!*            of a Discreet Function                 *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*   XA:    Table of abcissas  (N)                   *
!*   YA:    Table of ordinates (N)                   *
!*    N:    Number of points                         *
!*    X:    Interpolating abscissa value             *
!* OUTPUTS:                                          *
!*    Y:    Returned estimation of function for X    *
!*   DY:    Estimated error for Y                    *
!*****************************************************
PARAMETER(NMAX=25)
REAL*8 XA(N),YA(N), X,Y,DY
REAL*8 C(NMAX),D(NMAX)
REAL*8 DEN,DIF,DIFT,HO,HP,W
NS=1
DIF=DABS(X-XA(1))
DO I=1,N
  DIFT=DABS(X-XA(I))
  IF(DIFT.LT.DIF) THEN
    NS=I                 !index of closest table entry
	DIF=DIFT
  ENDIF
  C(I)=YA(I)             !Initialize the C's and D's
  D(I)=YA(I)
END DO
Y=YA(NS)                 !Initial approximation of Y
NS=NS-1
DO M=1,N-1
  DO I=1,N-M
    HO=XA(I)-X
	HP=XA(I+M)-X
    W=C(I+1)-D(I) 
    DEN=HO-HP
	IF(DEN.EQ.0.) Pause 'Error: two identical abcissas)'
	DEN=W/DEN
	D(I)=HP*DEN          !Update the C's and D's
	C(I)=HO*DEN
  END DO
  IF(2*NS.LT.N-M) THEN   !After each column in the tableau XA
    DY=C(NS+1)           !is completed, we decide which correction,
  ELSE                   !C or D, we want to add to our accumulating 
    DY=D(NS)             !value of Y, i.e. which path to take through 
	NS=NS-1              !the tableau, forking up or down. We do this
  ENDIF                  !in such a way as to take the most "straight 
  Y=Y+DY	             !line" route through the tableau to its apex,
END DO                   !updating NS accordingly to keep track of 
                         !where we are. This route keeps the partial
RETURN                   !approximations centered (insofar as possible)
ENDSUBROUTINE POLINT1    !on the target X.The last DY added is thus the
                         !error indication.
ENDMODULE Interpolation