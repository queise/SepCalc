	MODULE sepfunct

	IMPLICIT NONE

	PRIVATE
        PUBLIC :: F_largesep, F_smallsep, F_secdiff, F_slope

	CONTAINS

!----------------------

	FUNCTION F_largesep(f)

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: f

	INTEGER :: i,j,j2	
	DOUBLE PRECISION, DIMENSION(0:(ubound(f,1)-1),ubound(f,2)) :: F_largesep

	F_largesep = 0.
	DO j=1,(ubound(f,1))
	  j2 = j - 1			! ASSUMES lmin = 0
	  DO i=2,ubound(f,2)
	    IF((f(j,i)==0.).OR.(f(j,(i-1))==0.)) GOTO 310
	    F_largesep(j2,i) = f(j,i) - f(j,(i-1))
c            WRITE(*,'(2I5,3F12.6)')j2,i,f(j,i),f(j,(i-1)),F_largesep(j2,i)
310	  ENDDO
	ENDDO

	END FUNCTION F_largesep
!----------------------

        FUNCTION F_smallsep(f)

        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: f

        INTEGER :: i,j,j2
        DOUBLE PRECISION, DIMENSION(0:(ubound(f,1)-3),ubound(f,2)) :: F_smallsep

        F_smallsep = 0.
	DO j=1,(ubound(f,1)-2)
	  j2 = j - 1
	  DO i=2,ubound(f,2)
            IF((f(j,i)==0.).OR.(f((j+2),(i-1))==0.)) GOTO 320
            F_smallsep(j2,i)=f(j,i)-f((j+2),(i-1))
c	    WRITE(*,'(2I5,3F12.6)')j2,i,f(j,i),f(j,(i-1)),F_smallsep(j2,i)
320       ENDDO
        ENDDO

	END FUNCTION F_smallsep
!----------------------

	FUNCTION F_secdiff(f)

        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: f

        INTEGER :: i,j,j2
        DOUBLE PRECISION, DIMENSION(0:(ubound(f,1)-1),ubound(f,2)) :: F_secdiff

        F_secdiff = 0.

	DO j=1,(ubound(f,1))
	  j2 = j - 1
	  DO i=2,(ubound(f,2)-1)
	    IF((f(j,i)==0.).OR.(f(j,(i+1))==0.).OR.(f(j,(i-1))==0.)) GOTO 330
	    F_secdiff(j2,i)=f(j,(i+1))-2*f(j,i)+f(j,(i-1))
c	    WRITE(*,'(2I5,3F12.6)')j2,i,f(j,i),f(j,(i-1)),F_secdiff(j2,i)
330       ENDDO
        ENDDO

        END FUNCTION F_secdiff
!----------------------

	FUNCTION F_slope(x,y,p1,p2)

        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: x,y
	INTEGER, INTENT(IN) :: p1,p2
	INTEGER :: i,np2,np1
        DOUBLE PRECISION :: F_slope

	i=p2
	DO WHILE (y(i)==0.)
	  i=i-1
	ENDDO
	np2=i
	i=p1
	DO WHILE (y(i)==0.)
	  i=i+1
	ENDDO
	np1=i

	F_slope = (y(np2)-y(np1)) / (x(np2)-x(np1))	
  
c	WRITE(*,*)'F_slope ',y(np2),'-',y(np1),' / ',x(np2),'-',x(np1)

        END FUNCTION F_slope
!----------------------

	END MODULE sepfunct
