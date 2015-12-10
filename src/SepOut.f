	MODULE SepOut

	IMPLICIT NONE

	PRIVATE
	PUBLIC :: S_printfiles


	INTEGER :: i,j

	CONTAINS

!----------------------

	SUBROUTINE S_printfiles(fileP,Model_name,
	1	f,fmin,fmax,largesep,smallsep,f_r010,r010,mr_ls,mr_ls_T,mr_ss,mr_ss_T,
	2	m_r02,sl_r02,m_r010,sl_r010,sl_Lsr010,
	3	Chi_st,Chi_fref,Chi_Ls,Chi_ss0,Chi_r02,Chi_sl_r02,Chi_r010,Chi_sl_r010,Chi_sl_Lsr010,Chi_st8,
	4	age,rad,lum,teff,zsx,logg,rcz,
	5	Chi_stLsss0,factor4Like)

	LOGICAL, INTENT(IN), DIMENSION(:) :: fileP
	CHARACTER(LEN=*), INTENT(IN) :: Model_name
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: f,largesep,smallsep
	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: f_r010,r010,mr_ls,mr_ss
	DOUBLE PRECISION, INTENT(IN) :: fmin,fmax,mr_ls_T,mr_ss_T,m_r02,sl_r02,m_r010,sl_r010,sl_Lsr010
	DOUBLE PRECISION, INTENT(IN) ::	Chi_st,Chi_fref,Chi_Ls,Chi_ss0,Chi_r02,Chi_sl_r02,Chi_r010,
	1	Chi_sl_r010,Chi_sl_Lsr010,Chi_st8
	DOUBLE PRECISION, INTENT(IN) :: age,rad,lum,teff,zsx,logg,rcz
	 DOUBLE PRECISION, INTENT(IN) :: Chi_stLsss0, factor4Like

	IF(fileP(1)) THEN
	  OPEN(unit=101,access='append',file='Seismparams.dat')
	  WRITE(101,'(A,7F12.6)')TRIM(Model_name),mr_ls_T,mr_ss(lbound(f,1)),
	1	m_r02,sl_r02,m_r010,sl_r010,sl_Lsr010
	  CLOSE(unit=101)
c	  WRITE(*,*)'	Seismparams.dat 		... done'
	ENDIF

	IF(fileP(2)) THEN
	  OPEN(unit=102,form='formatted',status='unknown',access='append',file='Chi.dat')
          WRITE(102,'(A,10F8.3)')model_name,Chi_st,Chi_fref,Chi_Ls,Chi_ss0,Chi_r02,Chi_sl_r02,Chi_r010,
	1	Chi_sl_r010,Chi_sl_Lsr010,Chi_st8
	  CLOSE(unit=102)
c	  WRITE(*,*)'	Chi.dat				... done' 
	ENDIF

	IF(fileP(3)) THEN
	  OPEN(unit=103,form='formatted',status='unknown',access='append',file='StelChars.dat')
          WRITE(103,'(A,1p7E12.3)')model_name,age,rad,lum,teff,zsx,logg,rcz
          CLOSE(unit=103)
c         WRITE(*,*)'   StelChars.dat                         ... done'
        ENDIF

	IF(fileP(4)) THEN
	  OPEN(unit=104,form='formatted',status='unknown',access='append',file=TRIM(Model_name)//'.ech')
          DO j=lbound(f,1),ubound(f,1)                                        !l
            DO i=2,ubound(f,2)                                  !n  Separations start at n=2
              IF((f(j,i)>=fmin).AND.(f(j,i)<=fmax).AND.(f(j,i)/=0.)) WRITE(104,*)MODULO(f(j,i),mr_ls_T),f(j,i)
c	      IF((f(j,i)>=fmin).AND.(f(j,i)<=fmax).AND.(f(j,i)/=0.)) WRITE(104,*)MODULO(f(j,i),0.0878),f(j,i)
            ENDDO
	  ENDDO
	  CLOSE(unit=104)
c	  WRITE(*,*)'	',TRIM(Model_name),'.ech	... done'
        ENDIF

	IF(fileP(5)) THEN
	  OPEN(unit=105,form='formatted',status='unknown',access='append',file=TRIM(Model_name)//'.r010')
	  DO i=1,SIZE(r010)
	    IF(f_r010(i)/=0.) WRITE(105,'(2F16.6)')f_r010(i),r010(i)
	  ENDDO
	  CLOSE(unit=105)
	ENDIF

	IF(fileP(6)) THEN
          OPEN(unit=106,form='formatted',status='unknown',access='append',file=TRIM(Model_name)//'.Ls012')
	  DO j=0,2
            DO i=1,SIZE(largesep(j,:))
              IF(largesep(j,i)/=0.) WRITE(106,'(2F16.6)')f(j,i),largesep(j,i)
            ENDDO
	  ENDDO
          CLOSE(unit=106)
        ENDIF

	IF(fileP(7)) THEN
          OPEN(unit=107,form='formatted',status='unknown',access='append',file=TRIM(Model_name)//'.Chi2')
	  WRITE(107,'(2F20.10)')Chi_stLsss0,factor4Like
	  CLOSE(unit=107)
	ENDIF

	END SUBROUTINE S_printfiles
!---------------------

	END MODULE SepOut
