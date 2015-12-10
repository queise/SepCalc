	PROGRAM SepCalc
! *************************************************************************************
!
!
!	This program calculates several seismic parameters (large and small separations,
!	ratios...) from a given set of frequencies of a stellar model, and then compares
!	the results with the observed values.
!
!
!  *** 	Modules:	SepCalc.f 	: Main program
!			SepFunct.f	: Functions
!			SepOut.f	: Prints output files
!
!
!  ***	Structure of main program SepCalc.f:
!
!	1) Reads information from input.dat file, and the frequencies of the model and
!	of the star.
!	2) Calculates the large and small separations, means, second differences, scaledss02,
!	scaledss13, dr0213, r02, r13, r01, r10, r010, and slopes.
!	3) Compares the results with the observed values, calculating chi squareds.
!
!
! ***	Input files:	
!
!	- input.dat :
!	- Model_name.freqs : 	as the format of the output of ADIPLS program:
!				n,l,order,sigma**2,E,v.frequency,Ri.frequency,case
!	- Star.freqs:		observed frequencies, format:
!				l,n,frequency,error
!	- Star.mean-Err-Ls.ss0, .mean-Err-r02.r010 and .sl-Err-r02.r010.Lsr010:	with value error value error...
!
!
! ***	Output files:
!	 	fileP(1)	Seismparams.dat
!        	fileP(2)	Chi.dat
!        	fileP(3)	StelChars.dat
!        	fileP(4)    	echelle diagram
!        	fileP(5)	r010(f)
!		fileP(6)	Ls(f)
!		fileP(7)	Chi2 and factor4Likelihood 1/2piErr
!
!
! ***	References:
!	- Separations and general formulae: 
!	Tassoul ApJSS 43 (1980), Lopes & Turck-Chieze A&A 290 (1994), Christensen-Dalsgaard, J. (1988)
!	- 5-point separations, ratios and slopes:
!	Roxburgh & Vorontsov A&A 411 (2003), Cunha & Metcalfe ApJ 666 (2007), arXiv: 0705.3009
!	Brandao, Cunha & Christensen-Dalsgaard (2014), arXiv: 1311.7600
!
! ************************************************************************************
!	Author: Jordi Casanellas   
!	Date: November 2014
!	SepCalc version: 2.2.6		   
! **********************************
!
	USE SepFunct
	USE SepOut

	IMPLICIT NONE

	INTEGER :: i,j,k,l,n, iloop, n_linesS, n_lines, n_maxS, n_max, kmax
	INTEGER, ALLOCATABLE, DIMENSION(:) :: n_last, locfmin, locfmax, l_of_r010
        INTEGER :: lmin, lmax, locfminr010, locfmaxr010, locfref
! Seismic parameters:
	DOUBLE PRECISION :: fmin, fmax, fminsl, fmaxsl
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: fS, Err_fS, f, largesep, smallsep, secdiff
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m_ls, m_ss, mr_ls, mrsl_ls, mrr_ls, mr_ss, scaledss02, scaledss13, 
	1	dr0213,	ratio02, ratio13, d01, ratio01, d10, ratio10, r010, f_r010, Lsr010
	DOUBLE PRECISION :: m_ls_T, m_ss_T, mr_ls_T, mrsl_ls_T, mr_ss_T, sl_dr0213, m_r02, sl_r02, m_r01, sl_r01,
	1	m_r10, sl_r10, m_r010, sl_r010, sl_Lsr010
	DOUBLE PRECISION :: LsSTAR,Err_LsSTAR,ss0STAR,Err_ss0STAR,r02STAR,Err_r02STAR,r010STAR,Err_r010STAR,
	1	sl_r02STAR,Err_sl_r02STAR,sl_r010STAR,Err_sl_r010STAR,sl_Lsr010STAR,Err_sl_Lsr010STAR
! Stellar characteristics of Model and Star:
	DOUBLE PRECISION :: age,rad,lum,teff,zsx,logg,rcz
	DOUBLE PRECISION :: loggSTAR,Err_loggSTAR,teffSTAR,Err_teffSTAR,zsxSTAR,Err_zsxSTAR
! Chi for Model-Star comparision:
	DOUBLE PRECISION :: Chi_st,Chi_fref,Chi_Ls,Chi_ss0,Chi_r02,Chi_sl_r02,Chi_r010,Chi_sl_r010,
	1	Chi_sl_Lsr010,Chi_st8,Chi_sl_dr0213,Chi_stLs,Chi_stLsss0,Chi_stslr101,
	2	Chi_logg,Chi_teff,Chi_zsx, factor4Like
	TYPE lfS        ! format of a line in the Star.freqs (observed frequencies)
          INTEGER :: ls , ns
          DOUBLE PRECISION :: Fr, Err
        END TYPE
	TYPE(lfS), ALLOCATABLE, DIMENSION(:) :: matfS
	TYPE lfM	! format of a line in the Model.freqs file, output of ADIPLS code
	  INTEGER :: num, ls , ns 
	  DOUBLE PRECISION :: sig2, En, F1, F2
	  INTEGER :: num2
	END TYPE
	TYPE(lfM), ALLOCATABLE, DIMENSION(:) :: matfM

	LOGICAL :: fileP(7)

	CHARACTER(LEN=100) :: Star_name, Model_name

c       WRITE(*,*)'-------------------------------------------'
c       WRITE(*,*)'             SEPCALC   v2.2
c       WRITE(*,*)'-------------------------------------------'
c       WRITE(*,*)' Jordi Casanellas, April 2014' ; WRITE(*,*)''

! *********************************************************************
! *****************    1) INPUTS        *******************************
! *********************************************************************

! ... Reading parameters of the observed freqs for better further comparision in chicalc:
	OPEN(unit=100,file='input.dat')
        READ(100,*) ; READ(100,'(A)')Star_name
c	WRITE(*,*)'Star_name: ',Star_name
	READ(100,*) ; READ(100,*)fminsl,fmaxsl 
c	WRITE(*,*)'fminsl=',fminsl,' fmaxsl=',fmaxsl
	READ(100,*) ; READ(100,*)loggSTAR,Err_loggSTAR
        READ(100,*) ; READ(100,*)teffSTAR,Err_teffSTAR
        READ(100,*) ; READ(100,*)zsxSTAR,Err_zsxSTAR
	READ(100,*) ; READ(100,*) ; READ(100,*) ; READ(100,*)
	READ(100,*) ; READ(100,'(A)')Model_name
	WRITE(*,*)TRIM(Model_name),'.freqs'
	DO i=1,7; READ(100,*) ; READ(100,*) fileP(i) ; ENDDO
        CLOSE(unit=100)

! ... Reading data:
! ..from observed Star:
        OPEN(unit=100,form='formatted',status='unknown',file=TRIM(Star_name)//'.freqs')
        DO 991 iloop = 1,10000
        READ(100,*,END=992)
991     CONTINUE
992     n_linesS  = iloop - 1
        REWIND(100)
        ALLOCATE (matfS(n_linesS))
        READ(100,*)matfS
        CLOSE(unit=100)

! ..from Model (format as adipls output):
        OPEN(unit=101,form='formatted',status='unknown',file=TRIM(Model_name)//'.freqs')
        DO 993 iloop = 1,10000
        READ(101,*,END=994)
993     CONTINUE
994     n_lines  = iloop - 1
        IF(n_lines==0) STOP
        REWIND(101)
        ALLOCATE (matfM(n_lines))
        READ(101,*)matfM
        CLOSE(unit=101)

! ... Extracting relevant information from data and storing data:
! ..from observed Star:
	lmin = MINVAL(matfS%ls) ; lmax = MAXVAL(matfS%ls)
	IF(lmin/=0) THEN
           WRITE(*,*)' ATTENTION '
           WRITE(*,*)' lmin different than zero produces error in functions, lmin=',lmin
           STOP
        ENDIF
	fmin = MINVAL(matfS%Fr)/1000. ; fmax = MAXVAL(matfS%Fr)/1000.	!  from muHz to mHz
	n_maxS = MAXVAL(matfS%ns)
	ALLOCATE (fS(lmin:lmax,n_maxS),Err_fS(lmin:lmax,n_maxS)) ; fS=0.
        DO l=lmin,lmax
          DO i=1,n_linesS
            IF((matfS(i)%ls)==l) THEN
	      fS(l,matfS(i)%ns) = matfS(i)%Fr
	      Err_fS(l,matfS(i)%ns) = matfS(i)%Err
	    ENDIF     
          ENDDO
        ENDDO
! ..from Model:
	n_max = MAXVAL(matfM%ns)
	ALLOCATE (f(lmin:lmax,n_max)) ; f=0.
	DO l=lmin,lmax
	  DO i=1,n_lines
	    IF((matfM(i)%ls)==l) f(l,matfM(i)%ns) = matfM(i)%F2
	  ENDDO
	ENDDO
c	DO l=lmin,lmax
c          DO n=1,n_max
c	    WRITE(*,*)'f(',l,',',n,')=',f(l,n)
c	  ENDDO
c	ENDDO

	ALLOCATE(locfmin(lmin:lmax),locfmax(lmin:lmax))
	locfmin = MINLOC(ABS(f(:,:)-fminsl),2,(f(:,:)-fminsl)>=0.) ! location(n) of 1st freq > fminsl, for each l
	locfmax = MINLOC(ABS(f(:,:)-fmaxsl),2,(f(:,:)-fmaxsl)<=0.) ! location(n) of last freq < fmaxsl, for each l

! *********************************************************************
! ****************    2) SEPARATIONS   ********************************
! *********************************************************************
! Refs: Lopes & Turck-Chieze A&A 290 (1994)
!	Cunha & Metcalfe ApJ 666 (2007), arXiv: 0705.3009

! ... Calculating separations ...
	ALLOCATE (largesep(lmin:lmax,n_max),smallsep(lmin:(lmax-2),n_max),secdiff(lmin:lmax,n_max))
	largesep = 0. ; smallsep = 0. ; secdiff = 0.
c       WRITE(*,*)' Calculating separations...'
	largesep = F_largesep(f) !	; WRITE(*,*)'	large separation 		... done'
	smallsep = F_smallsep(f) !	; WRITE(*,*)' 	small separation		... done'
	secdiff  = F_secdiff(f)	 !	; WRITE(*,*)' 	2nd diffs  			... done'
	
! ... Mean separations 
	ALLOCATE (m_ls(lmin:lmax),m_ss(lmin:(lmax-2)),mr_ls(lmin:lmax),
	1	mrsl_ls(lmin:lmax),mrr_ls(lmin:lmax),mr_ss(lmin:(lmax-2)))
	m_ls  = SUM(largesep,2,largesep/=0) / COUNT(largesep/=0,2)
	m_ls_T = SUM(m_ls) / SIZE(m_ls)
	m_ss = SUM(smallsep,2,smallsep/=0) / COUNT(smallsep/=0)
	m_ss_T = SUM(m_ss) / SIZE(m_ss)
c	WRITE(*,*)'	means                 		... done'

! ... Mean separations within ranges
	mr_ls = SUM(largesep,2,((largesep/=0.).AND.(f>fmin).AND.(f<fmax))) / 
	1	COUNT((largesep/=0.).AND.(f>fmin).AND.(f<fmax),2)
	mr_ls_T = SUM(mr_ls) / SIZE(mr_ls)
	mrsl_ls = SUM(largesep,2,((largesep/=0.).AND.(f>fminsl).AND.(f<fmaxsl))) /
	1	COUNT((largesep/=0.).AND.(f>fminsl).AND.(f<fmaxsl),2)
	mrsl_ls_T = SUM(mrsl_ls) / SIZE(mrsl_ls)
	mr_ss = SUM(smallsep,2,((smallsep/=0.).AND.(f>fmin).AND.(f<fmax))) /
	1	COUNT((smallsep/=0.).AND.(f>fmin).AND.(f<fmax),2)
        mr_ss_T = SUM(mr_ss) / SIZE(mr_ss)
c	WRITE(*,*)'	means within range		... done'

! ... Other seismic parameters:
c	WRITE(*,*)' Calculating other seismic parameters...'
	ALLOCATE (scaledss02(n_max),scaledss13(n_max),dr0213(n_max),ratio02(n_max),ratio13(n_max),
	1	d01(n_max),ratio01(n_max),d10(n_max),ratio10(n_max),r010(2*n_max),f_r010(2*n_max),
	2	l_of_r010(2*n_max),Lsr010(2*n_max))
        scaledss02 = 0. ; scaledss13 = 0. ; dr0213 = 0. ; sl_dr0213 = 0. ; ratio02 = 0. ; m_r02 = 0. ; sl_r02 = 0.
	ratio13 = 0. ; d01 = 0. ; ratio01 = 0. ; m_r01 = 0. ; sl_r01 = 0. ; d10 = 0. ; ratio10 = 0. ; m_r10 = 0. ; sl_r10 = 0.
	r010 = 0. ; l_of_r010 = 0 ; f_r010 = 0. ; m_r010 = 0. ; sl_r010 = 0. ; Lsr010 = 0. ; sl_Lsr010 = 0.
! ... scaled small separations:
	scaledss02 = smallsep(0,:)/6.
	IF(lmax>2) scaledss13 = smallsep(1,:)/6.
c        WRITE(*,*)'	scaled small separations        ... done'

! ... Delta dr0213 (Cunha & Metcalfe 07)
	IF(lmax<3) GOTO 995	! needs modes of l=3 to be calculated
	DO i=2,n_max
	  IF((largesep(1,i-1)/=0.).AND.(largesep(0,i)/=0.)) dr0213(i)=
	1 (scaledss02(i)/largesep(1,i-1)) - (scaledss13(i)/largesep(0,i)) 
	ENDDO	
	sl_dr0213 = F_slope(f(0,:),dr0213,locfmin(0),locfmax(0)) ! slope between freq range
c	WRITE(*,*)'     dr0213 ...done'
995     CONTINUE

! ... Ratios (Roxburgh & Vorontsov 93)
! ..ratio02
	DO i=1,n_max
	  IF(largesep(1,i)/=0.) ratio02(i) = smallsep(0,i) / largesep(1,i)
	ENDDO
	m_r02 = SUM(ratio02,((f(0,:)>fminsl).AND.(f(0,:)<fmaxsl))) /
	1	COUNT((ratio02/=0.).AND.(f(0,:)>fminsl).AND.(f(0,:)<fmaxsl))
        sl_r02 = F_slope(f(0,:),ratio02,locfmin(0),locfmax(0))

! ..ratio13
	DO i=1,(n_max-1)
          IF ((smallsep(1,i)/=0.).AND.(largesep(0,i+1)/=0.)) ratio13(i)=smallsep(1,i)/largesep(0,i+1)
	ENDDO
! ..ratio01 and d01
	DO i=2,(n_max-1)
          IF ((f(0,(i-1))/=0.).AND.(f(1,(i-1))/=0.).AND.(f(0,i)/=0.).AND.(f(1,i)/=0.).AND.(f(0,(i+1))/=0.)) THEN
	    d01(i)=(f(0,(i-1))-4*f(1,(i-1))+6*f(0,i)-4*f(1,i)+f(0,(i+1)))/8.
	    IF(largesep(1,i)/=0.) ratio01(i) = d01(i) / largesep(1,i)
	  ENDIF
	ENDDO
	m_r01 = SUM(ratio01,((f(0,:)>fminsl).AND.(f(0,:)<fmaxsl))) /
	1	COUNT((ratio01/=0.).AND.(f(0,:)>fminsl).AND.(f(0,:)<fmaxsl))
        sl_r01 = F_slope(f(0,:),ratio01,locfmin(0),locfmax(0))	
! ..ratio10 and d10
	DO i=2,(n_max-1)
          IF ((f(1,(i-1))/=0.).AND.(f(0,i)/=0.).AND.(f(1,i)/=0.).AND.(f(0,(i+1))/=0.).AND.(f(1,(i+1))/=0.)) THEN
	    d10(i)=-(f(1,(i-1))-4*f(0,i)+6*f(1,i)-4*f(0,(i+1))+f(1,(i+1)))/8.
	    IF(largesep(0,(i+1))/=0.) ratio10(i)=d10(i)/largesep(0,(i+1))
	  ENDIF
        ENDDO
        m_r10 = SUM(ratio10,((f(1,:)>fminsl).AND.(f(1,:)<fmaxsl))) /
	1	COUNT((ratio10/=0.).AND.(f(1,:)>fminsl).AND.(f(1,:)<fmaxsl))
	sl_r10 = F_slope(f(1,:),ratio10,locfmin(1),locfmax(1))
! ..ratio010
        k=1
	DO i=2,(n_max-1)
          IF (ratio01(i)>0.) THEN	! avoid negative ratios
            r010(k)=ratio01(i)
            f_r010(k)=f(0,i)
            l_of_r010(k)=1
            k=k+1
          ENDIF
          IF (ratio10(i)>0.) THEN
            r010(k)=ratio10(i)
            f_r010(k)=f(1,i)
            l_of_r010(k)=0
            k=k+1
          ENDIF
        ENDDO
        kmax=k-1
	m_r010 = SUM(r010,((f_r010(:)>fminsl).AND.(f_r010(:)<fmaxsl))) /
	1	COUNT((r010/=0.).AND.(f_r010(:)>fminsl).AND.(f_r010(:)<fmaxsl))
	locfminr010 = MINLOC(ABS(f_r010(:)-fminsl),1,(f_r010(:)-fminsl)>=0.)
	locfmaxr010 = MINLOC(ABS(f_r010(:)-fmaxsl),1,(f_r010(:)-fmaxsl)<=0.)
        sl_r010 = F_slope(f_r010,r010,locfminr010,locfmaxr010)
!       slope of product Ls * r010 (Brandao, Cunha, CD 2014)
c	Lsr010 = mr_ls(l_of_r010(:)) * r010
	Lsr010 = mrsl_ls(l_of_r010(:)) * r010	! mean of Ls within the slope range
	sl_Lsr010 = F_slope(f_r010,Lsr010,locfminr010,locfmaxr010)
c	mrr_ls = SUM(largesep,2,((largesep/=0.).AND.(f>f_r010(kmax-6)).AND.(f<f_r010(kmax)))) /
c	1	COUNT((largesep/=0.).AND.(f>f_r010(kmax-6)).AND.(f<f_r010(kmax)),2)
c	Lsr010 = mrr_ls(l_of_r010(:)) * r010
c	sl_Lsr010 = F_slope(f_r010,Lsr010,kmax-6,kmax)
!	Absolute value of slope
	sl_Lsr010 = ABS(sl_Lsr010)

! *********************************************************************
! ***************    3) CHI2 COMPARISION STAR-MODEL    ****************
! *********************************************************************

! ... Reading STAR seismic parameters (previously calculated)
	OPEN(unit=91,file=TRIM(STAR_name)//'.mean-Err-Ls.ss0')
        READ(91,*)LsSTAR,Err_LsSTAR,ss0STAR,Err_ss0STAR
        CLOSE(unit=91)
        OPEN(unit=91,file=TRIM(STAR_name)//'.mean-Err-r02.r010')
        READ(91,*)r02STAR,Err_r02STAR,r010STAR,Err_r010STAR
        CLOSE(unit=91)
        OPEN(unit=91,file=TRIM(STAR_name)//'.sl-Err-r02.r010.Lsr010')
        READ(91,*)sl_r02STAR,Err_sl_r02STAR,sl_r010STAR,Err_sl_r010STAR,sl_Lsr010STAR,Err_sl_Lsr010STAR
        CLOSE(unit=91)
	sl_Lsr010STAR = ABS(sl_Lsr010STAR)
c	WRITE(*,*)'       STAR seismic parameters read 	... done'

! ... Reading model stellar characteristics:
	OPEN(unit=100,file=TRIM(model_name)//'.fin')
        READ(100,*)age,rad,lum,teff,zsx,logg,rcz
        CLOSE(unit=100)

! ... First freq of star at l=0 (fmin) as a reference to compare with model:
        locfref = MINLOC(ABS(f(0,:)-fmin),1)

! ... Calculating Chi2
	Chi_fref = ABS((fmin - f(0,locfref))/0.020)
	Chi_Ls   = ABS((LsSTAR - mr_ls_T)/Err_LsSTAR)	! ..the mean Large separation (total l=0,1,2...)
        Chi_ss0  = ABS((ss0STAR - mr_ss(0))/Err_ss0STAR)	! ..the mean small separation (only l=0)
        Chi_r02  = ABS((r02STAR - m_r02)/Err_r02STAR)
        Chi_r010 = ABS((r010STAR - m_r010)/Err_r010STAR)
c        Chi_sl_dr0213 = ABS( (sl_dr0213STAR - sl_dr0213)/Err_sl_dr0213STAR)
        Chi_sl_dr0213 =0.
        Chi_sl_r010   = ABS((sl_r010STAR - sl_r010)/Err_sl_r010STAR)
        Chi_sl_r02    = ABS((sl_r02STAR - sl_r02)/Err_sl_r02STAR)
        Chi_sl_Lsr010 = ABS((sl_Lsr010STAR - sl_Lsr010)/Err_sl_Lsr010STAR)	! ..slope of product Ls and r010
        Chi_teff = ABS(( teffSTAR - teff ) / Err_teffSTAR )
        Chi_zsx  = ABS(( zsxSTAR - zsx ) / Err_zsxSTAR )
        Chi_logg = ABS(( loggSTAR - logg ) / Err_loggSTAR )
        Chi_st = sqrt(Chi_teff**2. + Chi_zsx**2. + Chi_logg**2.) / 3.
        Chi_stLs = sqrt(Chi_teff**2. + Chi_zsx**2. + Chi_logg**2. + Chi_Ls**2.) / 4.		! ..the stellar characteristics and large separation together:
        Chi_stLsss0 = sqrt(Chi_teff**2. + Chi_zsx**2. + Chi_logg**2. + Chi_Ls**2. + Chi_ss0**2.) / 5.
        Chi_stslr101 = sqrt(Chi_teff**2. + Chi_zsx**2. + Chi_logg**2. + Chi_sl_r010**2.) / 4.
        Chi_st8 = sqrt(Chi_teff**2. + Chi_zsx**2. + Chi_logg**2. + Chi_fref**2. + Chi_Ls**2. + Chi_ss0**2. +
	1	Chi_r02**2. + Chi_sl_r02**2. + Chi_r010**2. + Chi_sl_r010**2. + Chi_sl_Lsr010**2.) / 11.
	factor4Like = ((2.*3.14159265)**(-5./2.)) / (Err_teffSTAR*Err_zsxSTAR*Err_loggSTAR*Err_LsSTAR*Err_ss0STAR)

! *********************************************************************
! ***************   4)  OUTPUT FILES    *******************************
! *********************************************************************
c	WRITE(*,*)' Printing output files...'
	CALL S_printfiles(fileP,
	2	TRIM(Model_name),
	1	f,fmin,fmax,largesep,smallsep,f_r010,r010,mr_ls,mr_ls_T,mr_ss,mr_ss_T,
	2	m_r02,sl_r02,m_r010,sl_r010,sl_Lsr010,
	3	Chi_st,Chi_fref,Chi_Ls,Chi_ss0,Chi_r02,Chi_sl_r02,Chi_r010,Chi_sl_r010,Chi_sl_Lsr010,Chi_st8,
	4	age,rad,lum,teff,zsx,logg,rcz,
	5	Chi_stLsss0, factor4Like)

	STOP
	END PROGRAM SepCalc

