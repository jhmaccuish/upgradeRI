MODULE SEARCH_A

USE PARAMS
USE GLOBALS
USE FUNCTIONS
USE NUMERICAL_LIBRARIES   
USE SEARCH_INT

IMPLICIT NONE

CONTAINS

!=========================================================
SUBROUTINE ASEARCH_Rx
!=========================================================
IMPLICIT NONE
!=========================================================
INTEGER:: amax_G
REAL(prec):: vmax
REAL(prec):: vtemp
!=========================================================

vmax	= -1.0E+15           
amax_G	= -1

askipR	= REAL(ac_max - ac_min)/REAL(intA)
askip	= CEILING(askipR)


DO acc0 = 1,intA+1  
	
	acc = ac_min + askipR * (acc0-1)  

    cons = ( x1 - grida2(acc) ) / onetauc

	IF (cons<0.0) GOTO 102

    vtemp = utilR(cons) + one_sv * utilB(acc)

    IF (vtemp>vmax) THEN      
		vmax = vtemp
        amax_G = acc
	ELSE
		GOTO 102
    ENDIF

ENDDO ! acc0


102 continue

IF (askip<2) GO TO 109

askip = askip/2

IF (amax_G>ac_min) THEN  

	acc = amax_G - askip
    cons = ( x1 - grida2(acc) ) / onetauc

    IF (cons<0.0) PRINT *,'WARNING: NEGATIVE CONS @ ASEARCH_Rx'

    vtemp = utilR(cons) + one_sv * utilB(acc)
    
    IF (vtemp>vmax) THEN
		vmax   = vtemp
		amax_G = acc
        GOTO 102
    ENDIF

ENDIF 


IF (amax_G < ac_max) THEN  

	acc  = amax_G + askip
    cons = ( x1 - grida2(acc) ) / onetauc

	IF (cons<0.0) GO TO 102

    vtemp = utilR(cons) + one_sv * utilB(acc)

    IF (vtemp>vmax) THEN
		vmax = vtemp
		amax_G = acc
    ENDIF

ENDIF 

GOTO 102

109 CONTINUE

vfunR(jc,ac,sc)   = vmax    ! jc=Nj
afunR_G(jc,ac,sc) = amax_G  ! jc=Nj

ENDSUBROUTINE ASEARCH_Rx



!=========================================================
SUBROUTINE ASEARCH_R
IMPLICIT NONE
INTEGER:: amax_G
REAL(prec):: vmax
REAL(prec):: vtemp

vmax	= -1.0E+15           
amax_G	= -1

askipR	= REAL(ac_max - ac_min)/REAL(intA)
askip	= CEILING(askipR)
                    		    
                 		    
DO acc0 = 1,intA+1  
	
	acc = ac_min + askipR * (acc0-1)  

    cons = ( x1 - grida2(acc) ) / onetauc

	IF (cons<0.0) GO TO 202

	CALL VINT_R  ! RETURNS vtemp99

    vtemp = utilR(cons) + bbbb_sv * vtemp99 +  one_sv * utilB(acc) 

    IF (vtemp>vmax) THEN      
		vmax = vtemp
        amax_G = acc
	ELSE
		GOTO 202
    ENDIF

ENDDO

202 continue

IF (askip<2) GO TO 209

askip = askip/2

IF (amax_G>ac_min) THEN  

	acc = amax_G - askip
    cons = ( x1 - grida2(acc) ) / onetauc

    IF (cons<0.0) PRINT *,'WARNING: NEGATIVE CONS @ ASEARCH_R'

	CALL VINT_R ! RETURNS vtemp99
    
    vtemp = utilR(cons) + bbbb_sv * vtemp99 +  one_sv * utilB(acc) 


	IF (vtemp>vmax) THEN
		vmax = vtemp
		amax_G = acc
        GOTO 202
    ENDIF

ENDIF 


IF (amax_G < ac_max) THEN  

	acc  = amax_G + askip
    cons = ( x1 - grida2(acc) ) / onetauc

	IF (cons<0.0) GO TO 202

    CALL VINT_R ! RETURNS vtemp99
   
    vtemp = utilR(cons) + bbbb_sv * vtemp99 +  one_sv * utilB(acc) 

    IF (vtemp>vmax) THEN
		vmax = vtemp
		amax_G = acc
    ENDIF

ENDIF 

GOTO 202


209 CONTINUE

vfunR(jc,ac,sc)   = vmax
afunR_G(jc,ac,sc) = amax_G


ENDSUBROUTINE ASEARCH_R


!=========================================================
SUBROUTINE ASEARCH_W 
IMPLICIT NONE
INTEGER:: amax_G
REAL(prec):: vmaxA
REAL(prec):: vtemp

vmaxA	= -1.0E+15           
amax_G	= -1

askipR	= REAL(ac_max - ac_min)/REAL(intA)
askip	= CEILING(askipR)


! EVALUATE 5 COARSE GRID POINTS

DO acc0 = 1,intA+1     
	
	acc  = ac_min + askipR * (acc0-1)  
    
	cons = (x0 - grida2(acc))/onetauc

    IF (cons<=0.0) GOTO 302  									
    
	IF (jc.eq.Njw) THEN 
		CALL VINT_W_R ! RETURNS vtemp99 GIVEN acc,ssbal
	ELSE 
		CALL VINT_W   ! RETURNS vtemp99 GIVEN acc,ssbal
	ENDIF
    
    vtemp = util(cons,lc,jc) + bbbb_sv * vtemp99 +  one_sv * utilB(acc) 

    IF (vtemp>vmaxA) THEN
		vmaxA	= vtemp
		amax_G	= acc
	ELSE
		GOTO 302
    ENDIF

ENDDO


!  SEARCH OVER REMAINING GRID POINTS

302 CONTINUE

IF (askip<2) GOTO 309

askip = askip/2

IF (amax_G>ac_min) THEN 

	acc  = amax_G - askip
	
	cons = (x0 - grida2(acc))/onetauc

    IF (cons<=0.0) PRINT *,'WARNING: NEGATIVE CONSUMPTION @ ASEARCH_W'
		! COULDN'T BE NEGATIVE (SINCE SAVING LESS THAN amax_G)!!									
    
	IF (jc.eq.Njw) THEN 
		CALL VINT_W_R ! RETURNS vtemp99 GIVEN acc,ssbal
	ELSE 
		CALL VINT_W   ! RETURNS vtemp99 GIVEN acc,ssbal
	ENDIF
    
    vtemp = util(cons,lc,jc) + bbbb_sv * vtemp99 +  one_sv * utilB(acc) 

    IF (vtemp>vmaxA) THEN
		vmaxA	= vtemp
		amax_G	= acc
		GOTO 302
    ENDIF
	
ENDIF


IF (amax_G < ac_max) THEN 

	acc  = amax_G + askip	
	    
	cons = (x0 - grida2(acc))/onetauc

    IF (cons<=0.0) GOTO 302  									
    
	IF (jc.eq.Njw) THEN 
		CALL VINT_W_R ! RETURNS vtemp99 GIVEN acc,ssbal
	ELSE 
		CALL VINT_W   ! RETURNS vtemp99 GIVEN acc,ssbal
	ENDIF
    
    vtemp = util(cons,lc,jc) + bbbb_sv * vtemp99 +  one_sv * utilB(acc) 

    IF (vtemp>vmaxA) THEN
		vmaxA	= vtemp
		amax_G	= acc
    ENDIF
	
ENDIF

GOTO 302

309 CONTINUE

! OUTPUT OF ASEARCH_W
vtempAA  = vmaxA
atemp_G  = amax_G


ENDSUBROUTINE ASEARCH_W



ENDMODULE SEARCH_A