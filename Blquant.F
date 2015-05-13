      PROGRAM GETREK
C
C     Take "k" and chord location, output Re_k
C
C
#include "precis.h"
      CHARACTER(80) :: GFILE,QFILE,OFILE
      LOGICAL :: EXIST
      INTEGER :: IUNIT,OUNIT,J,K,L,QUNIT,IND
      INTEGER :: NGRID,MGRID,NQ,NQC,N,IBLANK,SD
      INTEGER :: NI,NJ,NK,NIDUMP,NJDUMP,NKDUMP,IG
      INTEGER :: JJ,LL,KK,KS,JS,JIND,JCUR,KIN,KCUR
      REAL :: SDIST,DIST,DISTK,KHEIGHT,REY,REK,REFMACH,
     &         UKLOW,UKHIGH,REKHIGH,REKLOW,DISTLOW,
     &         REKAVG,UKLOW2,UKHIGH2,DISTLEN,DISTINT,
     &         UKLOWAVG,UKHIGHAVG,DS
      REAL,  DIMENSION(:,:,:), ALLOCATABLE :: X,Y,Z
      REAL,  DIMENSION(:,:,:,:), ALLOCATABLE :: Q
      INTEGER,  DIMENSION(:), ALLOCATABLE :: JD,KD,LD
   
C
C    Get the input filenames
C
      WRITE(*,1)
 1    FORMAT('Enter PLOT3D grid filename:  ',$)
      READ(*,2,END=100) GFILE
 2    FORMAT(A)
C
C
      WRITE(*,22)
22    FORMAT('Enter PLOT3D qfile filename:  ',$)
      READ(*,2,END=100) QFILE
C 
C 
      IUNIT = 1
      QUNIT = 2
      OPEN(1,FILE=GFILE,STATUS='OLD',FORM='UNFORMATTED')
      READ(1,ERR=30) NI,NJ,NK
      REWIND(1)
      NGRID = 1
      MGRID = 1
      GOTO 5
 30   REWIND(1)
      READ(1) NGRID
      WRITE(*,4) 'Total number of grids is ',NGRID
 4    FORMAT(A,I5)
      WRITE(*,2,ADVANCE='NO') 'Which grid to read? '
      READ(*,*,END=100) MGRID
 5    IF(NGRID.EQ.1) THEN
         READ(1,ERR=6) NI,NJ,NK
      ELSE
         READ(1,ERR=6) (NI,NJ,NK,N=1,MGRID)
      ENDIF      
      GOTO 7
 6    CONTINUE
      REWIND(1)
      IF(NGRID.EQ.1) THEN
         READ(1) NI,NJ,NK
      ELSE
         READ(1)
         READ(1) (NI,NJ,NK,N=1,MGRID)
      ENDIF
      NQ = 5
 7    CONTINUE      
      WRITE(*,101)
      WRITE(*,8) NI,NJ,NK
 8    FORMAT('Reading Grid file, dimensioned ',3I7)
C
C
C      WRITE(*,20,ADVANCE='NO') 'Enter output file name: '
C 20   FORMAT(A)
C      READ(*,*,END=100) OFILE
C      OUNIT = 801
C
C     Skipping Unwanted Grids
C
      DO N = 1,MGRID-1
         READ(1)
      ENDDO
C
C     Allocate X,Y,Z arrays and call RDG which returns loaded arrays  
C
      ALLOCATE( X(NI,NJ,NK), Y(NI,NJ,NK), Z(NI,NJ,NK))
C
      CALL RDG(NI,NJ,NK,NQ,NQC,IUNIT,X,Y,Z)
C
      CLOSE(IUNIT)
C
C
      OPEN(2,FILE=QFILE,STATUS='OLD',FORM='UNFORMATTED')
C
      READ(QUNIT,ERR=33) NI,NJ,NK
      REWIND(QUNIT)
      GOTO 50
 33   REWIND(QUNIT)
      READ(2) NGRID
 50   IF(NGRID.EQ.1) THEN
         READ(QUNIT,ERR=60) NI,NJ,NK,NQ,NQC
      ELSE
         READ(QUNIT,ERR=60) (NI,NJ,NK,N=1,MGRID),
     &     (NIDUMP,NJDUMP,NKDUMP,N=1,NGRID-MGRID),NQ,NQC
      ENDIF
      GOTO 70
 60   CONTINUE
      REWIND(QUNIT)
      IF(NGRID.EQ.1) THEN
         READ(QUNIT) NI,NJ,NK
      ELSE
         READ(QUNIT)
         READ(QUNIT) (NI,NJ,NK,N=1,MGRID)
      ENDIF
 70   CONTINUE
C
C
      DO N = 1,MGRID-1
         READ(QUNIT)
         READ(QUNIT)
      ENDDO
C
      ALLOCATE( Q(NI,NJ,NK,NQ) )
      CALL RDQ(NI,NJ,NK,NQ,NQC,QUNIT,QFILE,MGRID,Q,
     &          REFMACH,REY)
C
C
C
      CLOSE(QUNIT)    
      WRITE(*,2,ADVANCE='NO') 'Enter output file name: '
      READ(*,*,END=100) OFILE
      OUNIT = 801
      OPEN(OUNIT,FILE=OFILE,STATUS='REPLACE',FORM='FORMATTED')
C
C
C
      WRITE(*,101)
      WRITE(*,81,ADVANCE='NO')
   81 FORMAT('Input k/c (Roughness height/chord): ')
      READ(*,*,END=100) KS
C
C      WRITE(*,82,ADVANCE='NO')
C   82 FORMAT('Input s/c location:  ')
C      READ(*,*,END=100) SDIST
C
C
C
C
C
C
      DO IND = 1,NI

      JIND    = IND
      DISTK   = 0
      KHEIGHT = KS*(0.000001)
      KCUR    = 2
      DO WHILE (DISTK.LT.KHEIGHT)
         DISTK = 0
         DO KIN = 2,KCUR
            DISTK = DISTK+SQRT((X(JIND,KIN,2) - X(JIND,(KIN-1),2))**2 +
     &                         (Z(JIND,KIN,2) - Z(JIND,(KIN-1),2))**2 )
         ENDDO
         KCUR = KCUR + 1 
         IF (KCUR.GT.NJ) THEN
            DISTK = DISTK+1           
         ENDIF
      ENDDO 
      DISTLOW = DISTK - SQRT((X(JIND,KIN,2) - X(JIND,(KIN-1),2))**2 +
     &                         (Z(JIND,KIN,2) - Z(JIND,(KIN-1),2))**2 )
C      WRITE(*,84) JCUR,KCUR
C   
C     Obtain U at the two J indicies nearest the choosen coordinate.
C
      UKLOW  = Q((JIND),(KIN-1),2,2)  
      UKHIGH = Q((JIND),KIN,2,2)
C      
C     Linearly interpolate along the j-direction      
C
C 
      REKLOW = REY*KHEIGHT*UKLOW/REFMACH
      REKHIGH = REY*KHEIGHT*UKHIGH/REFMACH
C         
      REKAVG = ((REKHIGH-REKLOW)/(DISTK-DISTLOW))*(KHEIGHT-DISTK) 
     &           + REKHIGH
      WRITE(*,103) JIND,REKAVG
      WRITE(OUNIT,103) JIND,REKAVG
 103  FORMAT(I5,F15.7)
      ENDDO

C     Gradients are much more substancial in the K direction
C     Output the direct distance the velocity is extracted from.
C      
      WRITE(*,101)
      WRITE(*,87) DISTK*1000000
  87  FORMAT('Distance (z/c [E-6] ) at (high) = ', F7.3)      
      WRITE(*,88) DISTLOW*1000000
  88  FORMAT('Distance (z/c [E-6] ) at (low) =', F7.3)      
C
C 
      WRITE(*,101)
 101  FORMAT('')
      WRITE(*,86) REKHIGH 
  86  FORMAT('RE_k (high) = ',F8.4)  
      WRITE(*,85) REKLOW 
  85  FORMAT('RE_k (low) = ',F8.4)      
C
      WRITE(*,101)
      WRITE(*,89) REKAVG
  89  FORMAT('Weighted Average RE_k  = ',F8.4)   

C 
  84  FORMAT(2I5)
 100  CONTINUE
      STOP
      END
C
C
C
      SUBROUTINE RDG(NI,NJ,NK,NQ,NQC,IUNIT,X,Y,Z)
#include "precis.h"
      INTEGER, INTENT (IN) :: NI,NJ,NK,NQ,NQC,IUNIT
      __REAL,  DIMENSION(NI,NJ,NK), INTENT (OUT) :: X,Y,Z
C
      INTEGER :: JJ,LL,KK
C
C      
      READ(IUNIT) (((X(JJ,KK,LL),JJ=1,NI),KK=1,NJ),LL=1,NK), 
     &        (((Y(JJ,KK,LL),JJ=1,NI),KK=1,NJ),LL=1,NK),
     &        (((Z(JJ,KK,LL),JJ=1,NI),KK=1,NJ),LL=1,NK)
C 
      RETURN
      END
C
C
      SUBROUTINE RDQ(NI,NJ,NK,NQ,NQC,QUNIT,QFILE,MGRID,Q,
     &                REFMACH,REY)
#include "precis.h"
      INTEGER, INTENT (IN) :: NI,NJ,NK,NQ,NQC,QUNIT,MGRID
C
      CHARACTER(80), INTENT (IN) :: QFILE
      __REAL, INTENT (OUT) :: REY,REFMACH
      __REAL,  DIMENSION(NI,NJ,NK,NQ), INTENT (OUT) :: Q
      __REAL,  DIMENSION(:), ALLOCATABLE :: RGAS
      __REAL  :: ALFA,TIME,GAMINF,BETA,TINF,
     &           HTINF,HT1,HT2,FSMACH,TVREF,DTVREF
      INTEGER :: IGAM,N
C
      INTEGER :: IS,IE,II,I,JS,JE,JI,J,KS,KE,KI,K
      INTEGER :: NQS,NQE
      CHARACTER :: ANS
C
 2    FORMAT(A)
C
      ALLOCATE( RGAS(MAX(2,NQC)) )
      READ(QUNIT) REFMACH,ALFA,REY,TIME,GAMINF,BETA,TINF,
     &            IGAM,HTINF,HT1,HT2,RGAS(1:MAX(2,NQC)),
     &            FSMACH,TVREF,DTVREF
      WRITE(*,40) 'refmach=',REFMACH,'alpha=',ALFA,'Re=',REY
 40   FORMAT(A,F5.2,4X,A,F8.2,3X,A,G12.5,3X,A,G12.5,3X)
      DEALLOCATE( RGAS )
C
      READ(QUNIT) Q
      CLOSE(QUNIT)
C     
C     
      RETURN
      END
