      PROGRAM BLQUANT
C
C     Take "k" and chord location, output Re_k
C
C
#include "precis.h"
      CHARACTER(80) :: GFILE,QFILE,OFILE
      LOGICAL :: EXIST
      INTEGER :: IUNIT,OUNIT,J,K,L,QUNIT,IND,MAXVEL
      INTEGER :: NGRID,MGRID,NQ,NQC,N,IBLANK,SD,MAXIND
      INTEGER :: NI,NJ,NK,NIDUMP,NJDUMP,NKDUMP,IG
      INTEGER :: JJ,LL,KK,KS,JS,JIND,JCUR,KIN,KCUR
      REAL :: SDIST,DIST,DISTK,KHEIGHT,REY,REK,REFMACH,
     &         UKLOW,UKHIGH,REKHIGH,REKLOW,DISTLOW,H,
     &         REKAVG,UKLOW2,UKHIGH2,DISTLEN,DISTINT,MV,
     &         UKLOWAVG,UKHIGHAVG,DS,DX,DY,DV,DU,VORT,
     &         INVOR,MAXVOR,RHOE,MOMTHK,DISTHK,U,UE,RHO,
     &         MTHK,DTHK,HT,IDTK,IHK,IDTKM,IHKM,UM,RHOM,
     &         RETH  
      REAL,  DIMENSION(:,:,:), ALLOCATABLE :: X,Y,Z
      REAL,  DIMENSION(:,:,:,:), ALLOCATABLE :: Q
      INTEGER,  DIMENSION(:), ALLOCATABLE :: JD,KD,LD
      LOGICAL :: HITMAX 
C
C     Going to hard code input file names for now
C
C      WRITE(*,1)
C 1    FORMAT('Enter PLOT3D grid filename:  ',$)
C      READ(*,2,END=100) GFILE
 2    FORMAT(A)
C
C
C      WRITE(*,22)
C22    FORMAT('Enter PLOT3D qfile filename:  ',$)
C      READ(*,2,END=100) QFILE
C 
C
      GFILE = 'grid.in'
      QFILE = 'q.save'  
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
C      WRITE(*,101)
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
C     Allocate Q var array and read the q file
      ALLOCATE( Q(NI,NJ,NK,NQ) )
      CALL RDQ(NI,NJ,NK,NQ,NQC,QUNIT,QFILE,MGRID,Q,
     &          REFMACH,REY)
C
C
C
      CLOSE(QUNIT)    
C     
C      WRITE(*,2,ADVANCE='NO') 'Enter output file name: '
C      READ(*,*,END=100) OFILE
      OFILE = 'BL_info.txt' 
      OUNIT = 801
      OPEN(OUNIT,FILE=OFILE,STATUS='REPLACE',FORM='FORMATTED')
C
C     
      DO JIND = 2,NI-1
          INVOR  = 0.0
          MAXVOR = 0.0
          MAXIND = 1
          DO KIN = 2,NJ-1 
             DX = X(JIND+1,KIN,2) - X(JIND-1,KIN,2)
             DY = Z(JIND,KIN+1,2) - Z(JIND,KIN-1,2) 
             DU = Q(JIND,KIN+1,2,2) - Q(JIND,KIN-1,2,2)  
             DV = Q(JIND+1,KIN,2,3) - Q(JIND-1,KIN,2,3)  
             VORT = DV/DX - DU/DY 
             INVOR = INVOR + VORT*(Z(JIND,KIN,2)-Z(JIND,KIN-1,2))    
             IF (ABS(INVOR) .GT. MAXVOR) THEN
                MAXVOR = ABS(INVOR)
             END IF
          END DO
C         
          HITMAX = .False.
C     Compute again until you hit .995 max
          INVOR = 0.0 
          DO KIN = 2,NJ -1
             DX = X(JIND+1,KIN,2) - X(JIND-1,KIN,2)
             DY = Z(JIND,KIN+1,2) - Z(JIND,KIN-1,2) 
             DU = Q(JIND,KIN+1,2,2) - Q(JIND,KIN-1,2,2)  
             DV = Q(JIND+1,KIN,2,3) - Q(JIND-1,KIN,2,3)  
             VORT = DV/DX - DU/DY 
             INVOR = INVOR + VORT*(Z(JIND,KIN,2)-Z(JIND,KIN-1,2))   
             IF ( ABS(INVOR) .gt. MAXVOR*.9995 ) THEN 
                IF (HITMAX .eqv. .False.) THEN
                   HITMAX = .True.
                   MAXIND = KIN 
C                   WRITE(*,*) ABS(INVOR),MAXVOR
                END IF 
             END IF
          END DO
C         More hueristic approach where we just look for max velocity  
          MV = 0.0
          HITMAX = .False.
          DO KIN = 2,NJ
             IF ( Q(JIND,KIN,2,2) .gt. MV 
     &           .and. HITMAX .eqv. .False.) THEN
                MV      = Q(JIND,KIN,2,2) 
                MAXVEL  = KIN
             ELSE 
                HITMAX = .True.   
             END IF
          END DO 
C         Some limiter (BAD!!)          
C          MAXIND = MIN(MAXIND,80)
          RHOE   = Q(JIND,MAXIND,2,1)
          UE     = Q(JIND,MAXIND,2,2)    
C         Compute momentum and displacement thickness
          DISTHK = 0.0
          MOMTHK = 0.0  
          DO KIN = 2,MAXIND
             U    = Q(JIND,KIN,2,2) 
             RHO  = Q(JIND,KIN,2,1)  
             UM   = Q(JIND,KIN-1,2,2) 
             RHOM = Q(JIND,KIN-1,2,1)  
             DY   = SQRT((Z(JIND,KIN,2) - Z(JIND,KIN-1,2))**2 +
     &              (X(JIND,KIN,2) - X(JIND,KIN-1,2))**2)  
C             IDTK,IHK,IDTKM,IHKM
C    
             IDTK   = (1- (RHO*U/(RHOE*UE)))
             IDTKM  = (1- (RHOM*UM/(RHOE*UE)))
             IHK    = (RHO*U/(RHOE*UE))*(1-(U/UE))
             IHKM   = (RHOM*UM/(RHOE*UE))*(1-(UM/UE))
             DISTHK = DISTHK + 0.5*(IDTK+IDTKM)*DY
             MOMTHK = MOMTHK + 0.5*(IHK+IHKM)*DY
             RETH   = MOMTHK*REY 
             H      = DISTHK/MOMTHK  
          END DO
C         Try a slightly different way to comute momenthick
          KIN = 2
          U = Q(JIND,KIN,2,2)  
          DO KIN = 2,MAXVEL
             UE   = MV
             U    = Q(JIND,KIN,2,2) 
             RHO  = Q(JIND,KIN,2,1)  
             UM   = Q(JIND,KIN-1,2,2) 
             RHOM = Q(JIND,KIN-1,2,1)  
             DY   = SQRT((Z(JIND,KIN,2) - Z(JIND,KIN-1,2))**2 +
     &              (X(JIND,KIN,2) - X(JIND,KIN-1,2))**2)  
C             IDTK,IHK,IDTKM,IHKM
C    
             IDTK   = (1- (RHO*U/(RHOE*UE)))
             IDTKM  = (1- (RHOM*UM/(RHOE*UE)))
             IHK    = (RHO*U/(RHOE*UE))*(1-(U/UE))
             IHKM   = (RHOM*UM/(RHOE*UE))*(1-(UM/UE))
             DTHK   = DISTHK + 0.5*(IDTK+IDTKM)*DY
             MTHK   = MOMTHK + 0.5*(IHK+IHKM)*DY
             HT     = DISTHK/MOMTHK   
          END DO  
          WRITE(OUNIT,*) X(JIND,1,2),MAXIND,JIND,DTHK,MTHK,HT
          WRITE(*,*) X(JIND,1,2),MAXIND,JIND,DISTHK,MOMTHK,H,RETH 
      END DO 


C      WRITE(*,103) JIND
C      WRITE(OUNIT,103) JIND
C 103  FORMAT(I5)
       

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
