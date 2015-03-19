      INCLUDE 'penfield.f'
C
C  This program generates electron and positron trajectories in uniform
C  electromagnetic fields in vacuum.
C
C  To generate the executable binary file, compile and link TRAJECT.F
C  with PENFIELD.F.
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,LBL
C  ****  EM field.
      COMMON/UFIELD/EFX,EFY,EFZ,BFX,BFY,BFZ
C
      MAT=1
      IBODY=1
C
      OPEN(7,FILE='track')
      WRITE(6,*) 'Components of the electric field (V/cm)?'
      READ(5,*) EFX,EFY,EFZ
      WRITE(7,'('' # E-field (V/cm)  ='',1P,3E15.7)') EFX,EFY,EFZ
      WRITE(6,*) 'Magnitude of the magnetic field (gauss)?'
      READ(5,*) BFX,BFY,BFZ
      WRITE(7,'('' # M-field (gauss) ='',1P,3E15.7)') BFX,BFY,BFZ
C
      WRITE(6,*) 'Electron (-1) or positron (+1)?'
      READ(5,*) KPAR
      IF(KPAR.EQ.-1) THEN
        WRITE(7,*) '# Electron'
        KPAR=1
      ELSE
        WRITE(7,*) '# Positron'
        KPAR=3
      ENDIF
      X=0.0D0
      Y=0.0D0
      Z=0.0D0
      WRITE(7,'('' # Initial position ='',1P,3E15.7)') X,Y,Z
      WRITE(6,*) 'Angles of initial direction (deg)?'
      READ(5,*) TH,PH
      U=DSIN(TH*PI/180.0D0)*DCOS(PH*PI/180.0D0)
      V=DSIN(TH*PI/180.0D0)*DSIN(PH*PI/180.0D0)
      W=DCOS(TH*PI/180.0D0)
      WRITE(7,'('' # Initial direction (theta,phi) ='',1P,3E15.7)')
     1   TH,PH
      WRITE(7,'('' # Initial direction ='',1P,3E15.7)') U,V,W
C
      WRITE(6,*) 'Initial energy (eV)?'
      READ(5,*) E
      WRITE(7,*) '# Energy (eV) =',E
C
      WRITE(6,*) 'Trajectory length (cm)?'
      READ(5,*) TL
      WRITE(7,*) '# Trajectory length (cm) =',TL
C
      WRITE(6,*) 'Upper limit of direction change?'
      READ(5,*) ULDV
      WRITE(7,*) '# Upper limit of direction change =',ULDV
C
      WRITE(6,*) 'Upper limit of energy change?'
      READ(5,*) ULDE
      WRITE(7,*) '# Upper limit of energy change =',ULDE
C
      WRITE(6,*) 'Upper limit of field change?'
      READ(5,*) ULEM
      WRITE(7,*) '# Upper limit of field change =',ULEM
C
      WRITE(7,*) '   '
      WRITE(7,*) '   '
      TREM=TL
      write(7,2000) x,y,z,e
 2000 format(1x,1p,4e15.7)
C
  100 CONTINUE
      CALL TPEMF0(ULDV,ULDE,ULEM,DSMAX)
C  ----  Moving the particle a path length DSMAX.
      DS=DSMAX
      IF(DS.GT.TREM) DS=TREM
      CALL TPEMF1(DS,DSEF,NCROSS)
      TREM=TREM-DSEF
      write(7,2000) x,y,z,e
      IF(TREM.GT.1.0D-10) GO TO 100
C
      CLOSE(7)
      STOP
      END

C  *********************************************************************
C  The following is a simplified STEP subroutine that moves particles
C  within and infinite medium.
      SUBROUTINE STEP(DSR,DSEF,NCROSS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      DSEF=DSR
      X=X+U*DSR
      Y=Y+V*DSR
      Z=Z+W*DSR
      NCROSS=0
      RETURN
      END
