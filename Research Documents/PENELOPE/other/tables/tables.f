      INCLUDE 'penelope.f'  ! File included to simplify compilation.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C           TTTTTT    AA    BBBBB   L       EEEEEE   SSSS              C
C             TT     A  A   B    B  L       E       S    S             C
C             TT    A    A  B    B  L       E       SS                 C
C             TT    AAAAAA  BBBBB   L       EEEE      SSSS             C
C             TT    A    A  B    B  L       E       S    S             C
C             TT    A    A  BBBBB   LLLLLL  EEEEEE   SSSS              C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  This program gives interpolated values of cross sections, mean free
C  paths and ranges of particles in the _first_ material of the input
C  material data file. This file must have been generated previously by
C  running the code MATERIAL. The program TABLES also generates detailed
C  tables of fundamental quantities used in the simulation, which are
C  written in the output file TABLES.DAT and in a number of separate
C  files with the extension '.TAB'.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2005)                                     C
C  Copyright (c) 2001-2005                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*1 YN
      CHARACTER*24 MFNAME
      CHARACTER*24 STR(3,6)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C
      STR(1,1)='not an interaction   (?)'
      STR(1,2)='elastic scattering      '
      STR(1,3)='inelastic scattering    '
      STR(1,4)='bremsstrahlung emission '
      STR(1,5)='inner-shell ionization  '
      STR(1,6)='not an interaction   (?)'
C
      STR(2,1)='Rayleigh scattering     '
      STR(2,2)='Compton scattering      '
      STR(2,3)='photoelectric absorption'
      STR(2,4)='pair production         '
      STR(2,5)='not an interaction   (?)'
      STR(2,6)='not an interaction   (?)'
C
      STR(3,1)='not an interaction   (?)'
      STR(3,2)='elastic scattering      '
      STR(3,3)='inelastic scattering    '
      STR(3,4)='bremsstrahlung emission '
      STR(3,5)='inner-shell ionization  '
      STR(3,6)='annihilation            '
C
C  ****  Parameters (to tabulate the complete energy range and to switch
C        soft interactions off).
C
      EMAX=1.0D9
      DO M=1,MAXMAT
        EABS(1,M)=50.0D0
        EABS(2,M)=50.0D0
        EABS(3,M)=50.0D0
        C1(M)=0.0D0
        C2(M)=0.0D0
        WCC(M)=0.0D0
        WCR(M)=-10.0D0
      ENDDO
C
C  ****  Material data file.
C
      WRITE(6,*) '  '
      WRITE(6,*) 'Enter the name of the material data file...'
      READ(5,'(A24)') MFNAME
      WRITE(6,'('' Material data file:  '', A24)') MFNAME
      WRITE(6,*) '  '
      WRITE(6,*) 'Processing material data. Please, wait...'
C
C  ****  We start by initializing PENELOPE...
C
      OPEN(11,FILE=MFNAME)
      OPEN(26,FILE='tables.dat')
      CALL PEINIT(EMAX,1,11,26,5)
      CLOSE(11)
      CLOSE(26)
      CALL GTAB
      CALL EPTAB
C
      WRITE(6,*) '  '
      WRITE(6,*) 'The output file ''tables.dat'' has been generated.'
      WRITE(6,*) '  '
      WRITE(6,*) 'Do you wish to calculate quantities for specific'
      WRITE(6,*) 'particles and energies?   (y/n)'
      READ(5,'(A1)') YN
      IF(YN.EQ.'N'.OR.YN.EQ.'n') STOP
      WRITE(6,*) '  '
C
C  ************  Now we can compute ranges and MFPs...
C
      M=1
      WRITE(6,*) 'Range (1) or mfp (2)?      (Ctrl-C stops the program)'
      READ(5,*) IRMFP
      IF(IRMFP.NE.1.AND.IRMFP.NE.2) STOP
      IF(IRMFP.EQ.2) GO TO 2
C
C  ****  Ranges.
C
    1 CONTINUE
      WRITE(6,*) '  '
      WRITE(6,*) 'Enter KPAR, E ...'
      READ(5,*) KPAR,E
      IF(E.LT.50.0D0.OR.E.GT.EMAX) THEN
        WRITE(6,*) '   Energy outside the allowed interval...'
        GO TO 1
      ENDIF
      IF(KPAR.EQ.3) THEN
        WRITE(6,1101) E
 1101   FORMAT(/1X,'Positron.    E =',1P,E11.4,' eV')
        RANGE=PRANGE(E,KPAR,M)
        WRITE(6,1111) RANGE,RANGE*RHO(1)
 1111   FORMAT(1X,'         range =',1P,E11.4,' cm ',
     1        /1X,'               =',E11.4,' g/cm**2')
      ELSE IF(KPAR.EQ.2) THEN
        WRITE(6,1102) E
 1102   FORMAT(/1X,'Photon.        E =',1P,E11.4,' eV')
        RANGE=PRANGE(E,KPAR,M)
        WRITE(6,1112) RANGE,RANGE*RHO(1)
 1112   FORMAT(1X,'  mean free path =',1P,E11.4,' cm ',
     1        /1X,'                 =',E11.4,' g/cm**2')
        RMU=1.0D0/RANGE
        WRITE(6,1212) RMU,RMU/RHO(1)
 1212   FORMAT(1X,'     att. coeff. =',1P,E11.4,' 1/cm ',
     1        /1X,'mass att. coeff. =',E11.4,' cm**2/g')
      ELSE
        KPAR=1
        WRITE(6,1103) E
 1103   FORMAT(/1X,'Electron.    E =',1P,E11.4,' eV')
        RANGE=PRANGE(E,KPAR,M)
        WRITE(6,1113) RANGE,RANGE*RHO(1)
 1113   FORMAT(1X,'         range =',1P,E11.4,' cm ',
     1        /1X,'               =',E11.4,' g/cm**2')
      ENDIF
      GO TO 1
C
C  ****  Mean free paths.
C
    2 CONTINUE
      WRITE(6,*) '  '
      WRITE(6,*) 'Enter KPAR, E, ICOL ...'
      READ(5,*) KPAR,E,ICOL
      IF(E.LT.50.0D0.OR.E.GT.EMAX) THEN
        WRITE(6,*) '   Energy outside the table interval...'
        GO TO 2
      ENDIF
      IF(ICOL.GT.6) ICOL=1
      IF(KPAR.EQ.3) THEN
        WRITE(6,1201) E,STR(KPAR,ICOL)
 1201   FORMAT(/1X,'Positron.    E =',1P,E11.4,' eV.',5X,A24)
      ELSE IF(KPAR.EQ.2) THEN
        WRITE(6,1202) E,STR(KPAR,ICOL)
 1202   FORMAT(/1X,'Photon.      E =',1P,E11.4,' eV.',5X,A24)
      ELSE
        KPAR=1
        WRITE(6,1203) E,STR(KPAR,ICOL)
 1203   FORMAT(/1X,'Electron.    E =',1P,E11.4,' eV.',5X,A24)
      ENDIF
C
      CMFP=PHMFP(E,KPAR,M,ICOL)
      WRITE(6,1204) CMFP,CMFP*RHO(1)
 1204 FORMAT(1X,'mean free path =',1P,E11.4,' cm',
     1      /1X,'               =',E11.4,' g/cm**2')
      WRITE(6,1205) 1.0D0/CMFP,1.0D0/(CMFP*RHO(1))
 1205 FORMAT(1X,'   inverse mfp =',1P,E11.4,' 1/cm',
     1      /1X,'               =',E11.4,' cm**2/g')
      XS=(1.0D0/CMFP)/VMOL(1)
      WRITE(6,1206) XS
 1206 FORMAT(1X,' cross section =',1P,E11.4,' cm**2')
      GO TO 2
C
      END
C  *********************************************************************
C                       SUBROUTINE GTAB
C  *********************************************************************
      SUBROUTINE GTAB
C
C  This subroutine prints tables of interaction properties of photons.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
      PARAMETER (NDIM=1500)
      DIMENSION ER(NDIM),XSR(NDIM),X1(NDIM),Y1(NDIM),X2(NDIM),Y2(NDIM)
C
      M=1
      OPEN(11,FILE='gxsceil.tab')
      WRITE(11,1101)
 1101 FORMAT(1X,'# PENELOPE >>>  Photelectric cross-section ceiling',
     1  ' (ph_ceil).')
      WRITE(11,1102)
 1102 FORMAT(1X,'#',/1X,'#   Energy',7X,'ph_ceil',
     1  /1X,'#    (keV)',7X,'(cm**2)',/1X,'#',26('-'))
      DO IE=1,NEGP-1
        WRITE(11,'(1X,1P,3E13.5)') ET(IE)*1.0D-3,SGPH(M,IE)/VMOL(M)
        WRITE(11,'(1X,1P,3E13.5)') ET(IE+1)*1.0D-3,SGPH(M,IE)/VMOL(M)
      ENDDO
      CLOSE(11)
C
      OPEN(12,FILE='gxs.tab')
      WRITE(12,1201)
 1201 FORMAT(1X,'# PENELOPE >>>  Photon cross sections (xs).')
      WRITE(12,1202)
 1202 FORMAT(1X,'#',/1X,'#   Energy',8X,'xs_Ra',8X,'xs_Co',8X,'xs_ph',
     1  8X,'xs_pp',7X,'xs_tot',/1X,'#    (keV)',7X,'(cm**2)',6X,
     2  '(cm**2)',6X,'(cm**2)',6X,'(cm**2)',6X,'(cm**2)',/1X,'#',
     3  78('-'))
C
      OPEN(13,FILE='gmacs.tab')
      WRITE(13,1301)
 1301 FORMAT(1X,'# PENELOPE >>>  Photon mass attenuation coefficie',
     1  'nts (mac).')
      WRITE(13,1302)
 1302 FORMAT(1X,'#',/1X,'#   Energy',7X,'mac_Ra',7X,'mac_Co',7X,
     1   'mac_ph',7X,'mac_pp',7X,'mac_tot',/1X,'#    (keV)',6X,
     2  '(cm**2/g)',4X,'(cm**2/g)',4X,'(cm**2/g)',4X,'(cm**2/g)',4X,
     3  '(cm**2/g)',/1X,'#',78('-'))
C
      OPEN(14,FILE='gmfp.tab')
      WRITE(14,1401)
 1401 FORMAT(1X,'# PENELOPE >>>  Photon mean free paths (mfp).')
      WRITE(14,1402)
 1402 FORMAT(1X,'#',/1X,'#   Energy',7X,'mfp_Ra',7X,'mfp_Co',7X,
     1   'mfp_ph',7X,'mfp_pp',7X,'mfp_tot',/1X,'#    (keV)',8X,
     2  '(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',/1X,'#',78('-'))
C
      IZZ=IZ(M,1)
      IST=IPHF(IZZ)
      LST=IPHL(IZZ)
      NPHD=0
      DO I=IST,LST
        NPHD=NPHD+1
        ER(NPHD)=EXP(EPH(I))
        XSR(NPHD)=STF(M,1)*EXP(XPH(I,1))
      ENDDO
      IF(NELEM(M).GT.1) THEN
        DO IEL=2,NELEM(M)
          N1=NPHD
          DO I=1,N1
            X1(I)=ER(I)
            Y1(I)=XSR(I)
          ENDDO
          IZZ=IZ(M,IEL)
          IST=IPHF(IZZ)
          LST=IPHL(IZZ)
          N2=0
          DO I=IST,LST
            N2=N2+1
            X2(N2)=EXP(EPH(I))
            Y2(N2)=STF(M,IEL)*EXP(XPH(I,1))
          ENDDO
          CALL MERGE2(X1,Y1,X2,Y2,ER,XSR,N1,N2,NPHD)
        ENDDO
      ENDIF
      DO I=1,NPHD
        E=ER(I)
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        XSRA=EXP(SGRA(M,KE)+(SGRA(M,KE+1)-SGRA(M,KE))*XEK)/VMOL(M)
        RFPRA=XSRA*VMOL(M)/RHO(M)
        FPRA=1.0D0/(RHO(M)*RFPRA)
        XSCO=EXP(SGCO(M,KE)+(SGCO(M,KE+1)-SGCO(M,KE))*XEK)/VMOL(M)
        RFPCO=XSCO*VMOL(M)/RHO(M)
        FPCO=1.0D0/(RHO(M)*RFPCO)
        RFPPH=XSR(I)*VMOL(M)/RHO(M)
        FPPH=1.0D0/(RHO(M)*RFPPH)
        IF(E.LT.1.023D6) THEN
          XSPP=1.0D-35
          RFPPP=1.0D-35
          FPPP=1.0D35
        ELSE
          XSPP=EXP(SGPP(M,KE)+(SGPP(M,KE+1)-SGPP(M,KE))*XEK)/VMOL(M)
          RFPPP=XSPP*VMOL(M)/RHO(M)
          FPPP=1.0D0/(RHO(M)*RFPPP)
        ENDIF
        XST=XSRA+XSCO+XSR(I)+XSPP
        RFPT=XST*VMOL(M)/RHO(M)
        FPT=1.0D0/(RHO(M)*RFPT)
        WRITE(12,'(1X,1P,6E13.5)') E*1.0D-3,XSRA,XSCO,XSR(I),XSPP,XST
        WRITE(13,'(1X,1P,6E13.5)') E*1.0D-3,RFPRA,RFPCO,RFPPH,RFPPP,
     1    RFPT
        WRITE(14,'(1X,1P,6E13.5)') E*1.0D-3,FPRA,FPCO,FPPH,FPPP,FPT
      ENDDO
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EPTAB
C  *********************************************************************
      SUBROUTINE EPTAB
C
C  This subroutine prints tables of interaction properties of electrons
C  and positrons.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Total and transport cross sections.
      COMMON/CEEL00/EIT(NEGP),XE0(NEGP),XE1(NEGP),XE2(NEGP),XP0(NEGP),
     1  XP1(NEGP),XP2(NEGP),T1E0(NEGP),T2E0(NEGP),T1P0(NEGP),T2P0(NEGP),
     2  EITL(NEGP),FL(NEGP),A(NEGP),B(NEGP),C(NEGP),D(NEGP)
C
      M=1
C
C  ****  Electrons.
C
      OPEN(11,FILE='exs.tab')
      WRITE(11,1101)
 1101 FORMAT(1X,'# PENELOPE >>>  Electron cross sections (xs).',/1X,'#',
     1  /1X,'# NOTE: The xs for inner shell ionization (isi) is listed',
     2  ' only for',/1X,'#',7X,'completeness. The xs for inelastic',
     3  ' collisions (in) has been',/1X,'#',7X,'calculated by',
     4  ' considering all inelastic events, including isi.',
     5  /1X,'#',7X,'The bremss (br) xs is for emission of photons with',
     6  ' W >10 eV.')
      WRITE(11,1102)
 1102 FORMAT(1X,'#',/1X,'#   Energy',8X,'xs_el',8X,'xs_in',8X,'xs_br',
     1  7X,'xs_tot',7X,'xs_isi',/1X,'#    (keV)',7X,'(cm**2)',6X,
     2  '(cm**2)',6X,'(cm**2)',6X,'(cm**2)',6X,'(cm**2)',/1X,'#',
     3  78('-'))
C
      OPEN(12,FILE='emfps.tab')
      WRITE(12,1201)
 1201 FORMAT(1X,'# PENELOPE >>>  Electron mean free paths (mfp).',
     1  /1X,'#',/1X,'# NOTE: The xs for inner shell ionization (isi)',
     2  ' is listed only for',/1X,'#',7X,'completeness. The xs for ',
     3  'inelastic collisions (in) has been',/1X,'#',7X,'calculated by',
     4  ' considering all inelastic events, including isi.',
     5  /1X,'#',7X,'The bremss (br) xs is for emission of photons with',
     6  ' W >10 eV.')
      WRITE(12,1202)
 1202 FORMAT(1X,'#',/1X,'#   Energy',7X,'mfp_el',7X,'mfp_in',7X,
     1  'mfp_br',6X,'mfp_tot',6X,'mfp_isi',/1X,'#    (keV)',8X,
     2  '(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',/1X,'#',78('-'))
C
      OPEN(13,FILE='emfprange.tab')
      WRITE(13,1301)
 1301 FORMAT(1X,'# PENELOPE >>>  Electron mean free paths and range.')
      WRITE(13,1302)
 1302 FORMAT(1X,'#',/1X,'#   Energy',7X,'mfp_tot',5X,'1st tmfp',
     1  7X,'range',/1X,'#    (keV)',6X,'(g/cm**2)',4X,'(g/cm**2)',
     1  4X,'(g/cm**2)',/1X,'#',52('-'))
C
      OPEN(14,FILE='estp.tab')
      WRITE(14,1401)
 1401 FORMAT(1X,'# PENELOPE >>>  Electron stopping powers (stp).')
      WRITE(14,1402)
 1402 FORMAT(1X,'#',/1X,'#   Energy',7X,'stp_col',6X,'stp_rad',
     1  6X,'stp_tot',/1X,'#    (keV)',5X,'(eV cm**2/g)',1X,
     2  '(eV cm**2/g)',1X,'(eV cm**2/g)',/1X,'#',52('-'))
C
      DO IE=1,NEGP
        E=ET(IE)
        XSEL=EXP(SEHEL(M,IE))/VMOL(M)
        XSIN=EXP(SEHIN(M,IE))/VMOL(M)
        XSBR=EXP(SEHBR(M,IE))/VMOL(M)
        XSSI=EXP(SEISI(M,IE))/VMOL(M)
        XSTOT=XSEL+XSIN+XSBR
        WRITE(11,'(1X,1P,7E13.5)') E*1.0D-3,XSEL,XSIN,XSBR,XSTOT,XSSI
        FPEL=1.0D0/(XSEL*VMOL(M))
        FPIN=1.0D0/(XSIN*VMOL(M))
        FPBR=1.0D0/(XSBR*VMOL(M))
        FPSI=1.0D0/(XSSI*VMOL(M))
        FPTOT=1.0D0/(XSTOT*VMOL(M))
        WRITE(12,'(1X,1P,7E13.5)') E*1.0D-3,FPEL,FPIN,FPBR,FPTOT,FPSI
        FPTOT=RHO(M)/(XSTOT*VMOL(M))
        FPEL1=RHO(M)/(XE1(IE)*VMOL(M))
        STPC=CSTPE(M,IE)/RHO(M)
        STPR=RSTPE(M,IE)/RHO(M)
        STP=STPC+STPR
        RANGE=PRANGE(E,1,M)*RHO(M)
        WRITE(13,'(1X,1P,7E13.5)') E*1.0D-3,FPTOT,FPEL1,RANGE
        WRITE(14,'(1X,1P,7E13.5)') E*1.0D-3,STPC,STPR,STP
      ENDDO
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
C
C  ****  Positrons.
C
      OPEN(21,FILE='pxs.tab')
      WRITE(21,2101)
 2101 FORMAT(1X,'# PENELOPE >>>  Positron cross sections (xs).',/1X,'#',
     1  /1X,'# NOTE: The xs for inner shell ionization (isi) is listed',
     2  ' only for',/1X,'#',7X,'completeness. The xs for inelastic',
     3  ' collisions (in) has been',/1X,'#',7X,'calculated by',
     4  ' considering all inelastic events, including isi.',
     5  /1X,'#',7X,'The bremss (br) xs is for emission of photons with',
     6  ' W >10 eV.')
      WRITE(21,2102)
 2102 FORMAT(1X,'#',/1X,'#   Energy',8X,'xs_el',8X,'xs_in',8X,'xs_br',
     1  8X,'xs_an',7X,'xs_tot',7X,'xs_isi',/1X,'#    (keV)',7X,
     2  '(cm**2)',6X,'(cm**2)',6X,'(cm**2)',6X,'(cm**2)',6x,'(cm**2)',
     3  6X,'(cm**2)',/1X,'#',91('-'))
C
      OPEN(22,FILE='pmfps.tab')
      WRITE(22,2201)
 2201 FORMAT(1X,'# PENELOPE >>>  Positron mean free paths (mfp).',
     1  /1X,'#',/1X,'# NOTE: The xs for inner shell ionization (isi)',
     2  ' is listed only for',/1X,'#',7X,'completeness. The xs for ',
     3  'inelastic collisions (in) has been',/1X,'#',7X,'calculated by',
     4  ' considering all inelastic events, including isi.',
     5  /1X,'#',7X,'The bremss (br) xs is for emission of photons with',
     6  ' W >10 eV.')
      WRITE(22,2202)
 2202 FORMAT(1X,'#',/1X,'#   Energy',7X,'mfp_el',7X,'mfp_in',7X,
     1  'mfp_br',7X,'mfp_an',7X,'mfp_tot',6X,'mfp_isi',/1X,'#    (keV)',
     2  8X,'(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',9X,'(cm)',
     3  /1X,'#',91('-'))
C
      OPEN(23,FILE='pmfprange.tab')
      WRITE(23,2301)
 2301 FORMAT(1X,'# PENELOPE >>>  Positron mean free paths and range.')
      WRITE(23,2302)
 2302 FORMAT(1X,'#',/1X,'#   Energy',7X,'mfp_tot',5X,'1st tmfp',
     1  7X,'range',/1X,'#    (keV)',6X,'(g/cm**2)',4X,'(g/cm**2)',
     1  4X,'(g/cm**2)',/1X,'#',52('-'))
C
      OPEN(24,FILE='pstp.tab')
      WRITE(24,2401)
 2401 FORMAT(1X,'# PENELOPE >>>  Positron stopping powers (stp).')
      WRITE(24,2402)
 2402 FORMAT(1X,'#',/1X,'#   Energy',7X,'stp_col',6X,'stp_rad',
     1  6X,'stp_tot',/1X,'#    (keV)',5X,'(eV cm**2/g)',1X,
     2  '(eV cm**2/g)',1X,'(eV cm**2/g)',/1X,'#',52('-'))
C
      DO IE=1,NEGP
        E=ET(IE)
        XSEL=EXP(SPHEL(M,IE))/VMOL(M)
        XSIN=EXP(SPHIN(M,IE))/VMOL(M)
        XSBR=EXP(SPHBR(M,IE))/VMOL(M)
        XSSI=EXP(SPISI(M,IE))/VMOL(M)
        XSAN=EXP(SPAN(M,IE))/VMOL(M)
        XSTOT=XSEL+XSIN+XSBR+XSAN
        WRITE(21,'(1X,1P,7E13.5)') E*1.0D-3,XSEL,XSIN,XSBR,XSAN,XSTOT,
     1    XSSI
        FPEL=1.0D0/(XSEL*VMOL(M))
        FPIN=1.0D0/(XSIN*VMOL(M))
        FPBR=1.0D0/(XSBR*VMOL(M))
        FPAN=1.0D0/(XSAN*VMOL(M))
        FPSI=1.0D0/(XSSI*VMOL(M))
        FPTOT=1.0D0/(XSTOT*VMOL(M))
        WRITE(22,'(1X,1P,7E13.5)') E*1.0D-3,FPEL,FPIN,FPBR,FPAN,FPTOT,
     1    FPSI
        FPTOT=RHO(M)/(XSTOT*VMOL(M))
        FPEL1=RHO(M)/(XP1(IE)*VMOL(M))
        STPC=CSTPP(M,IE)/RHO(M)
        STPR=RSTPP(M,IE)/RHO(M)
        STP=STPC+STPR
        RANGE=PRANGE(E,3,M)*RHO(M)
        WRITE(23,'(1X,1P,7E13.5)') E*1.0D-3,FPTOT,FPEL1,RANGE
        WRITE(24,'(1X,1P,7E13.5)') E*1.0D-3,STPC,STPR,STP
      ENDDO
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)
C
      RETURN
      END
