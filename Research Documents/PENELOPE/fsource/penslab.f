      INCLUDE 'penelope.f'  ! Inserted to simplify compilation.
      INCLUDE 'penvared.f'
      INCLUDE 'timer.f'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C       PPPPP   EEEEEE  N    N   SSSS   L         AA    BBBBB          C
C       P    P  E       NN   N  S    S  L        A  A   B    B         C
C       P    P  E       N N  N  SS      L       A    A  B    B         C
C       PPPPP   EEEE    N  N N    SSSS  L       AAAAAA  BBBBB          C
C       P       E       N   NN  S    S  L       A    A  B    B         C
C       P       EEEEEE  N    N   SSSS   LLLLLL  A    A  BBBBB          C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
C  PENELOPE's main program for simulation of electron-photon showers   C
C  in a single slab.                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
C  The material system consists of a homogeneous layer limited by the
C  planes Z=0 and Z=thickness. Primary particles of a given kind, KPARP,
C  are emitted from a point source, either with fixed energy SE0 or with
C  a specified (histogram-like) energy spectrum. The initial direction
C  of the primary particles is sampled uniformly inside a cone of semi-
C  aperture SALPHA and with central axis in the direction (STHETA,SPHI).
C  Thus, SALPHA=0 defines a monodirectional source and SALPHA=180 deg
C  corresponds to an isotropic source.
C
C  ************  The input data file.
C
C  Each line in the input data file consists of a 6-character keyword
C  (columns 1-6) followed either by numerical data (in free format) or
C  by a character string, which start at the 8th column. Keywords are
C  explicitly used/verified by the program (which is case sensitive!).
C  Notice also that the order of the data lines is important. The
C  keyword '______' (6 blanks) indicates comment lines, these can be
C  placed anywhere after the end of the geometry definition list. The
C  program ignores any text following the first blank after the last
C  numerical datum, or after the character string, in each line (thus,
C  in the table given below, the comments in square brackets are
C  ignored by the program). Lines with certain keywords (e.g., 'SPECTR')
C  can appear an arbitrary number of times, limited only by the allo-
C  cated amount of memory.
C
C  The program assigns default values to many input variables; lines
C  that declare default values may be eliminated from the input file.
C  Notice that the default source is a pencil beam that moves upwards
C  along the Z-axis.
C
C  The structure of the input file is the following,
C
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C  TITLE  Title of the job, up to 120 characters.
C
C         >>>>>>>> Source definition.
C  SKPAR  KPARP    [Primary particles: 1=electron, 2=photon, 3=positron]
C  SENERG SE0              [Initial energy (monoenergetic sources only)]
C  SPECTR Ei,Pi                 [E bin: lower-end and total probability]
C  SPOSIT SX0,SY0,SZ0                 [Coordinates of the source centre]
C  SDIREC STHETA,SPHI               [Beam axis direction angles, in deg]
C  SAPERT SALPHA                                 [Beam aperture, in deg]
C
C         >>>>>>>> Material data and simulation parameters.
C  SIMPAR EABS(1:3),C1,C2,WCC,WCR                [Simulation parameters]
C  PFNAME mat_filename.ext     [Material definition file, 20 characters]
C
C         >>>>>>>> Geometry definition.
C  THICKN THICK                                  [Slab thickness, in cm]
C  DSMAX  DSMAX                             [Maximum step length, in cm]
C
C         >>>>>>>> Interaction forcing.
C  IFORCE KPAR,ICOL,FORCER,WLOW,WHIG               [Interaction forcing]
C
C         >>>>>>>> Counter array dimensions and pdf ranges.
C  NBE    EMIN,EMAX,NBE              [Energy interval and no. of E-bins]
C  NBTH   NBTH                   [No. of bins for the polar angle THETA]
C  NBZ    NBZ                         [No. of bins for the Z-coordinate]
C  NBTL   TLMIN,TLMAX,NBTL    [Track-length interval and no. of TL-bins]
C
C         >>>>>>>> Job properties.
C  RESUME dumpfile_1.dat          [Resume from this dump file, 20 chars]
C  DUMPTO dumpfile_2.dat             [Generate this dump file, 20 chars]
C  DUMPP  DUMPP                                 [Dumping period, in sec]
C
C  NSIMSH DSHN                     [Desired number of simulated showers]
C  RSEED  ISEED1,ISEED2           [Seeds of the random number generator]
C  TIME   TIMEA                       [Allotted simulation time, in sec]
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C
C
C  The following listing describes the function of each of the keywords,
C  the accompanying data and their default values. For clarity, blanks
C  in keywords are indicated as '_'.
C
C  TITLE_ : Title of the job (up to 120 characters).
C             DEFAULT: none (the input file must start with this line)
C
C           The TITLE string is used to mark dump files. To prevent the
C           improper use of wrong resuming files, change the title each
C           time you modify basic parameters of your problem. The code
C           will then be able to identify the inconsistency and to print
C           an error message before stopping.
C
C  >>>>>>>> Source definition.
C
C  SKPAR_ : Kind of primary particle (1=electrons, 2=photons or
C           3=positrons).
C             DEFAULT: KPARP=1
C
C  SENERG : For a monoenergetic source, initial energy SE0 of primary
C           particles.
C             DEFAULT: SE0=1.0E6
C
C  SPECTR : For a source with continuous (stepwise constant) spectrum,
C           each 'SPECTR' line gives the lower end-point of an energy
C           bin of the source spectrum (Ei) and the associated relative
C           probability (Pi), integrated over the bin. Up to NSEM=200
C           lines, in arbitrary order. The upper end of the spectrum is
C           defined by entering a line with Ei equal to the upper energy
C           value and with Pi set to a negative value.
C             DEFAULT: none
C
C  SPOSIT : Coordinates of the centre of the source
C             DEFAULT: SX0=SY0=0.0, SZ0=0.0
C
C  SDIREC : Polar and azimuthal angles of the source beam axis direc-
C           tion, in deg.
C             DEFAULTS: STHETA=0.0, SPHI=0.0
C
C           NOTE: The angular distribution of emerging particles depends
C           on both the polar angle THETA and the azimuthal angle PHI.
C           PENSLAB generates only the probability density function of
C           THETA, that is, the angular distribution averaged over PHI.
C           This averaged distribution represents the true distribution
C           only when the source is symmetrical about the Z-axis (i.e.
C           when the incident beam axis is perpendicular to the surface
C           of the slab).
C
C  SAPERT : Angular aperture of the source beam, in deg.
C             DEFAULT: SALPHA=0.0
C
C  >>>>>>>> Material data and simulation parameters.
C
C  SIMPAR : Simulation parameters; absorption energies, EABS(1:3),
C           elastic scattering parameters, C1 and C2, and cutoff energy
C           losses for inelastic collisions and bremsstrahlung emission,
C           WCC and WCR. One line for each material.
C             DEFAULTS: EABS(1,M)=EABS(3,M)=0.01*EPMAX,
C                       EABS(2,M)=0.001*EPMAX
C                       C1(M)=C2(M)=0.1, WCC=EABS(1,M), WCR=EABS(2,M)
C             EPMAX is the maximum energy of all particles found in the
C                   simulation (depends on the source energies).
C
C  PFNAME : Name of PENELOPE's input material data file (20 characters).
C             DEFAULT: none
C
C  >>>>>>>> Geometry definition.
C
C  THICKN : Slab thickness, in cm.
C             DEFAULT: THICK=1.0E0
C
C  DSMAX_ : Maximum step length for electrons and positrons, in cm.
C             DEFAULT: DSMAX=1.0E35
C
C  >>>>>>>> Interaction forcing.
C
C  IFORCE : Activates forcing of interactions of type ICOL of particles
C           KPAR. FORCER is the forcing factor, which must be larger
C           than unity. WLOW and WHIG are the lower and upper limits of
C           the weight window where interaction forcing is applied.
C             DEFAULT: no interaction forcing
C
C           If the mean free path for real interactions of type ICOL is
C           MFP, the program will simulate interactions of this type
C           (real or forced) with an effective mean free path equal to
C           MFP/FORCER.
C
C           TRICK: a negative input value of FORCER, -FN, is  assumed to
C           mean that each particle should interact, on average and
C           approximately, +FN times in a path length equal to the range
C           of that kind of particle with E=EPMAX. This is very useful
C           to generate x-ray spectra from bulk samples.
C
C  The real effect of interaction forcing on the efficiency is not easy
C  to predict. Please, do tentative runs with different FORCER values
C  and check the efficiency gain (or loss!).
C
C  >>>>>>>> Counter array dimensions and tallied ranges.
C
C  NBE___ : Limits of the interval where energy distributions of
C           emerging particles are tallied and number of energy bins.
C             DEFAULT: EMIN=0.0, EMAX=EPMAX, NBE=100
C
C  NBTH__ : Number of bins for the polar angle THETA.
C             DEFAULT: NBTH=90
C
C           WARNING: In the output files, the terms 'transmitted' and
C           'backscattered' are used to denote particles that leave the
C           material system moving upwards (THETA>0) and downwards
C           (THETA<0), respectively. Notice that this agrees with the
C           usual meaning of these terms only when primary particles
C           impinge on the system coming from below (i.e. with THETA>0).
C
C  NBZ___ : Number of bins for the Z-coordinate.
C             DEFAULT: NBZ=100
C
C  NBTL__ : Limits of the interval where track-length distributions of
C           primary particles are tallied. Number of track-length bins.
C             DEFAULT: TLMIN=0, TLMAX=2.5*RANGE(EPMAX,KPARP,1), NBTL=100
C
C  >>>>>>>> Job properties.
C
C  RESUME : The program will read the dump file named 'dumpfile_1.dat'
C           (20 characters) and resume the simulation from the point
C           where it was left. Use this option very, _VERY_ carefully.
C           Make sure that the input data file is fully consistent with
C           the one used to generate the dump file.
C             DEFAULT: off
C
C  DUMPTO : Generate a dump file named 'dumpfile_2.ext' (20 characters)
C           after completing the simulation run. This allows resuming
C           the simulation later on to improve statistics.
C             DEFAULT: off
C
C           NOTE: If the file 'dumpfile_2.ext' already exists, it is
C           overwritten.
C
C  DUMPP_ : When the DUMPTO option is activated, simulation results are
C           written in the output files every DUMPP seconds. This option
C           is useful to check the progress of long simulations. It also
C           allows running the program with a long execution time and
C           stopping it when the required statistical uncertainty has
C           been reached.
C             DEFAULT: DUMPP=1.0D15
C
C  NSIMSH : Desired number of simulated showers.
C             DEFAULT: DSHN=2.0E9
C
C  RSEED_ : Seeds of the random number generator.
C             DEFAULT: ISEED1=12345; ISEED2=54321
C
C  TIME__ : Allotted simulation time, in sec.
C             DEFAULT: TIMEA=100.0E0
C
C  The program is aborted when an incorrect input datum is found. The
C  conflicting quantity usually appears in the last line of the output
C  file. If the trouble is with arrays having dimensions smaller than
C  required, the program indicates how the problem can be solved (this
C  usually requires editing the source file, be careful).
C
C  The clock subroutine (TIMER) may have to be adapted to your specific
C  computer-compiler configuration; standard FORTRAN 77 does not provide
C  timing tools. However, the routines in module TIMER.F do work for
C  many FORTRAN compilers.
C
C  ************  Generating the executable PENSLAB and running it.
C
C  To generate the executable binary file PENSLAB.EXE, compile and link
C  the FORTRAN 77 source files PENSLAB.F, PENELOPE.F, PENVARED.F and
C  TIMER.F. For example, if you are using the g77 compiler under
C  Windows, place these four files in the same directory, open a DOS
C  window and from that directory enter the command
C    `g77 -Wall -O PENSLAB.F -o PENSLAB.EXE'
C  (The same, with file names in lowercase, should work under Linux).
C
C  To run PENSLAB, you have to generate an input data file, let's call
C  it PENSLAB.IN, and the corresponding material data file. Place these
C  two files and the binary file PENSLAB.EXE in the same directory and,
C  from there, issue the command
C    'PENSLAB.EXE < PENSLAB.IN'
C
C  The calculated distributions are written in separate files, whose
C  names start with 'ps_' (for PenSlab) and have the extension '.dat'.
C  These files are in a format suited for direct visualisation with
C  GNUPLOT (version 4.0).
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*2 LIT
      CHARACTER*20 PFILE,PFILED,PFILER
      CHARACTER*23 DATE23
      CHARACTER*120 TITLE,TITLE2,BUFFER
      CHARACTER*6 KWORD,
     1  KWTITL,KWKPAR,KWSENE,KWSPEC, KWSPOS,KWSDIR,KWSAPE,KWSIMP,
     1  KWPFNA,KWTHCK,KWSMAX,KWIFOR, KWNBE ,KWNBTH,KWNBZ ,KWNBTL,
     1  KWRESU,KWDUMP,KWDMPP,KWNSIM, KWRSEE,KWTIME,KWCOMM
      PARAMETER(
     1  KWTITL='TITLE ',KWKPAR='SKPAR ',KWSENE='SENERG',KWSPEC='SPECTR',
     1  KWSPOS='SPOSIT',KWSDIR='SDIREC',KWSAPE='SAPERT',KWSIMP='SIMPAR',
     1  KWPFNA='PFNAME',KWTHCK='THICKN',KWSMAX='DSMAX ',KWIFOR='IFORCE',
     1  KWNBE ='NBE   ',KWNBTH='NBTH  ',KWNBZ ='NBZ   ',KWNBTL='NBTL  ',
     1  KWRESU='RESUME',KWDUMP='DUMPTO',KWDMPP='DUMPP ',KWNSIM='NSIMSH',
     1  KWRSEE='RSEED ',KWTIME='TIME  ',KWCOMM='      ')
C
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI,
     1  RA2DE=180.0D0/PI, DE2RA=PI/180.0D0)
C
C  ************  Main-PENELOPE commons.
C
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      DIMENSION RHOI(MAXMAT)
C
C  ****  Source.
C
C  ****  Source energy spectrum.
      PARAMETER (NSEBM=200)
      DIMENSION ES(NSEBM),PTS(NSEBM),IAS(NSEBM),FS(NSEBM),SHIST(NSEBM)
      DATA SHIST/NSEBM*0.0D0/
C
C  ************  Counters.
C
C  ****  Primary particle counters.
      DIMENSION
     1  PRIM(3),PRIM2(3),DPRIM(3),     ! Numbers of IEXIT particles.
     1  SEC(3,3),SEC2(3,3),DSEC(3,3),  ! Generated secondary particles.
     1  AVTL(3),AVTL2(3),     ! Total track length (in material media).
     1  AVW(2),AVW2(2),DAVW(2),        ! Final polar director cosine.
     1  AVA(2),AVA2(2),DAVA(2),        ! Final polar angle.
     1  AVE(2),AVE2(2),DAVE(2),        ! Final energy.
     1  AVI(8),AVI2(8),DAVI(8)         ! Number of interaction events.
C
      DATA PRIM,PRIM2,SEC,SEC2,AVTL,AVTL2,AVW,AVW2,AVA,AVA2,AVE,AVE2,
     1  AVI,AVI2/58*0.0D0/
C
      DIMENSION WSEC(3,3),WSEC2(3,3),WAVTL(3),WAVTL2(3),WAVW(2),
     1  WAVW2(2),WAVA(2),WAVA2(2),WAVE(2),WAVE2(2),WAVI(8),WAVI2(8)
C
C  ****  Continuous distributions.
C
      PARAMETER (NBEM=200,NBTHM=180,NBZM=200,NBTLM=200)
      DIMENSION BSE(3)
C  ----  Track length distributions.
      PARAMETER (NBTL6=6*NBTLM)
      DIMENSION PDTL(3,NBTLM),PDTL2(3,NBTLM)
      DATA PDTL,PDTL2/NBTL6*0.0D0/
C  ----  Angular distributions of emerging particles.
      PARAMETER (NBTH12=12*NBTHM)
      DIMENSION PDA(3,NBTHM),PDA2(3,NBTHM),PDAP(3,NBTHM),LPDA(3,NBTHM)
      DATA PDA,PDA2,PDAP,LPDA/NBTH12*0.0D0/
C  ----  Energy distributions of emerging particles.
      PARAMETER (NBE12=12*NBEM)
      DIMENSION PDE(3,2,NBEM),PDE2(3,2,NBEM),PDEP(3,2,NBEM),
     1          LPDE(3,2,NBEM)
      DATA PDE,PDE2,PDEP,LPDE/NBE12*0.0D0,NBE12*0.0D0/
C  ----  Depth-distributions of deposited energy (dose) and charge.
      PARAMETER (NDZTT=4*NBZM)
      DIMENSION DOSZ(NBZM),DOSZ2(NBZM),DOSZP(NBZM),LDOSZ(NBZM)
      DATA DOSZ,DOSZ2,DOSZP,LDOSZ/NDZTT*0.0D0/
      DIMENSION CHRZ(NBZM),CHRZ2(NBZM),CHRZP(NBZM),LCHRZ(NBZM)
      DATA CHRZ,CHRZ2,CHRZP,LCHRZ/NDZTT*0.0D0/
C  ----  Deposited energy distributions.
      DIMENSION EDEP(NBEM),EDEP2(NBEM)
      DATA EDEP,EDEP2/NBEM*0.0D0,NBEM*0.0D0/
C
C  ****  Interaction forcing parameters.
C
      PARAMETER (NB=1250)
      DIMENSION IFORCE(NB,3),WLOW(NB,3),WHIG(NB,3)
      COMMON/CFORCE/FORCE(NB,3,8)
C
      EXTERNAL RAND
C
C  ****  Time counter initialisation.
C
      CALL TIME0
C
C  ------------------------  Read input data file.
C
      OPEN(26,FILE='penslab.dat')
      WRITE(26,1100)
 1100 FORMAT(//3X,44('*'),/3X,'**    PROGRAM PENSLAB. Input data ',
     1  'file.   **',/3X,44('*'))
C
      CALL PDATET(DATE23)
      WRITE(26,1101) DATE23
 1101 FORMAT(/3X,'Date and time: ',A23)
C
C  ****  Title.
C
      READ(5,'(A6,1X,A120)') KWORD,TITLE
      WRITE(26,'(/3X,A120)') TITLE
C
C  ************  Source description.
C
      WRITE(26,1200)
 1200 FORMAT(//3X,72('-'),/3X,'>>>>>>  Source description.')
C
   21 CONTINUE
      READ(5,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.KWCOMM) GO TO 21
C
      IF(KWORD.EQ.KWKPAR) THEN
        READ(BUFFER,*) KPARP
   22   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 22
      ELSE
        KPARP=1
      ENDIF
      IF(KPARP.LT.1.OR.KPARP.GT.3) KPARP=1
      IF(KPARP.EQ.1) WRITE(26,1210)
 1210 FORMAT(3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(26,1211)
 1211 FORMAT(3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(26,1212)
 1212 FORMAT(3X,'Primary particles: positrons')
C
C  ****  Initial energy of primary particles.
C
      ISPEC=0
      IF(KWORD.EQ.KWSENE) THEN
        NSEB=1
        READ(BUFFER,*) E0
        WRITE(26,1220) E0
 1220   FORMAT(3X,'Initial energy = ',1P,E13.6,' eV')
   23   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 23
      ELSE IF(KWORD.EQ.KWSPEC) THEN
        ISPEC=1
        NSEB=0
   24   CONTINUE
        NSEB=NSEB+1
        READ(BUFFER,*) ES(NSEB),PTS(NSEB)
   25   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 25
        IF(KWORD.EQ.KWSPEC) GO TO 24
      ELSE
        E0=1.0D6
        WRITE(26,1220) E0
      ENDIF
      IF(ISPEC.EQ.1) THEN
        IF(NSEB.GT.1) THEN
          CALL SORT2(ES,PTS,NSEB)
          WRITE(26,1221)
 1221     FORMAT(/3X,'Spectrum:',7X,'I',4X,'E_low(eV)',4x,'E_high(eV)',
     1      5X,'P_sum(E)',/16X,45('-'))
          DO I=1,NSEB-1
            WRITE(26,'(16X,I4,1P,3E14.6)') I,ES(I),ES(I+1),PTS(I)
          ENDDO
          WRITE(26,*) '  '
          E0=ES(NSEB)
          NSEB=NSEB-1
          CALL IRND0(PTS,FS,IAS,NSEB)
        ELSE
          WRITE(26,*) 'The source energy spectrum is not defined.'
          STOP 'The source energy spectrum is not defined.'
        ENDIF
      ENDIF
      IF(E0.LT.100.0D0) THEN
        WRITE(26,*) 'The initial energy E0 is too small.'
        STOP 'The initial energy E0 is too small.'
      ENDIF
      EPMAX=E0
C  ----  Positrons eventually give annihilation gamma-rays. The maximum
C        energy of annihilation photons is .lt. 1.21*(E0+me*c**2).
      IF(KPARP.EQ.3) EPMAX=1.21D0*(E0+5.12D5)
C
C  ****  Source position.
C
      IF(KWORD.EQ.KWSPOS) THEN
        READ(BUFFER,*) SX0,SY0,SZ0
   28   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 28
      ELSE
        SX0=0.0D0
        SY0=0.0D0
        SZ0=-1.0D15
      ENDIF
      WRITE(26,1232) SX0,SY0,SZ0
 1232 FORMAT(3X,'Source coordinates:        SX0 = ',1P,E13.6,
     1  ' cm',/30X,'SY0 = ',E13.6,' cm',/30X,'SZ0 = ',E13.6,' cm')
C
C  ****  Source beam direction and aperture.
C
      IF(KWORD.EQ.KWSDIR) THEN
        READ(BUFFER,*) STHETA,SPHI
   29   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 29
      ELSE
        STHETA=0.0D0
        SPHI=0.0D0
      ENDIF
      WRITE(26,1233) STHETA,SPHI
 1233 FORMAT(3X,'Beam direction angles:   THETA = ',1P,E13.6,' deg',/
     1  30X,'PHI = ',E13.6,' deg')
C
      IF(KWORD.EQ.KWSAPE) THEN
        READ(BUFFER,*) SALPHA
   30   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 30
      ELSE
        SALPHA=0.0D0
      ENDIF
      WRITE(26,1234) SALPHA
 1234 FORMAT(3X,'Beam aperture:',11X,'ALPHA = ',1P,E13.6,' deg')
      CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
C
C  ************  Material data and simulation parameters.
C
      WRITE(26,1300)
 1300 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Material data and simulation parameters.')
C
C  ****  Simulation parameters.
C
      NMAT=1
      DO M=1,NMAT
        EABS(1,M)=0.010D0*EPMAX
        EABS(2,M)=0.001D0*EPMAX
        EABS(3,M)=0.010D0*EPMAX
        C1(M)=0.10D0
        C2(M)=0.10D0
        WCC(M)=EABS(1,M)
        WCR(M)=EABS(2,M)
      ENDDO
C
      IF(KWORD.EQ.KWSIMP) THEN
        READ(BUFFER,*) EABS(1,1),EABS(2,1),EABS(3,1),C1(1),C2(1),
     1    WCC(1),WCR(1)
   32   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 32
      ENDIF
C
      DO M=1,NMAT
        IF(M.EQ.1) LIT='st'
        IF(M.EQ.2) LIT='nd'
        IF(M.EQ.3) LIT='rd'
        IF(M.GT.3) LIT='th'
        WRITE(26,1320) M,LIT
 1320   FORMAT(/3X,'**** ',I2,A2,' material')
        IF(EABS(1,M).LT.1.0D2) EABS(1,M)=1.0D2
        IF(EABS(2,M).LT.1.0D2) EABS(2,M)=1.0D2
        IF(EABS(3,M).LT.1.0D2) EABS(3,M)=1.0D2
        WRITE(26,1321) EABS(1,M)
 1321 FORMAT(3X,'Electron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1322) EABS(2,M)
 1322 FORMAT(3X,'  Photon absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1323) EABS(3,M)
 1323 FORMAT(3X,'Positron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1324) C1(M),C2(M),WCC(M),WCR(M)
 1324 FORMAT(3X,'Electron-positron simulation parameters:',
     1 /4X,'C1 =',1P,E13.6,',      C2 =',E13.6,/3X,'Wcc =',E13.6,
     1  ' eV,  Wcr =',E13.6,' eV')
      ENDDO
C
C  ****  Initialisation of PENELOPE.
C
      IF(KWORD.EQ.KWPFNA) THEN
        READ(BUFFER,'(A20)') PFILE
        WRITE(26,1330) PFILE
 1330   FORMAT(/3X,'PENELOPE''s material definition file: ',A20)
   33   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 33
      ELSE
        WRITE(26,*) 'You have to specify a material file.'
        STOP 'You have to specify a material file.'
      ENDIF
C
      OPEN(15,FILE=PFILE)
      OPEN(16,FILE='material.dat')
        INFO=3
        CALL PEINIT(EPMAX,NMAT,15,16,INFO)
      CLOSE(UNIT=15)
      CLOSE(UNIT=16)
C  ----  Inverse densities are used to score the local absorbed dose.
      DO M=1,NMAT
        RHOI(M)=1.0D0/RHO(M)
      ENDDO
C
C  ****  Slab thickness.
C
      IF(KWORD.EQ.KWTHCK) THEN
        READ(BUFFER,*) THICK
   34   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 34
        WRITE(26,1080) THICK
 1080   FORMAT(/3X,'Slab thickness = ',1P,E13.6,' cm')
        IF(THICK.LT.0.0D0) STOP
        IF(SZ0.LT.0.0D0) SZ0=0.0D0
        IF(SZ0.GT.THICK) SZ0=THICK
      ELSE
        THICK=1.0D0
      ENDIF
C
C  ****  Maximum step length of electrons and positrons.
C
      DSMAX=1.0D35
      IF(KWORD.EQ.KWSMAX) THEN
        READ(BUFFER,*) DSMAX
        IF(DSMAX.LT.1.0D-7) DSMAX=1.0D20
   35   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 35
      ENDIF
      WRITE(26,1351) DSMAX
 1351 FORMAT(3X,'Maximum step length, DSMAX = ',1P,E13.6,' cm')
C
C  ************  Variance reduction (only interaction forcing, for the
C                time being...).
C
      DO KB=1,NB
        DO ICOL=1,8
          DO KPAR=1,3
            FORCE(KB,KPAR,ICOL)=1.0D0
          ENDDO
        ENDDO
        DO KPAR=1,3
          IFORCE(KB,KPAR)=0
          WLOW(KB,KPAR)=0.0D0
          WHIG(KB,KPAR)=1.0D6
        ENDDO
      ENDDO
C
      IF(KWORD.EQ.KWIFOR) THEN
        WRITE(26,1600)
 1600   FORMAT(//3X,72('-'),/
     1    3X,'>>>>>>  Interaction forcing: FORCE(KPAR,ICOL)')
   61   CONTINUE
        READ(BUFFER,*) KPAR,ICOL,FORCER,WWLOW,WWHIG
C  ****  Negative FORCER values are re-interpreted, as described in the
C        heading comments above.
        IF(FORCER.LT.-1.0D-6) THEN
          MM=1
          EVENTS=MAX(ABS(FORCER),1.0D0)
          PLT=PRANGE(E0,KPAR,MM)
          RMFP=PHMFP(E0,KPAR,MM,ICOL)
          FORCER=EVENTS*RMFP/PLT
        ENDIF
        IF(WWLOW.LT.1.0D-6) WWLOW=1.0D-6
        IF(WWHIG.GT.1.0D6) WWHIG=1.0D6
        IF(KPAR.LT.1.OR.KPAR.GT.3) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of KPAR.'
          STOP 'Incorrect value of KPAR.'
        ENDIF
        IF(ICOL.LT.1.OR.ICOL.GT.8) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of ICOL.'
          STOP 'Incorrect value of ICOL.'
        ENDIF
        KB=1
        WLOW(KB,KPAR)=MAX(WLOW(KB,KPAR),WWLOW)
        WHIG(KB,KPAR)=MIN(WHIG(KB,KPAR),WWHIG)
        IF(WLOW(KB,KPAR).GT.WHIG(KB,KPAR)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect weight window limits.'
          STOP 'Incorrect weight window limits.'
        ENDIF
        IF(FORCER.LT.1.0D0) STOP 'FORCER must be greater than unity.'
        IFORCE(KB,KPAR)=1
        FORCE(KB,KPAR,ICOL)=FORCER
        WRITE(26,1610) KPAR,ICOL,FORCER,WLOW(KB,KPAR),WHIG(KB,KPAR)
 1610   FORMAT(3X,'FORCE(',I1,',',I1,') =',1P,E13.6,
     1    ',  weight window = (',E9.2,',',E9.2,')')
   62   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 62
        IF(KWORD.EQ.KWIFOR) GO TO 61
      ENDIF
C
C  ************  Tallied distributions.
C
      WRITE(26,1400)
 1400 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Dimensions of counter arrays.')
C
      IF(KWORD.EQ.KWNBE) THEN
        READ(BUFFER,*) EMIN,EMAX,NBE
   41   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 41
      ELSE
        NBE=NBEM
        EMIN=0.0D0
        EMAX=EPMAX
      ENDIF
      EMIN=MAX(EMIN,0.0D0)
      EMAX=MAX(EMAX,EPMAX)
      WRITE(26,1410) NBE,EMIN,EMAX
 1410 FORMAT(3X,'E:       NBE = ',I3,
     1  ',  EMIN =',1P,E13.6,' eV,  EMAX =',E13.6,' eV')
      IF(NBE.LT.1) THEN
        WRITE(26,*) 'Wrong number of energy bins.'
        STOP 'Wrong number of energy bins.'
      ELSE IF (NBE.GT.NBEM) THEN
        WRITE(26,*) 'NBE is too large.'
        WRITE(26,*) 'Set the parameter NBEM equal to ',NBE
        STOP 'NBE is too large.'
      ENDIF
      IF(EMIN.GT.EMAX) THEN
        WRITE(26,*) 'Energy interval is too narrow.'
        STOP 'Energy interval is too narrow.'
      ENDIF
C
      IF(KWORD.EQ.KWNBTH) THEN
        READ(BUFFER,*) NBTH
   42   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 42
      ELSE
        NBTH=NBTHM
      ENDIF
      WRITE(26,1420) NBTH
 1420 FORMAT(3X,'Theta:  NBTH = ',I3)
      IF(NBTH.LT.1) THEN
        WRITE(26,*) 'Wrong number of THETA bins.'
        STOP 'Wrong number of THETA bins.'
      ELSE IF (NBTH.GT.NBTHM) THEN
        WRITE(26,*) 'NBTH is too large.'
        WRITE(26,*) 'Set the parameter NBTHM equal to ',NBTH
        STOP 'NBTH is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBZ) THEN
        READ(BUFFER,*) NBZ
   44   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 44
      ELSE
        NBZ=NBZM
      ENDIF
      WRITE(26,1440) NBZ
 1440 FORMAT(3X,'Z:       NBZ = ',I3)
      IF(NBZ.LT.1) THEN
        WRITE(26,*) 'Wrong number of Z bins.'
        STOP 'Wrong number of Z bins.'
      ELSE IF (NBZ.GT.NBZM) THEN
        WRITE(26,*) 'NBZ is too large.'
        WRITE(26,*) 'Set the parameter NBZM equal to ',NBZ
        STOP 'NBZ is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBTL) THEN
        READ(BUFFER,*) TLMIN,TLMAX,NBTL
   46   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 46
      ELSE
        NBTL=NBTLM
        TLMIN=0.0D0
        TLMAX=2.5D0*PRANGE(EPMAX,KPARP,1)
      ENDIF
      TLMIN=MAX(TLMIN,0.0D0)
      TLMAX=MAX(TLMAX,2.5D0*PRANGE(EPMAX,KPARP,1))
      WRITE(26,1460) NBTL,TLMIN,TLMAX
 1460 FORMAT(3X,'TL:     NBTL = ',I3,
     1  ', TLMIN =',1P,E13.6,' cm, TLMAX =',E13.6,' cm')
      IF(NBTL.LT.1) THEN
        WRITE(26,*) 'Wrong number of track length bins.'
        STOP 'Wrong number of track length bins.'
      ENDIF
      IF(TLMIN.GE.TLMAX) THEN
        WRITE(26,*) 'Track length interval is too narrow.'
        STOP 'Track length interval is too narrow.'
      ENDIF
C
C  ****  Bin sizes.
C
      BSDE=EMAX/DBLE(NBE)
      BSE(1)=(EMAX-EMIN)/DBLE(NBE)
      BSE(2)=(EMAX-EMIN)/DBLE(NBE)
      BSE(3)=(EMAX-EMIN)/DBLE(NBE)
      BSTH=180.0D0/DBLE(NBTH)
      BSTL=(TLMAX-TLMIN)/DBLE(NBTL)
      BSZ=THICK/DBLE(NBZ)
C  ----  The factor 1.0000001 serves to ensure that the upper limit of
C  the tallied interval is within the last channel (otherwise, the array
C  dimensions could be exceeded).
      BSDE=1.0000001D0*BSDE
      BSE(1)=1.0000001D0*BSE(1)
      BSE(2)=1.0000001D0*BSE(2)
      BSE(3)=1.0000001D0*BSE(3)
      BSTH=1.0000001D0*BSTH
      BSTL=1.0000001D0*BSTL
      BSZ=1.0000001D0*BSZ
C
C  ************  Tallied distributions (selected by the user).
C
C  ************  Job characteristics.
C
      WRITE(26,1700)
 1700 FORMAT(//3X,72('-'),/
     1  3X,'>>>>>>  Job characteristics.')
C
      IRESUM=0
      IF(KWORD.EQ.KWRESU) THEN
        READ(BUFFER,'(A20)') PFILER
        WRITE(26,1710) PFILER
 1710   FORMAT(3X,'Resume simulation from previous dump file: ',A20)
        IRESUM=1
   71   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 71
      ENDIF
C
      IDUMP=0
      DUMPP=1.0D15
      IF(KWORD.EQ.KWDUMP) THEN
        READ(BUFFER,'(A20)') PFILED
        WRITE(26,1720) PFILED
 1720   FORMAT(3X,'Write final counter values on the dump file: ',A20)
        IDUMP=1
   72   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 72
      ENDIF
C
      IF(KWORD.EQ.KWDMPP) THEN
        READ(BUFFER,*) DUMPP
        IF(IDUMP.EQ.1) THEN
          IF(DUMPP.LT.15.0D0) DUMPP=15.0D0
          IF(DUMPP.GT.86400.0D0) DUMPP=86400.0D0
          WRITE(26,1730) DUMPP
 1730     FORMAT(3X,'Dumping period: DUMPP =',1P,E13.6)
        ENDIF
   73   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 73
      ENDIF
C
      IF(KWORD.EQ.KWNSIM) THEN
        READ(BUFFER,*) DSHN
        IF(DSHN.LT.1.0D0) DSHN=2.0D9
   74   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 74
      ELSE
        DSHN=2.0D9
      ENDIF
      WRITE(26,1740) DSHN
 1740 FORMAT(/3X,'Number of showers to be simulated =',1P,E13.6)
C
      IF(KWORD.EQ.KWRSEE) THEN
        READ(BUFFER,*) ISEED1,ISEED2
        WRITE(26,1750) ISEED1,ISEED2
 1750   FORMAT(3X,'Random number generator seeds = ',I10,', ',I10)
   75   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 75
      ELSE
        ISEED1=12345
        ISEED2=54321
      ENDIF
C
      IF(KWORD.EQ.KWTIME) THEN
        READ(BUFFER,*) TIMEA
      ELSE
        TIMEA=100.0D0
      ENDIF
      IF(TIMEA.LT.1.0D0) TIMEA=100.0D0
      WRITE(26,1760) TIMEA
 1760 FORMAT(3X,'Computation time available = ',1P,E13.6,' sec')
C
      CALL TIMER(TSECIN)
      CPUT0=CPUTIM()
      TSECA=TIMEA+TSECIN
      TSECAD=TSECIN
      WRITE(26,1770)
 1770 FORMAT(//3X,72('-'))
C
C  ************  If 'RESUME' is active, read previously generated
C                counters...
C
      SHNA=0.0D0
      CPUTA=0.0D0
      IRETRN=0
      TDEP=0.0D0
      TDEP2=0.0D0
C
      IF(IRESUM.EQ.1) THEN
        OPEN(9,FILE=PFILER)
        READ (9,*,ERR=1800,END=1800) SHNA,CPUTA
        READ (9,'(A120)') TITLE2
        IF(TITLE2.NE.TITLE) THEN
          WRITE(26,*)
     1      'The dump file is corrupted (the TITLE does not match).'
          STOP 'The dump file is corrupted (the TITLE does not match).'
        ENDIF
        READ (9,*) ISEED1,ISEED2
        READ (9,*) (SHIST(I),I=1,NSEB)
        READ (9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        READ (9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        READ (9,*) (AVTL(I),I=1,3),(AVTL2(I),I=1,3)
        READ (9,*) (AVW(I),I=1,2),(AVW2(I),I=1,2)
        READ (9,*) (AVA(I),I=1,2),(AVA2(I),I=1,2)
        READ (9,*) (AVE(I),I=1,2),(AVE2(I),I=1,2)
        READ (9,*) (AVI(I),I=1,8),(AVI2(I),I=1,8)
        READ (9,*) ((PDTL(I,J),J=1,NBTL),I=1,3),
     1             ((PDTL2(I,J),J=1,NBTL),I=1,3)
        READ (9,*) ((PDA(I,J),J=1,NBTH),I=1,3),
     1             ((PDA2(I,J),J=1,NBTH),I=1,3)
        READ (9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        READ (9,*) TDEP,TDEP2
C
        READ (9,*) (DOSZ(J),J=1,NBZ),(DOSZ2(J),J=1,NBZ)
        READ (9,*) (CHRZ(J),J=1,NBZ),(CHRZ2(J),J=1,NBZ)
        READ (9,*) (EDEP(J),J=1,NBE),(EDEP2(J),J=1,NBE)
        CLOSE(9)
        GO TO 1802
 1800   CONTINUE
        WRITE(26,1801)
 1801   FORMAT(/3X,'WARNING: Could not resume from dump file...',/)
      ENDIF
 1802 CONTINUE
C
C  ************  Initialise constants.
C
      WGHT0=1.0D0   ! Primary particle weight.
      INTFOR=0      ! No interaction forcing as default.
      IBODY=1
C
      SHN=SHNA      ! Shower counter, including the dump file.
      DSHN=DSHN+SHNA
      N=MOD(SHN,2.0D9)+0.5D0
      IF(SHN.GE.DSHN) GO TO 105
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  ------------------------  Shower simulation starts here.
C
  101 CONTINUE  ! The simulation loop starts here.
C
C  **********  We set the initial state of the primary particle.
C
      SHN=SHN+1.0D0
      N=N+1
      IF(N.GT.2000000000) N=N-2000000000
      KPAR=KPARP
      WGHT=WGHT0
      MAT=1
C  ----  Initial direction ...
      CALL GCONE(U,V,W)
C  ----  Initial position ...
      X=SX0
      Y=SY0
      IF(SZ0.LT.0.0D0) THEN
        IF(W.GT.0.0D0) THEN
          Z=0.0D0
        ELSE
          GO TO 105
        ENDIF
      ELSE IF(SZ0.LT.THICK) THEN
        Z=SZ0
      ELSE
        IF(W.LT.0.0D0) THEN
          Z=THICK
        ELSE
          GO TO 105
        ENDIF
      ENDIF
C  ----  Initial energy ...
      IF(ISPEC.EQ.0) THEN
        E=E0    ! Monoenergetic source.
        SHIST(1)=SHIST(1)+1.0D0
      ELSE      ! Continuous spectrum. E sampled by Walker method.
        RN=RAND(4.0D0)*NSEB+1
        K=INT(RN)
        RNF=RN-K
        IF(RNF.GT.FS(K)) THEN
          KE=IAS(K)
        ELSE
          KE=K
        ENDIF
        E=ES(KE)+RAND(5.0D0)*(ES(KE+1)-ES(KE))
        SHIST(KE)=SHIST(KE)+1.0D0
      ENDIF
      IF(E.LT.EABS(KPAR,MAT)) GO TO 105
C
C  ************  Initialisation of primary particle counters.
C
      ILB(1)=1          ! Identifies primary particles.
      ILB(2)=0
      ILB(3)=0
      ILB(4)=0
      ILB(5)=0
      DO I=1,3
        DPRIM(I)=0.0D0
        DO K=1,3
          DSEC(K,I)=0.0D0
        ENDDO
      ENDDO
      DO I=1,8
        DAVI(I)=0.0D0   ! Numbers of collisions.
      ENDDO
      DO I=1,2
        DAVW(I)=0.0D0
        DAVA(I)=0.0D0
        DAVE(I)=0.0D0
      ENDDO
      DEPP=E0*WGHT      ! Deposited energy within the slab.
      TL=0.0D0          ! Accumulated track length.
C  ---------------------------------------------------------------------
C  ------------------------  Track simulation begins here.
C
C     EBAL=E*WGHT       ! Used to verify that energy is conserved.
C     IF(KPARP.EQ.3) EBAL=EBAL+TREV*WGHT
C
      CALL CLEANS       ! Cleans secondary stack.
  102 CONTINUE
      CALL START        ! Starts simulation in current medium.
C
C  ----  Free path length to the next interaction event.
C
  103 CONTINUE
C     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,i3)')
C    1    mod(n,100),kpar,ilb(1),x,y,z,w,e,ibody
      IF((IFORCE(IBODY,KPAR).EQ.1).AND.
     1   ((WGHT.GT.WLOW(IBODY,KPAR)).AND.
     1    (WGHT.LT.WHIG(IBODY,KPAR)))) THEN
        CALL JUMPF(DSMAX,DS)  ! Interaction forcing.
        INTFOR=1
      ELSE
        CALL JUMP(DSMAX,DS)   ! Analogue simulation.
        INTFOR=0
      ENDIF
C  ****  New position.
      X=X+U*DS
      Y=Y+V*DS
      Z=Z+W*DS
C
      TL=TL+DS   ! Accumulated track length.
C  ****  Transmitted particle.
      IF(Z.GE.THICK) THEN
        IEXIT=1  ! Labels transmitted particles.
        DSCOR=(Z-THICK)/(ABS(W)+1.0D-30)
        Z=THICK
        X=X-U*DSCOR
        Y=Y-V*DSCOR
        TL=TL-DSCOR
        DEPP=DEPP-E*WGHT
        GO TO 104
      ENDIF
C  ****  Backscattered particle.
      IF(Z.LE.0.0D0) THEN
        IEXIT=2  ! Labels backscattered particles.
        DSCOR=ABS(Z)/(ABS(W)+1.0D-30)
        Z=0.0D0
        X=X-U*DSCOR
        Y=Y-V*DSCOR
        TL=TL-DSCOR
        DEPP=DEPP-E*WGHT
        GO TO 104
      ENDIF
C  ----  Simulate next event.
      IF(INTFOR.EQ.0) THEN
        CALL KNOCK(DE,ICOL)   ! Analogue simulation.
      ELSE
        CALL KNOCKF(DE,ICOL)  ! Interaction forcing is active.
      ENDIF
      IF(ILB(1).EQ.1) DAVI(ICOL)=DAVI(ICOL)+WGHT
C  ----  Energy is locally deposited in the material.
      DEP=DE*WGHT
C
C  ****  Depth distributions of dose and deposited charge.
C
      KZ=1.0D0+Z/BSZ          ! Depth channel number.
      IF(KZ.GT.0.AND.KZ.LE.NBZ) THEN
        IF(N.NE.LDOSZ(KZ)) THEN
          DOSZ(KZ)=DOSZ(KZ)+DOSZP(KZ)
          DOSZ2(KZ)=DOSZ2(KZ)+DOSZP(KZ)**2
          DOSZP(KZ)=DEP*RHOI(MAT)
          LDOSZ(KZ)=N
        ELSE
          DOSZP(KZ)=DOSZP(KZ)+DEP*RHOI(MAT)
        ENDIF
C
        IF(E.LT.EABS(KPAR,MAT).AND.KPAR.NE.2) THEN
          IF(N.NE.LCHRZ(KZ)) THEN
            CHRZ(KZ)=CHRZ(KZ)+CHRZP(KZ)
            CHRZ2(KZ)=CHRZ2(KZ)+CHRZP(KZ)**2
            CHRZP(KZ)=(KPAR-2)*WGHT
            LCHRZ(KZ)=N
          ELSE
            CHRZP(KZ)=CHRZP(KZ)+(KPAR-2)*WGHT
          ENDIF
        ENDIF
      ENDIF
C
C  ----  Check if the particle has been absorbed.
      IF(E.GT.EABS(KPAR,MAT)) GO TO 103
      IEXIT=3  ! Labels absorbed particles.
C  ------------------------  The simulation of the track ends here.
C  ---------------------------------------------------------------------
  104 CONTINUE
C
C     EBAL=EBAL-E*WGHT
C     IF(KPAR.EQ.3.AND.IEXIT.NE.3) EBAL=EBAL-TREV*WGHT
C
C  ************  Counters.
C
C  ****  Primary particles.
C
      THETA=ACOS(W)
      IF(ILB(1).EQ.1) THEN
        DPRIM(IEXIT)=DPRIM(IEXIT)+WGHT0
        AVTL(IEXIT)=AVTL(IEXIT)+TL*WGHT0
        AVTL2(IEXIT)=AVTL2(IEXIT)+(TL*WGHT0)**2
        IF(IEXIT.LT.3) THEN
          DAVW(IEXIT)=DAVW(IEXIT)+W*WGHT0
          DAVA(IEXIT)=DAVA(IEXIT)+THETA*WGHT0
          DAVE(IEXIT)=DAVE(IEXIT)+E*WGHT0
        ENDIF
C  ****  Track-length distribution.
        KT=1.0D0+(TL-TLMIN)/BSTL
        IF(KT.GT.0.AND.KT.LE.NBTL) THEN
          PDTL(IEXIT,KT)=PDTL(IEXIT,KT)+WGHT0
          PDTL2(IEXIT,KT)=PDTL2(IEXIT,KT)+WGHT0**2
        ENDIF
      ELSE
        DSEC(KPAR,IEXIT)=DSEC(KPAR,IEXIT)+WGHT
      ENDIF
C
C  ****  Other counters.
C
      IF(IEXIT.LT.3) THEN
C  ****  Angular distribution of emerging particles.
        KTH=1.0D0+THETA*RA2DE/BSTH
        IF(N.NE.LPDA(KPAR,KTH)) THEN
          PDA(KPAR,KTH)=PDA(KPAR,KTH)+PDAP(KPAR,KTH)
          PDA2(KPAR,KTH)=PDA2(KPAR,KTH)+PDAP(KPAR,KTH)**2
          PDAP(KPAR,KTH)=WGHT
          LPDA(KPAR,KTH)=N
        ELSE
          PDAP(KPAR,KTH)=PDAP(KPAR,KTH)+WGHT
        ENDIF
C  ****  Energy distribution of emerging particles.
        K=1.0D0+(E-EMIN)/BSE(KPAR)
        IF(K.GT.0.AND.K.LE.NBE) THEN
          IF(N.NE.LPDE(KPAR,IEXIT,K)) THEN
            PDE(KPAR,IEXIT,K)=PDE(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)
            PDE2(KPAR,IEXIT,K)=PDE2(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)**2
            PDEP(KPAR,IEXIT,K)=WGHT
            LPDE(KPAR,IEXIT,K)=N
          ELSE
            PDEP(KPAR,IEXIT,K)=PDEP(KPAR,IEXIT,K)+WGHT
          ENDIF
        ENDIF
      ENDIF
C
C  ************  Any secondary left?
C
      CALL SECPAR(LEFT)
      IF(LEFT.GT.0) THEN
C     write(6,'(/''new secondary'')')
C     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,i3)')
C    1    mod(n,100),kpar,ilb(1),x,y,z,w,e,ibody
C  ****  Subtract E and charge from the tallied distributions to avoid
C        double-counting.
        KZ=1.0D0+Z/BSZ          ! Depth channel number.
        IF(N.NE.LDOSZ(KZ)) THEN
          DOSZ(KZ)=DOSZ(KZ)+DOSZP(KZ)
          DOSZ2(KZ)=DOSZ2(KZ)+DOSZP(KZ)**2
          DOSZP(KZ)=-E*WGHT*RHOI(MAT)
          LDOSZ(KZ)=N
        ELSE
          DOSZP(KZ)=DOSZP(KZ)-E*WGHT*RHOI(MAT)
        ENDIF
C
        IF(N.NE.LCHRZ(KZ)) THEN
          CHRZ(KZ)=CHRZ(KZ)+CHRZP(KZ)
          CHRZ2(KZ)=CHRZ2(KZ)+CHRZP(KZ)**2
          CHRZP(KZ)=(2-KPAR)*WGHT
          LCHRZ(KZ)=N
        ELSE
          CHRZP(KZ)=CHRZP(KZ)+(2-KPAR)*WGHT
        ENDIF
C
        GO TO 102
      ENDIF
C
C  ------------------------  The simulation of the shower ends here.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
        PRIM2(I)=PRIM2(I)+DPRIM(I)**2
        DO K=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)**2
        ENDDO
      ENDDO
      DO I=1,2
        AVW(I)=AVW(I)+DAVW(I)
        AVW2(I)=AVW2(I)+DAVW(I)**2
        AVA(I)=AVA(I)+DAVA(I)
        AVA2(I)=AVA2(I)+DAVA(I)**2
        AVE(I)=AVE(I)+DAVE(I)
        AVE2(I)=AVE2(I)+DAVE(I)**2
      ENDDO
      DO I=1,8
        AVI(I)=AVI(I)+DAVI(I)
        AVI2(I)=AVI2(I)+DAVI(I)**2
      ENDDO
C
      TDEP=TDEP+DEPP
      TDEP2=TDEP2+DEPP**2
      IF(DEPP.GT.1.0D0) THEN
        KE=1.0D0+DEPP/BSDE
        IF(KE.GT.0.AND.KE.LE.NBEM) THEN
          EDEP(KE)=EDEP(KE)+WGHT0
          EDEP2(KE)=EDEP2(KE)+WGHT0**2
        ENDIF
      ENDIF
C
C     IF(ABS(EBAL).GT.0.01D0) THEN
C       WRITE(26,'(//3X,''Energy is not conserved:'')')
C       WRITE(26,'(5X,''EBAL ='',1P,E13.6,'' eV,  N ='',I9)') EBAL,N
C       WRITE(26,'(3X,''Correct the error, or switch off the test.'')')
C       STOP 'Verify energy conservation.'
C     ENDIF
C
  105 CONTINUE
      CALL TIMER(TSEC)
C
C  ----  End the simulation after the allotted time or after completing
C        DSHN showers.
C
      IRETRN=0
      IF(TSEC.LT.TSECA.AND.SHN.LT.DSHN) THEN
C  ****  Write partial results after each dumping period.
        IF(TSEC-TSECAD.GT.DUMPP) THEN
          IF(IDUMP.EQ.1) THEN
            TSECAD=TSEC
            IRETRN=1
            GO TO 203
          ENDIF
        ENDIF
        GO TO 101
      ENDIF
  203 CONTINUE
      TSIM=CPUTA+CPUTIM()-CPUT0
C
C  ****  Transfer contents of partial counters to global counters.
C
      DO KPAR=1,3
        DO K=1,NBTH
          PDA(KPAR,K)=PDA(KPAR,K)+PDAP(KPAR,K)
          PDA2(KPAR,K)=PDA2(KPAR,K)+PDAP(KPAR,K)**2
          PDAP(KPAR,K)=0.0D0
          LPDA(KPAR,K)=0
        ENDDO
        DO K=1,NBE
          PDE(KPAR,1,K)=PDE(KPAR,1,K)+PDEP(KPAR,1,K)
          PDE2(KPAR,1,K)=PDE2(KPAR,1,K)+PDEP(KPAR,1,K)**2
          PDEP(KPAR,1,K)=0.0D0
          LPDE(KPAR,1,K)=0
          PDE(KPAR,2,K)=PDE(KPAR,2,K)+PDEP(KPAR,2,K)
          PDE2(KPAR,2,K)=PDE2(KPAR,2,K)+PDEP(KPAR,2,K)**2
          PDEP(KPAR,2,K)=0.0D0
          LPDE(KPAR,2,K)=0
        ENDDO
      ENDDO
C
      DO J=1,NBZ
        DOSZ(J)=DOSZ(J)+DOSZP(J)
        DOSZ2(J)=DOSZ2(J)+DOSZP(J)**2
        DOSZP(J)=0.0D0
        LDOSZ(J)=0
        CHRZ(J)=CHRZ(J)+CHRZP(J)
        CHRZ2(J)=CHRZ2(J)+CHRZP(J)**2
        CHRZP(J)=0.0D0
        LCHRZ(J)=0
      ENDDO
C
C  ************  If 'DUMPTO' is active, write counters in a dump file.
C
      IF(IDUMP.EQ.1) THEN
        OPEN(9,FILE=PFILED)
        WRITE(9,*) SHN,TSIM
        WRITE(9,'(A120)') TITLE
        WRITE(9,*) ISEED1,ISEED2
        WRITE(9,*) (SHIST(I),I=1,NSEB)
        WRITE(9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        WRITE(9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        WRITE(9,*) (AVTL(I),I=1,3),(AVTL2(I),I=1,3)
        WRITE(9,*) (AVW(I),I=1,2),(AVW2(I),I=1,2)
        WRITE(9,*) (AVA(I),I=1,2),(AVA2(I),I=1,2)
        WRITE(9,*) (AVE(I),I=1,2),(AVE2(I),I=1,2)
        WRITE(9,*) (AVI(I),I=1,8),(AVI2(I),I=1,8)
        WRITE(9,*) ((PDTL(I,J),J=1,NBTL),I=1,3),
     1             ((PDTL2(I,J),J=1,NBTL),I=1,3)
        WRITE(9,*) ((PDA(I,J),J=1,NBTH),I=1,3),
     1             ((PDA2(I,J),J=1,NBTH),I=1,3)
        WRITE(9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        WRITE(9,*) TDEP,TDEP2
C
        WRITE(9,*) (DOSZ(J),J=1,NBZ),(DOSZ2(J),J=1,NBZ)
        WRITE(9,*) (CHRZ(J),J=1,NBZ),(CHRZ2(J),J=1,NBZ)
        WRITE(9,*) (EDEP(J),J=1,NBE),(EDEP2(J),J=1,NBE)
        WRITE(9,'(/3X,''*** END ***'')')
        CLOSE(9)
      ENDIF
C
C  ------------------------  Print simulation results.
C
C     IEXIT: 1=transmitted, 2=backscattered, 3=absorbed.
C
      TOTN=SHN*WGHT0
      WRITE(26,3000)
 3000 FORMAT(///3X,35('*')/3X,'**   Program PENSLAB. Results.   **',
     1  /3X,35('*'))
C
      WRITE(26,3001) TSIM
 3001 FORMAT(/3X,'Simulation time ......................... ',
     1  1P,E13.6,' sec')
      TAVS=TOTN/TSIM
      WRITE(26,3002) TAVS
 3002 FORMAT(3X,'Simulation speed ........................ ',
     1  1P,E13.6,' showers/sec')
      WRITE(26,3003) TOTN
 3003 FORMAT(//3X,'Simulated primary particles ............. ',
     1  1P,E13.6)
C
      WRITE(26,3004) PRIM(1)
 3004 FORMAT(/3X,'Transmitted primary particles ........... ',
     1  1P,E13.6)
      WRITE(26,3005) PRIM(2)
 3005 FORMAT(3X,'Backscattered primary particles ......... ',
     1  1P,E13.6)
      WRITE(26,3006) PRIM(3)
 3006 FORMAT(3X,'Absorbed primary particles .............. ',
     1  1P,E13.6)
C
      FNT=1.0D0/TOTN
      FT=(PRIM(1)+SEC(KPARP,1))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(1)-PRIM(1)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,1)-SEC(KPARP,1)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(26,3007) FT,ERR
 3007 FORMAT(/3X,'Fractional transmission ............ ',
     1  1P,E13.6,' +-',E8.1)
      FB=(PRIM(2)+SEC(KPARP,2))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(2)-PRIM(2)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,2)-SEC(KPARP,2)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(26,3008) FB,ERR
 3008 FORMAT(3X,'Fractional backscattering .......... ',
     1  1P,E13.6,' +-',E8.1)
      FA=PRIM(3)*FNT
      ERR=3.0D0*FNT*SQRT(ABS(PRIM2(3)-PRIM(3)**2*FNT))
      WRITE(26,3009) FA,ERR
 3009 FORMAT(3X,'Fractional absorption .............. ',
     1  1P,E13.6,' +-',E8.1)
C
      DO K=1,3
        DO I=1,3
          WSEC2(K,I)=3.0D0*FNT*SQRT(ABS(SEC2(K,I)-SEC(K,I)**2*FNT))
          WSEC(K,I)=SEC(K,I)*FNT
        ENDDO
      ENDDO
      WRITE(26,3010)
     1  WSEC(1,1),WSEC(2,1),WSEC(3,1),WSEC2(1,1),WSEC2(2,1),WSEC2(3,1),
     1  WSEC(1,2),WSEC(2,2),WSEC(3,2),WSEC2(1,2),WSEC2(2,2),WSEC2(3,2),
     1  WSEC(1,3),WSEC(2,3),WSEC(3,3),WSEC2(1,3),WSEC2(2,3),WSEC2(3,3)
 3010 FORMAT(/3X,'Secondary-particle generation probabilities:',
     1 /19X,46('-'),
     1 /19X,'|  electrons   |   photons    |  positrons   |',1P,
     1 /3X,62('-')/3X,'|  transmitted  |',3(E13.6,1X,'|'),
     1 /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1 /3X,62('-')/3X,'| backscattered |',3(E13.6,1X,'|'),
     1 /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1 /3X,62('-')/3X,'|   absorbed    |',3(E13.6,1X,'|'),
     1 /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1 /3X,62('-'))
C
C  ************  Average values for primary particles.
C
      DO KC=1,8
        IF(AVI(KC).GT.1.0D-16) THEN
          WAVI2(KC)=3.0D0*FNT*SQRT(ABS(AVI2(KC)-AVI(KC)**2*FNT))
          WAVI(KC)=AVI(KC)*FNT
        ELSE
          WAVI2(KC)=0.0D0
          WAVI(KC)=0.0D0
        ENDIF
      ENDDO
      DO I=1,3
        DF=1.0D0/MAX(PRIM(I),1.0D0)
        WAVTL2(I)=3.0D0*DF*SQRT(ABS(AVTL2(I)-AVTL(I)**2*DF))
        WAVTL(I)=AVTL(I)*DF
        IF(I.LT.3) THEN
          WAVE2(I)=3.0D0*DF*SQRT(ABS(AVE2(I)-AVE(I)**2*DF))
          WAVE(I)=AVE(I)*DF
          WAVW2(I)=3.0D0*DF*SQRT(ABS(AVW2(I)-AVW(I)**2*DF))
          WAVW(I)=AVW(I)*DF
          WAVA2(I)=3.0D0*DF*RA2DE*SQRT(ABS(AVA2(I)-AVA(I)**2*DF))
          WAVA(I)=AVA(I)*RA2DE*DF
        ENDIF
      ENDDO
C
      WRITE(26,3011)
 3011 FORMAT(//3X,'Mean number of events per primary track:')
      IF(KPARP.EQ.2) THEN
        WRITE(26,3013) WAVI(1),WAVI2(1)
 3013   FORMAT(6X,'Coherent (Rayleigh) ............. ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3014) WAVI(2),WAVI2(2)
 3014   FORMAT(6X,'Incoherent (Compton) ............ ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3015) WAVI(3),WAVI2(3)
 3015   FORMAT(6X,'Photoelectric absorption ........ ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3016) WAVI(4),WAVI2(4)
 3016   FORMAT(6X,'Pair production ................. ',1P,E13.6,
     1    ' +-',E8.1)
      ELSE
        WRITE(26,3017) WAVI(1),WAVI2(1)
 3017   FORMAT(6X,'Hinges (soft events) ............ ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3018) WAVI(2),WAVI2(2)
 3018   FORMAT(6X,'Hard elastic collisions ......... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3019) WAVI(3),WAVI2(3)
 3019   FORMAT(6X,'Hard inelastic collisions ....... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3020) WAVI(4),WAVI2(4)
 3020   FORMAT(6X,'Hard bremsstrahlung emissions ... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3021) WAVI(5),WAVI2(5)
 3021   FORMAT(6X,'Inner-shell ionizations ......... ',1P,E13.6,
     1    ' +-',E8.1)
        WRITE(26,3022) WAVI(7),WAVI2(7)
 3022   FORMAT(6X,'Delta interactions .............. ',1P,E13.6,
     1    ' +-',E8.1)
      ENDIF
C
      WRITE(26,3025)
 3025 FORMAT(/3X,'Average final energy:')
      WRITE(26,3026) WAVE(1),WAVE2(1)
 3026 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1,' eV')
      WRITE(26,3027) WAVE(2),WAVE2(2)
 3027 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1,' eV')
C
      WRITE(26,3028)
 3028 FORMAT(/3X,'Average track length:')
      WRITE(26,3029) WAVTL(1),WAVTL2(1)
 3029 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1,' cm')
      WRITE(26,3030) WAVTL(2),WAVTL2(2)
 3030 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1,' cm')
      WRITE(26,3031) WAVTL(3),WAVTL2(3)
 3031 FORMAT(6X,'Absorbed particles .............. ',1P,E13.6,
     1  ' +-',E8.1,' cm')
C
      WRITE(26,3032)
 3032 FORMAT(/3X,'Mean value of the polar cosine of the exit',
     1  ' direction:')
      WRITE(26,3033) WAVW(1),WAVW2(1)
 3033 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1)
      WRITE(26,3034) WAVW(2),WAVW2(2)
 3034 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1)
C
      WRITE(26,3035)
 3035 FORMAT(/3X,'Mean value of the polar angle of the exit dir',
     1  'ection:')
      WRITE(26,3036) WAVA(1),WAVA2(1)
 3036 FORMAT(6X,'Transmitted particles ........... ',1P,E13.6,
     1  ' +-',E8.1,' deg')
      WRITE(26,3037) WAVA(2),WAVA2(2)
 3037 FORMAT(6X,'Backscattered particles ......... ',1P,E13.6,
     1  ' +-',E8.1,' deg')
C
C  ****  Deposited energy.
C
      DF=1.0D0/TOTN
      QER=3.0D0*DF*SQRT(ABS(TDEP2-TDEP**2*DF))
      QAV=TDEP*DF
      IF(QER.GT.1.0D-10*ABS(QAV)) THEN
        EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
      ELSE
        EFFIC=0.0D0
      ENDIF
      WRITE(26,3039) QAV,QER,EFFIC
 3039 FORMAT(/3X,'Average deposited energy ........... ',1P,E13.6,
     1  ' +-',E8.1,' eV',/41X,'(efficiency =',E9.2,')')
C
      WRITE(26,3040) ISEED1,ISEED2
 3040 FORMAT(/3X,'Last random seeds = ',I10,' , ',I10)
C
C  ------------------------  Print tallied distributions.
C
      IF(ISPEC.EQ.1) THEN
C
C  ************  Energy spectrum of the source.
C
        OPEN(9,FILE='psource.dat')
        WRITE(9,9000)
 9000 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Source energy spectrum.',
     1 /1X,'#  1st column: E (eV). 2nd column: spectrum (1/eV).',
     1 /1X,'#  3rd and 4th columns: simul. pdf limits (eV/particle).',/)
        PTOT=0.0D0
        WRITE(9,'(1X,1P,4E14.6)') ES(1),PTOT,PTOT,PTOT
        DO KE=1,NSEB
          PTOT=PTOT+PTS(KE)
        ENDDO
        DO KE=1,NSEB
          YAV=SHIST(KE)*DF
          YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV)*DF))
          EINTL=ES(KE+1)-ES(KE)
          IF(EINTL.GT.1.0D-15) THEN
            FACT=1.0D0/EINTL
          ELSE
            FACT=1.0D15
          ENDIF
          WRITE(9,'(1X,1P,4E14.6)') ES(KE),PTS(KE)*FACT/PTOT,
     1      (YAV-YERR)*FACT,(YAV+YERR)*FACT
          WRITE(9,'(1X,1P,4E14.6)') ES(KE+1),PTS(KE)*FACT/PTOT,
     1      (YAV-YERR)*FACT,(YAV+YERR)*FACT
        ENDDO
        PTOT=0.0D0
        WRITE(9,'(1X,1P,4E14.6)') ES(NSEB+1),PTOT,PTOT,PTOT
        CLOSE(9)
      ENDIF
C
C  ************  Depth-dose.
C
      OPEN(9,FILE='psddose.dat')
      WRITE(9,9010)
 9010 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Depth dose distribution.',
     1 /1X,'#  1st column: Z (cm). 4th column: corresponding layer.',
     1 /1X,'#  2nd and 3rd columns: dose and STU (eV/',
     1      '((g/cm**2)*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBZ
        XX=(K-0.5D0)*BSZ
        YERR=3.0D0*SQRT(ABS(DOSZ2(K)-DOSZ(K)**2*DF))
        YAV=DOSZ(K)*DF/BSZ
        YERR=YERR*DF/BSZ
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Depth distribution of deposited charge.
C
      OPEN(9,FILE='psdchar.dat')
      WRITE(9,9020)
 9020 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Depth distribution of deposited charge.',
     1 /1X,'#  1st column: Z (cm). 4th column: corresponding layer.',
     1 /1X,'#  2nd and 3rd columns:',
     1      ' charge density and STU (e/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBZ
        XX=(K-0.5D0)*BSZ
        YERR=3.0D0*SQRT(ABS(CHRZ2(K)-CHRZ(K)**2*DF))
        YAV=CHRZ(K)*DF/BSZ
        YERR=YERR*DF/BSZ
        WRITE(9,'(1X,1P,3E14.6)') XX,YAV,YERR
      ENDDO
      CLOSE(9)
C
C  ************  Track length distributions of primary particles.
C
C  ****  Transmitted particles.
      OPEN(9,FILE='pstltr.dat')
      WRITE(9,9030)
 9030 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Track length distribution of transmitted particles.',
     1 /1X,'#  1st column: track length (cm).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBTL
        XX=TLMIN+(K-0.5D0)*BSTL
        YERR=3.0D0*SQRT(ABS(PDTL2(1,K)-PDTL(1,K)**2*DF))
        YAV=PDTL(1,K)*DF/BSTL
        YERR=YERR*DF/BSTL
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered particles.
      OPEN(9,FILE='pstlbk.dat')
      WRITE(9,9040)
 9040 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Track length distribution of backscattered particles.',
     1 /1X,'#  1st column: track length (cm).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBTL
        XX=TLMIN+(K-0.5D0)*BSTL
        YERR=3.0D0*SQRT(ABS(PDTL2(2,K)-PDTL(2,K)**2*DF))
        YAV=PDTL(2,K)*DF/BSTL
        YERR=YERR*DF/BSTL
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Absorbed particles.
      OPEN(9,FILE='pstlab.dat')
      WRITE(9,9050)
 9050 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Track length distribution of absorbed particles.',
     1 /1X,'#  1st column: track length (cm).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(cm*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBTL
        XX=TLMIN+(K-0.5D0)*BSTL
        YERR=3.0D0*SQRT(ABS(PDTL2(3,K)-PDTL(3,K)**2*DF))
        YAV=PDTL(3,K)*DF/BSTL
        YERR=YERR*DF/BSTL
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Angular distributions of emerging particles.
C
      OPEN(9,FILE='psanel.dat')
      WRITE(9,9060)
 9060 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Angular distribution of emerging electrons.',
     1 /1X,'#  1st column: THETA (deg).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*TWOPI
        YERR=3.0D0*SQRT(ABS(PDA2(1,K)-PDA(1,K)**2*DF))
        YAV=PDA(1,K)*DF/DSANG
        YERR=YERR*DF/DSANG
        WRITE(9,'(1X,1P,6E14.6)')
     1       XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
      OPEN(9,FILE='psanga.dat')
      WRITE(9,9070)
 9070 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Angular distribution of emerging photons.',
     1 /1X,'#  1st column: THETA (deg).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*TWOPI
        YERR=3.0D0*SQRT(ABS(PDA2(2,K)-PDA(2,K)**2*DF))
        YAV=PDA(2,K)*DF/DSANG
        YERR=YERR*DF/DSANG
        WRITE(9,'(1X,1P,6E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
      OPEN(9,FILE='psanpo.dat')
      WRITE(9,9080)
 9080 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Angular distribution of emerging positrons.',
     1 /1X,'#  1st column: THETA (deg).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*TWOPI
        YERR=3.0D0*SQRT(ABS(PDA2(3,K)-PDA(3,K)**2*DF))
        YAV=PDA(3,K)*DF/DSANG
        YERR=YERR*DF/DSANG
        WRITE(9,'(1X,1P,6E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Energy distributions of emerging particles.
C
C  ****  Transmitted electrons.
      OPEN(9,FILE='pseneltr.dat')
      WRITE(9,9090)
 9090 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Energy distribution of transmitted electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(1)
        YERR=3.0D0*SQRT(ABS(PDE2(1,1,K)-PDE(1,1,K)**2*DF))
        YAV=PDE(1,1,K)*DF/BSE(1)
        YERR=YERR*DF/BSE(1)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered electrons.
      OPEN(9,FILE='psenelbk.dat')
      WRITE(9,9100)
 9100 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Energy distribution of backscattered electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(1)
        YERR=3.0D0*SQRT(ABS(PDE2(1,2,K)-PDE(1,2,K)**2*DF))
        YAV=PDE(1,2,K)*DF/BSE(1)
        YERR=YERR*DF/BSE(1)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Transmitted photons.
      OPEN(9,FILE='psengatr.dat')
      WRITE(9,9110)
 9110 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Energy distribution of transmitted photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(2)
        YERR=3.0D0*SQRT(ABS(PDE2(2,1,K)-PDE(2,1,K)**2*DF))
        YAV=PDE(2,1,K)*DF/BSE(2)
        YERR=YERR*DF/BSE(2)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered photons.
      OPEN(9,FILE='psengabk.dat')
      WRITE(9,9120)
 9120 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Energy distribution of backscattered photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(2)
        YERR=3.0D0*SQRT(ABS(PDE2(2,2,K)-PDE(2,2,K)**2*DF))
        YAV=PDE(2,2,K)*DF/BSE(2)
        YERR=YERR*DF/BSE(2)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Transmitted positrons.
      OPEN(9,FILE='psenpotr.dat')
      WRITE(9,9130)
 9130 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Energy distribution of transmitted positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(3)
        YERR=3.0D0*SQRT(ABS(PDE2(3,1,K)-PDE(3,1,K)**2*DF))
        YAV=PDE(3,1,K)*DF/BSE(3)
        YERR=YERR*DF/BSE(3)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Backscattered positrons.
      OPEN(9,FILE='psenpobk.dat')
      WRITE(9,9140)
 9140 FORMAT(
     1  1X,'#  Results from PENSLAB. ',
     1 /1X,'#  Energy distribution of backscattered positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DF=1.0D0/TOTN
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE(3)
        YERR=3.0D0*SQRT(ABS(PDE2(3,2,K)-PDE(3,2,K)**2*DF))
        YAV=PDE(3,2,K)*DF/BSE(3)
        YERR=YERR*DF/BSE(3)
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Distribution of deposited energy.
C
      OPEN(9,FILE='psedepm.dat')
      WRITE(9,9150)
 9150 FORMAT(
     1  1X,'#  Results from PENSLAB.',
     1 /1X,'#  Distribution of deposited energy.',
     1 /1X,'#  WARNING: May be strongly biased if interaction ',
     1  'forcing is used!',
     1 /1X,'#  1st column: deposited energy (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      IF(KPARP.EQ.3) WRITE(9,9151)
 9151 FORMAT(1X,'#  WARNING: For positrons, tallied only for E<E0.',/)
      DF=1.0D0/TOTN
      DO J=1,NBE
        XX=(J-0.5D0)*BSDE
        YERR=3.0D0*SQRT(ABS(EDEP2(J)-EDEP(J)**2*DF))
        YAV=EDEP(J)*DF/BSDE
        YERR=YERR*DF/BSDE
        WRITE(9,'(1X,1P,3E14.6)')
     1    XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ****  Continue if the DUMPTO option is active and IRETRN=1.
C
      IF(IRETRN.EQ.1) THEN
        WRITE(26,'(/3X,72(''-''))')
        WRITE(6,3045) SHN
 3045   FORMAT(2X,'Number of simulated showers =',1P,E14.7)
        GO TO 101
      ENDIF
C
      WRITE(26,3050)
 3050 FORMAT(//3X,'***  END  ***')
      CLOSE(26)
      STOP
      END
C  *********************************************************************
C                       SUBROUTINE GCONE
C  *********************************************************************
      SUBROUTINE GCONE(UF,VF,WF)
C
C  This subroutine samples a random direction uniformly within a cone
C  with central axis in the direction (THETA,PHI) and aperture ALPHA.
C  Parameters are initialised by calling subroutine GCONE0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      EXTERNAL RAND
C  ****  Define a direction relative to the z-axis.
      WT=CAPER+(1.0D0-CAPER)*RAND(1.0D0)
      DF=TWOPI*RAND(2.0D0)
      SUV=SQRT(1.0D0-WT*WT)
      UT=SUV*COS(DF)
      VT=SUV*SIN(DF)
C  **** Rotate to the beam axis direction.
      UF=CPCT*UT-SPHI*VT+CPST*WT
      VF=SPCT*UT+CPHI*VT+SPST*WT
      WF=-STHE*UT+CTHE*WT
C  ****  Ensure normalisation.
      DXY=UF*UF+VF*VF
      DXYZ=DXY+WF*WF
      IF(ABS(DXYZ-1.0D0).GT.1.0D-14) THEN
        FNORM=1.0D0/SQRT(DXYZ)
        UF=FNORM*UF
        VF=FNORM*VF
        WF=FNORM*WF
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GCONE0
C  *********************************************************************
      SUBROUTINE GCONE0(THETA,PHI,ALPHA)
C
C  This subroutine defines the parameters for sampling random directions
C  uniformly within a cone with axis in the direction (THETA,PHI) and
C  aperture ALPHA (in rad).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      CPCT=COS(PHI)*COS(THETA)
      CPST=COS(PHI)*SIN(THETA)
      SPCT=SIN(PHI)*COS(THETA)
      SPST=SIN(PHI)*SIN(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      STHE=SIN(THETA)
      CTHE=COS(THETA)
      CAPER=COS(ALPHA)
      RETURN
      END
