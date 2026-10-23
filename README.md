# Generating toponium events with MG5_aMC

This guide documents the MG5_aMC modifications for colour-singlet, S-wave top-antitop production with Green's-function re-weighting. It follows section 3 of [2411.18962](https://arxiv.org/abs/2411.18962), updated for the single-grid format supplied in this directory and for a consistent choice of the top mass and width.

The prescribed event-selection region for this implementation is **`-6 <= E <= 4 GeV` and `p* <= 50 GeV`**.

**The modifications described here have been tested with MG5_aMC version 3.5.7.** As the generated Fortran sources can differ between MG5_aMC versions, this guide provides the changes to apply to your own generated process.

## 1. Initialisation

For the di-leptonic case, generate the process with your MG5_aMC installation according to the instructions given in [2411.18962](https://arxiv.org/abs/2411.18962), then edit its generated `driver.f` and matrix-element source file `matrix1_orig.f` as described below. Copy and unpack the Green's function grid [Gmap_S_1725_131_coulomb.tar.gz](Gmap_S_1725_131_coulomb.tar.gz) (complex S-wave Green's function at `mt = 172.5 GeV`, `Gamma_t = 1.315 GeV`) into the directory `SubProcesses/P1_gg_wpbwmbx_wp_lvl_wm_lvl`. For the single leptonic case, the procedure is similar but two occurrences of the files are concerned, located in the directories `SubProcesses/P1_gg_wpbwmbx_wp_lvl_wm_qq` and `SubProcesses/P2_gg_wpbwmbx_wp_qqwm_lvl`. More information can be found in [2512.03220](https://arxiv.org/abs/2512.03220).


## 2. Fortran modifications

### 2.1 Load the grid in `driver.f`

In the generated `driver.f`, add this module before `PROGRAM DRIVER`:

```fortran
      MODULE GREENDATA
      IMPLICIT NONE
      INTEGER NMOM, NENE
      PARAMETER (NENE=1201,NMOM=801)
      REAL*8 EDATA(NENE,NMOM),PDATA(NENE,NMOM)
      REAL*8 RGDATA(NENE,NMOM),IGDATA(NENE,NMOM)
      SAVE
      END MODULE GREENDATA
```

Inside `PROGRAM DRIVER`, add `USE GREENDATA` before `IMPLICIT NONE`, then declare `INTEGER IGRN,JGRN` with the other declarations. Insert the following loading block after all declarations, before integration begins:

```fortran
      open(unit=99, file="../Gmap_S_1725_131_coulomb.txt")
      READ(99,*)
      READ(99,*)
      read(99,*) ((EDATA(IGRN,JGRN), PDATA(IGRN,JGRN),
     $   RGDATA(IGRN,JGRN), IGDATA(IGRN,JGRN), JGRN=1,NMOM,1),
     $   IGRN=1,NENE,1)
      close(99)
```

This minimal reader skips exactly two records and trusts all following data. 


### 2.2 Project the colour matrix in `MATRIX1`

In `matrix1_orig.f`, locate the function `MATRIX1` and add `USE GREENDATA` before its `IMPLICIT NONE` statement:

```fortran
      REAL*8 FUNCTION MATRIX1(P,NHEL,IC,IHEL)
      USE GREENDATA
      IMPLICIT NONE
```

Within the same function, locate the `COLOR DATA` block. The generated colour matrix includes both singlet and octet contributions. Taking the example of the di-lepton case, the original block reads

```fortran
      DATA (CF(I,1),I=1,2) /
     $ 5.333333333333333D+00,-6.666666666666666D-01/
C     1 T(1,2,5,8)
      DATA (CF(I,2),I=1,2) /
     $ -6.666666666666666D-01,5.333333333333333D+00/
C     1 T(2,1,5,8)
```

We should replace this block with the following to retain only the colour-singlet contribution:

```fortran
      DATA (CF(I,1),I=1,2) /
     $ 6.666666666666666D-01,6.666666666666666D-01/
C     1 T(1,2,5,8)
      DATA (CF(I,2),I=1,2) /
     $ 6.666666666666666D-01,6.666666666666666D-01/
C     1 T(2,1,5,8)
```
In the single-leptonic case, these numbers are all multiplied by 3.


### 2.3 Reconstruct the kinematics

Add the following declarations inside `MATRIX1`, before the executable code:

```fortran
      REAL*8 PMOM(0:3,3),GXOM,GYOM,GZOM,GEOM
      REAL*8 GMOM,GTHE,GPHI,GBET,GGAM,GMAS
      REAL*8 EXXX,MXXX
      COMPLEX*16 GREEN
      REAL*8 EMIN,EMAX,PMAX
      PARAMETER (EMIN=-6.D0,EMAX=4.D0,PMAX=5.D1)
```

Before the colour contraction starting with `MATRIX1 = 0.D0`, insert the
following reconstruction and rest-frame transformation from the paper:

```fortran
      DO M = 0, 3
         PMOM(M, 2) = 0
         PMOM(M, 3) = 0
         DO I = 3, 5
            PMOM(M, 2) = PMOM(M, 2) + P(M, I)
            PMOM(M, 3) = PMOM(M, 3) + P(M, I+3)
         ENDDO
         PMOM(M, 1) = PMOM(M, 2) + PMOM(M, 3)
      ENDDO
C     
C     Lorentz transformation factors
      GMOM = SQRT(PMOM(1, 1)**2 + PMOM(2, 1)**2 + PMOM(3, 1)**2)
      GMAS = SQRT(PMOM(0, 1)**2 - GMOM**2)
      GTHE = DACOS(PMOM(3, 1)/(GMOM+1.e-8))
      GPHI = ATAN2(PMOM(2, 1), PMOM(1, 1))
      GBET = GMOM/PMOM(0, 1)
      GGAM = PMOM(0, 1)/GMAS
C     
C     Rotating around Z-axis
      DO M = 1, 3
         GXOM = PMOM(1, M)
         GYOM = PMOM(2, M)
         PMOM(1, M) = COS(GPHI)*GXOM + SIN(GPHI)*GYOM
         PMOM(2, M) = COS(GPHI)*GYOM - SIN(GPHI)*GXOM
      ENDDO
C     
C     Rotating around Y-axis
      DO M = 1, 3
         GXOM = PMOM(1, M)
         GZOM = PMOM(3, M)
         PMOM(1, M) = COS(GTHE)*GXOM - SIN(GTHE)*GZOM
         PMOM(3, M) = COS(GTHE)*GZOM + SIN(GTHE)*GXOM
      ENDDO
C     
C     Boost along Z-axis
      DO M = 1, 3
         GEOM = PMOM(0, M)
         GZOM = PMOM(3, M)
         PMOM(0, M) = GGAM*GEOM - GGAM*GBET*GZOM
         PMOM(3, M) = GGAM*GZOM - GGAM*GBET*GEOM
      ENDDO
```

Here `PMOM(:,1)` is the pair, `PMOM(:,2)` the top and `PMOM(:,3)` the antitop. 

### 2.4 Interpolate G and construct G/G0

Append the following new routine after the existing routines in the matrix source. It imports the grid from `GREENDATA` and takes the model mass and width as arguments. 

```fortran
      SUBROUTINE CALCGREEN(EOMX,MOMX,MT,GT,GREEN)
      USE GREENDATA
      IMPLICIT NONE

C     Arguments
      REAL*8 EOMX, MOMX, MT, GT
      COMPLEX*16 GREEN

C     Variables
      REAL*8 DENE, DMOM, AA, BB
      INTEGER II, JJ
      COMPLEX*16 INVG0, GINT

C     Outside the tabulated region
      GREEN = DCMPLX(1D-16,0D0)
      IF (EOMX.LT.EDATA(1,1) .OR. EOMX.GT.EDATA(NENE,1)) RETURN
      IF (MOMX.LT.PDATA(1,1) .OR. MOMX.GT.PDATA(1,NMOM)) RETURN

C     Uniform grid: infer spacings from file
      DENE = EDATA(2,1)-EDATA(1,1)
      DMOM = PDATA(1,2)-PDATA(1,1)
      II = MAX(1,MIN(INT((EOMX-EDATA(1,1))/DENE)+1,NENE-1))
      JJ = MAX(1,MIN(INT((MOMX-PDATA(1,1))/DMOM)+1,NMOM-1))


C     Bilinear interpolation of the complex Green function
      AA = (EOMX-EDATA(II,1))/DENE
      BB = (MOMX-PDATA(1,JJ))/DMOM
      GINT =
     $ (1D0-AA)*(1D0-BB)*DCMPLX(RGDATA(II,JJ),IGDATA(II,JJ))
     $ + AA*(1D0-BB)*DCMPLX(RGDATA(II+1,JJ),IGDATA(II+1,JJ))
     $ + (1D0-AA)*BB*DCMPLX(RGDATA(II,JJ+1),IGDATA(II,JJ+1))
     $ + AA*BB*DCMPLX(RGDATA(II+1,JJ+1),IGDATA(II+1,JJ+1))

C     Ratio of Green's functions
      INVG0 = DCMPLX(EOMX-MOMX**2/MT,GT)
      GREEN = GINT*INVG0
      END
```

### 2.5 Apply the weight and define the selected region

After the colour contraction has filled `MATRIX1`, and before the
`SDE_STRAT` bookkeeping, insert:

```fortran
      EXXX = PMOM(0,1)-2D0*MDL_MT
      MXXX = SQRT(PMOM(1,2)**2+PMOM(2,2)**2+PMOM(3,2)**2)
      IF (EXXX.LT.EMIN .OR. EXXX.GT.EMAX .OR. MXXX.GT.PMAX) THEN
          MATRIX1 = MATRIX1*1.D-16
      ELSE
         CALL CALCGREEN(EXXX,MXXX,MDL_MT,MDL_WT, GREEN)
         MATRIX1 = DBLE(GREEN*DCONJG(GREEN))*MATRIX1
      ENDIF
```

`MDL_MT` and `MDL_WT` are the model parameters exposed by the generated
`coupl.inc`.


## 3. Build and launch

The module must be compiled before any matrix source using it. From the generated process directory, run:

```sh
./bin/generate_events
```

Select hard-scattering generation first, with showering disabled. Confirm that
the edited routines are compiled and the grid is loaded in each worker.



## 4. LHE preparation before parton showering

Changing `CF` changes the matrix-element weight, but it does not by itself establish the event record's singlet colour flow and the LHE events must therefore be post-processed (specifically before invoking Pythia). This is achieved by means of the [`reprocess.py`](reprocess.py) script.
