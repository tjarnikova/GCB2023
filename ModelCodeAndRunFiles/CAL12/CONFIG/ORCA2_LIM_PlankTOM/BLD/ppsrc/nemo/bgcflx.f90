/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */


/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */

/* We do support the IEC 559 math functionality, real and complex.  */

/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */

/* We do not support C11 <threads.h>.  */

      SUBROUTINE bgcflx
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE bgcflx
!!!                     ******************
!!!
!!!
!!     PURPOSE.
!!     --------
!!          calculates gas exchange and chemistry at sea surface 
!!
!!     METHOD.
!!     -------
!!          solving system of two non-linear simultaneous equations for [H2CO3]i
!!          and [H+] 
!!
!!     EXTERNALS.
!!     ----------
!!          none.
!!
!!     REFERENCE.
!!     ----------
!!
!!          Broecker, W.S., and T.-H. Pend (1982) Tracers in the sea.
!!          Eldigio press, Lamont-Doherty Geological Observatory,
!!          Palisades, N.Y., 690 PP..
!!
!!          Mook, W.G. (1986), 13C in atmospheric CO2.
!!          Netherlands Journal of Sea Research, 20(2/3): 211-223.
!!
!!          Scarborough, J. (1958) Numerical Mathematical Analysis.
!!          Oxford Univertiry Press, London, 4TH ED., 576 PP..
!!
!!      VARIABLE           TYPE    PURPOSE.
!!      --------           ----    --------
!!
!!      alka              REAL    actual [ALK] [EQV/L]
!!      akb               REAL    1. dissoc. constant of boric acid
!!      krorr          INTEGER    counts iterations in Newton-Raphson algorithm
!!      ct1               REAL    actual valua of inorg. C
!!      x1                REAL    [H+] expressed as sqrt(AK1*AK2)/[H+]
!!      zpa               REAL    atmospheric pco2 [ppm]
!!      fld               REAL    [CO2] flux atmosphere -> ocean in [mol/m2/s]
!!      flu               REAL    [CO2] flux ocean -> atmosphere in [mol/m2/s]
!!      relw13            REAL    ratio [sum((13C)O2)]/[sum((12C)O2)] in surface layer
!!      rela13            REAL    ratio [sum((13C)O2)]/[sum((12C)O2)] in the atmosphere
!!
!!   MODIFICATIONS:
!!   --------------
!!      original      : 1988-07 E. Maier-Reimer
!!      additions     : 1998    O. Aumont
!!      modifications : 1999    C. Le Quere
!!      modifications : 2004    C. Le Quere
!!     -----------------------------------------------------------------
!!  parameters and commons
!! ======================
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE iom
      USE wrk_nemo
      IMPLICIT NONE
!!----------------------------------------------------------------------
!! local declarations
!! ==================
!
!
      INTEGER ji, jj, krorr,yearrc
      CHARACTER (len=20) :: diaadd
      REAL ttc,ttc2,ttc3,ttc4, ws
      REAL ttcdms
      REAL za, zpa, zvapor, ztkel
      REAL oxy16
      REAL zph,ah2,bot
      REAL ct1,alka
      REAL schmico2, schmio2
      REAL schmidms
      REAL piza,piph,piah2,piic
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
      dpco2=0.0
      pco2=0.0
      cflx=0.0
      flu16=0.0
      fld = 0.0
      flu=0.0

!
! surface chemistry (pCO2 AND [H+] in surface layer)
! --------------------------------------------------
      yearrc=ndastp/10000-int(yrcfc(1))+1 !1 record per year in atmco2cfc.dat
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
! total borate
! ------------
!
            bot  = borat(ji,jj,1)
            ct1  = trn(ji,jj,1,jpdic)
            alka = trn(ji,jj,1,jptal)
            zph  = amax1(hi(ji,jj,1),1.E-10)
            piic = trn(ji,jj,1,jppiic)
            piph = amax1(pihi(ji,jj,1),1.E-10)
!
! calculate [alk]([CO3--], [HCO3-])
! ---------------------------------
!
          DO krorr = 1, 15
            za = alka- &
     &          (akw3(ji,jj,1)/zph-zph+bot/(1.+zph/akb3(ji,jj,1)))
!
! calculate [H+] AND [H2CO3]
! --------------------------
!
            ah2 = sqrt((ct1-za)**2+4*(za*ak23(ji,jj,1)/ak13(ji,jj,1)) &
     &           * (2*ct1-za)) 
            zph = 0.5*ak13(ji,jj,1)/za*((ct1-za)+ah2)
            piza = alka- &
     &          (akw3(ji,jj,1)/piph-piph+bot/(1.+piph/akb3(ji,jj,1)))
            piah2 = sqrt((piic-piza)**2+4*(piza*ak23(ji,jj,1)/ak13(ji,jj,1)) &
     &           * (2*piic-piza))
            piph = 0.5*ak13(ji,jj,1)/piza*((piic-piza)+piah2)
          END DO
            hi(ji,jj,1)  = zph
            h2co3(ji,jj) = (2*ct1-za)/(2.+ak13(ji,jj,1)/zph)
            pihi(ji,jj,1)  = piph
            pih2co3(ji,jj) = (2*piic-piza)/(2.+ak13(ji,jj,1)/piph)
        END DO
      END DO
!
! compute fluxes
! --------------
!
! first compute gas exchange coefficients
! ---------------------------------------
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
          ttc      = min(39.,tsn(ji,jj,1,1))
          ttc2     = ttc**2
          ttc3     = ttc**3
          ttc4     = ttc**4
!
! use wind speed in m/s
! ---------------------
!
          ws = wndm(ji,jj)
!
! this is Wanninkhof (1992) equation 8 (with chemical enhancement), in cm/h
! -------------------------------------------------------------------------
!
          kgwanin(ji,jj) = (0.26*ws*ws + 2.5*(0.5246+ttc*(0.016256+ &
     &        ttc*0.00049946)))
!
! convert from cm/h to m/s and apply ice cover
! --------------------------------------------
!
          kgwanin(ji,jj) = kgwanin(ji,jj) /100./3600.  &
     &                   * (1-fr_i(ji,jj))
!
! compute Schmitt number for CO2 (Wanninkhof 1992)
! ------------------------------------------------
!
          schmico2 = 2073.1-125.62*ttc+3.6276*ttc2-0.043126*ttc3
          zkgco2(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmico2)
!
! this is Wanninkhof (1992) equation 3 (steady winds)
!
          kgwanin(ji,jj) = 0.27*ws*ws
          kgwanin(ji,jj) = kgwanin(ji,jj) /100./3600.  &
     &                   * (1-fr_i(ji,jj))
!
! compute Schmitt number for O2 (Wanninkhof 1992)
! -----------------------------------------------
!
          schmio2 = 1953.4-128.0*ttc+3.9918*ttc2-0.050091*ttc3
          zkgo2(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmio2)
! correct atmospheric CO2 (pco2at in ppm) for 100% water vapor 
! following Sarmiento et al., JGR 1992
! ------------------------------------------------------------
!
           ztkel   = tsn(ji,jj,1,1)+temzer
           zvapor = exp( 20.1050 - 0.0097982 * ztkel - 6163.10/ztkel)
           zpa = pco2at * (1. - zvapor)
!
! compute CO2 flux for the sea and air, 
! fld and flu are in mol/m2/s
! ---------------------------------------------
!
           fld(ji,jj) = zpa*chemc(ji,jj,1)*1.e3*zkgco2(ji,jj)*tmask(ji,jj,1)
           flu(ji,jj) = h2co3(ji,jj)      *1.e3*zkgco2(ji,jj)*tmask(ji,jj,1)
           dpco2(ji,jj)=min(h2co3(ji,jj)/chemc(ji,jj,1)-zpa,4000.)
!
! add tendency in mol/L/s
! -----------------------
!
           tra(ji,jj,1,jpdic)= tra(ji,jj,1,jpdic)+(fld(ji,jj)-flu(ji,jj))/1000./e3t_0(ji,jj,1)
           cflx(ji,jj) = (278.*(1.-zvapor)*chemc(ji,jj,1)               &
             -pih2co3(ji,jj))*zkgco2(ji,jj)*tmask(ji,jj,1)
           tra(ji,jj,1,jppiic)= tra(ji,jj,1,jppiic)+cflx(ji,jj)/e3t_0(ji,jj,1)
           cflx(ji,jj) = cflx(ji,jj)*1.e3
!
! Calculate cumulative flux mol/timestep
!
           qcumul(jpdic)=qcumul(jpdic)+(fld(ji,jj)-flu(ji,jj))*rfact*e1t(ji,jj) &
     &         *e2t(ji,jj)
!
! Compute O2 flux 
! ---------------
!
          oxy16 = trn(ji,jj,1,jpoxy)
          flu16(ji,jj) = (atcox*chemc(ji,jj,3)*(1.-zvapor) - oxy16) &
     &      * zkgo2(ji,jj) * tmask(ji,jj,1)
          tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy)+flu16(ji,jj)/e3t_0(ji,jj,1)
          qcumul(jpoxy)=qcumul(jpoxy)+flu16(ji,jj)*rfact*e1t(ji,jj)*e2t(ji,jj)
         END DO
       END DO
!      IF(lwp) WRITE(numout,*) 'fluch4 ',fluch4(20,10),pch4at,chemc(20,10,5),trn(20,10,1,jpch4),zkgch4(20,10)
!      IF(lwp) WRITE(numout,*) 'fluch4 ',pch4at,chemc(20,10,5),trn(20,10,1,jpch4)
!      CALL FLUSH(numout)
! Save diagnostics
! ---------------- 
 
 
!
 
       CALL iom_put("PICflx", cflx )
       cflx(2:nlci-1,2:nlcj-1) = (fld(2:nlci-1,2:nlcj-1)-flu(2:nlci-1,2:nlcj-1))
       pco2(2:nlci-1,2:nlcj-1) = min(h2co3(2:nlci-1,2:nlcj-1)/chemc(2:nlci-1,2:nlcj-1,1),4000.)
       flu16(2:nlci-1,2:nlcj-1) =  flu16(2:nlci-1,2:nlcj-1)*1000.
!       where (tmask(:,:,1) .eq. 0 )
!         cflx = ncf_fill
!         dpco2 = ncf_fill
!         flu16 = ncf_fill
!       end where
       diaadd="Cflx"
       CALL iom_put(diaadd, cflx(:,:) )
       diaadd="dpCO2"
       CALL iom_put(diaadd, dpco2(:,:) )
       CALL iom_put("pCO2", pco2(:,:) )
       diaadd="Oflx"
       CALL iom_put(diaadd, flu16(:,:) )
      RETURN
      END

