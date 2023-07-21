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

      SUBROUTINE bgcint(kt)
!
      USE trp_trc
      USE sms_planktom
      USE trc
      USE oce_trc
      USE dom_oce , ONLY :   gphit      =>   gphit      !: latitude  of t-point (degre)
      USE dom_oce , ONLY :   nmonth    =>   nmonth      !: Current month
      USE sbc_oce , ONLY :   fr_i      =>   fr_i        !: ice fraction (between 0 to 1)
      IMPLICIT NONE
      INTEGER kt
      INTEGER ji, jj, jk, jl
      INTEGER ipdtant, ipdtmo, nvit1t, nvit2t, lecvit
      REAL    zpdtan, zt, ztmax  
      REAL    zice1, zice2, zfacti, zfactl, zfactd, zmonth, ztemp
!
      lecvit = 12
      zt     = ( float ( kt) + zdemi) / zpdtmo
      
!  recherche de l'indice des enregistrements
!  du modele dynamique encadrant le pas de temps kt.
!  --------------------------------------------------
!
      xtvit = zt - int ( zt)
      nvit1t = int (( float ( kt) + zdemi)/ zpdtmo)
      nvit2t = nvit1t+1
      nvit1t = MOD ( nvit1t-1, lecvit)+1
      nvit2t = MOD ( nvit2t-1, lecvit)+1
!
!
          DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
              dust(ji,jj) = (1.-xtvit)*dustmo(ji,jj,nvit1t) &
     &            +xtvit*dustmo(ji,jj,nvit2t)
            END DO
          END DO
!
! For all pft's - including bacteria
!
       DO jl = jpbac, jpdia+jppft-1
         DO jk = 1, jpkm1
           DO jj = 2, nlcj-1
             DO ji = 2, nlci-1
               tgfunc(ji,jj,jk,jl) = exp(-1.*(tsn(ji,jj,jk,1)-rn_mutpft(jl))**2./rn_mudpft(jl)**2.)
             END DO
           END DO
         END DO
       END DO
!
! coccolithophorids growth is reduced below ztemp:
!
       ztemp = 10.
       DO jj = 2, nlcj-1
         DO ji = 2, nlci-1
             ztmax = (1./5. + 4./5.*max(tsn(ji,jj,jk,1),-1.8)/(ztemp+rtrn))*tgfunc(ji,jj,jk,jpcoc)
             tgfunc(ji,jj,jk,jpcoc) = min(ztmax,tgfunc(ji,jj,jk,jpcoc))
!
! increase the growth rate of mac when ice is between 0.1-0.3 to represent enhanced recruitment
! as in Wiedenmann et al., Limnology and Oceanography 2009
!
! zfacti is 1 when ice is greater then 0.1 and less than 0.3 and 0 otherwise 
!
             zice1  = fr_i(ji,jj)-0.1
             zice2  = fr_i(ji,jj)-0.3
             zfacti = max(0.,zice1)/(zice1+rtrn) * min(0.,zice2)/(zice2+rtrn)
!
! zfactl is 1 during Aug-Nov in the SH and Feb-May in the NH and 0 otherwise
!
             zmonth = float(nmonth)
             zfactl =  min(gphit(ji,jj),0.)/(gphit(ji,jj)+rtrn) *          & 
     &                    (max(0.,(zmonth-7.))/(zmonth-7.+rtrn)) *         &
     &                    (min(0.,(zmonth-12.))/(zmonth-12.+rtrn)) +       &
     &                 max(gphit(ji,jj),0.)/(gphit(ji,jj)+rtrn) *          &
     &                    (max(0.,(zmonth-1.))/(zmonth-1.+rtrn)) *         &
     &                    (min(0.,(zmonth-6.))/(zmonth-6.+rtrn))                    
!
! zfactd is 1 when the depth does not exceed 600 m and 0 otherwise
!
             zfactd = min(0.,mdept(ji,jj)-600.)/(mdept(ji,jj)-600.+rtrn)
!
! multiply growth by a factor rn_icemac (from Figure 4 in Wiedeman) when the above conditions are true
!
           DO jk = 1, jpkm1
             tgfunc(ji,jj,jk,jpmac) = max(1.,zfacti*zfactl*zfactd*rn_icemac)*    &
     &                                tgfunc(ji,jj,jk,jpmac)
!
           END DO
         END DO
       END DO
!
! assign susceptibility to advection (rn_trnmac) in coastal or ice-covered areas
!
       DO jj = 1, nlcj
         DO ji = 1, nlci
!
! zfacti is rn_trnmac when ice is greater then 0.1 and 1 otherwise 
!
           zice1  = fr_i(ji,jj)-0.1
           zfacti = 1. - max(0.,zice1)/(zice1+rtrn)*rn_trnmac 
!
! trnmac is rn_trnmac in shelf regions (<600 m), 1 otherwise
!
           trnmac(ji,jj) = 1. - min(0.,mdept(ji,jj)-600.)/(mdept(ji,jj)-600.+rtrn)*rn_trnmac
!
! combine ice and shelf effect
!
           trnmac(ji,jj) = min(trnmac(ji,jj),zfacti)
!
         END DO
       END DO
!
      RETURN
      END
