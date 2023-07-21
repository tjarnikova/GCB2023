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

       SUBROUTINE bgcatm

!
!-----------------------------------------------------------------------
!
!      Calculates the atmospheric CO2 concentration
!
!-----------------------------------------------------------------------
!
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      IMPLICIT NONE
!
      INTEGER iyy,is
      REAL zyr,slope,inter
!
      iyy = ndastp/10000
      zyr = float(iyy) + (float(nday_year)-0.5)/365.
!
      IF (zyr .lt. yrco2(1) .OR. zyr .gt. yrco2(nmaxrec)) THEN
        write(numout,*)'Caution: The date is outside tabulated values.'
        write(numout,*)'zyr, zyr(min) zyr(max) = ', zyr, yrco2(1),yrco2(nmaxrec)
        write(numout,*)'approximate pCO2'
        pco2at = 278+0.00009*max(zyr-1870.,0.)**2.847
      ENDIF
      DO is = 1, nmaxrec-1
        IF(zyr.ge.yrco2(is) .AND. zyr.le.yrco2(is+1)) THEN
          slope = (sipco2(is+1) - sipco2(is)) /(yrco2(is+1) - yrco2(is))
          inter = sipco2(is)-slope*yrco2(is)
          pco2at = slope*zyr + inter
          EXIT
        ENDIF
      END DO
      IF(lwp) WRITE(numout,*) 'Atmospheric xCO2 xN2O xCH4 concentration: ',&
     &  ndastp,zyr,pco2at,pn2oat,pch4at
      RETURN
      END
