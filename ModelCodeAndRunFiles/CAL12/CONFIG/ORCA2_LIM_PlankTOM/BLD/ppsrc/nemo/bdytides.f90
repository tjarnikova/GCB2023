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

MODULE bdytides
   !!======================================================================
   !!                       ***  MODULE  bdytides  ***
   !! Ocean dynamics:   Tidal forcing at open boundaries
   !!======================================================================
   !! History :  2.0  !  2007-01  (D.Storkey)  Original code
   !!            2.3  !  2008-01  (J.Holt)  Add date correction. Origins POLCOMS v6.3 2007
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D.Storkey and E.O'Dea)  bug fixes
   !!            3.4  !  2012-09  (G. Reffray and J. Chanut) New inputs + mods
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module         NO Unstruct Open Boundary Conditions for tides
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdytide_init             ! Empty routine
      WRITE(*,*) 'bdytide_init: You should not have seen this print! error?'
   END SUBROUTINE bdytide_init
   SUBROUTINE bdytide_update( kt, jit )   ! Empty routine
      WRITE(*,*) 'bdytide_update: You should not have seen this print! error?', kt, jit
   END SUBROUTINE bdytide_update
   SUBROUTINE bdy_dta_tides( kt, kit, time_offset )     ! Empty routine
      INTEGER, INTENT( in )            ::   kt          ! Dummy argument empty routine      
      INTEGER, INTENT( in ),OPTIONAL   ::   kit         ! Dummy argument empty routine
      INTEGER, INTENT( in ),OPTIONAL   ::   time_offset ! Dummy argument empty routine
      WRITE(*,*) 'bdy_dta_tides: You should not have seen this print! error?', kt, jit
   END SUBROUTINE bdy_dta_tides

   !!======================================================================
END MODULE bdytides

