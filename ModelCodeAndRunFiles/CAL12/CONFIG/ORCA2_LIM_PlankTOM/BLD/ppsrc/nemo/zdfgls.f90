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

MODULE zdfgls
   !!======================================================================
   !!                       ***  MODULE  zdfgls  ***
   !! Ocean physics:  vertical mixing coefficient computed from the gls 
   !!                 turbulent closure parameterization
   !!======================================================================
   !! History :   3.0  !  2009-09  (G. Reffray)  Original code
   !!             3.3  !  2010-10  (C. Bricaud)  Add in the reference
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO TKE scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfgls = .FALSE.   !: TKE flag
CONTAINS
   SUBROUTINE zdf_gls_init           ! Empty routine
      WRITE(*,*) 'zdf_gls_init: You should not have seen this print! error?'
   END SUBROUTINE zdf_gls_init
   SUBROUTINE zdf_gls( kt )          ! Empty routine
      WRITE(*,*) 'zdf_gls: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_gls
   SUBROUTINE gls_rst( kt, cdrw )          ! Empty routine
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      WRITE(*,*) 'gls_rst: You should not have seen this print! error?', kt, cdrw
   END SUBROUTINE gls_rst

   !!======================================================================
END MODULE zdfgls

