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

MODULE c1d
   !!======================================================================
   !!                     ***  MODULE  c1d  ***
   !! Ocean domain  :  1D configuration
   !!=====================================================================
   !! History :   2.0  !  2004-09 (C. Ethe)     Original code
   !!             3.0  !  2008-04 (G. Madec)    adaptation to SBC
   !!             3.5  !  2013-10 (D. Calvert)  add namelist
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                           No use of 1D configuration
   !!----------------------------------------------------------------------
   USE par_kind         ! kind parameters

   LOGICAL, PUBLIC, PARAMETER ::   lk_c1d = .FALSE.   !: 1D config. flag de-activated
   REAL(wp)                   ::   rn_lat1d, rn_lon1d
   LOGICAL , PUBLIC           ::   ln_c1d_locpt = .FALSE. 

CONTAINS

   SUBROUTINE c1d_init               ! Dummy routine
   END SUBROUTINE c1d_init


   !!======================================================================
END MODULE c1d
