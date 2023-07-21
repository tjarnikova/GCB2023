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

MODULE bdyice_lim
   !!======================================================================
   !!                       ***  MODULE  bdyice_lim  ***
   !! Unstructured Open Boundary Cond. :  Open boundary conditions for sea-ice (LIM2 and LIM3)
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!             3.4  !  2011    (D. Storkey)  rewrite in preparation for OBC-BDY merge
   !!              -   !  2012-01 (C. Rousset)  add lim3 and remove useless jk loop 
   !!----------------------------------------------------------------------
   !!---------------------------------------------------------------------------------
   !!   Default option                                                    Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_ice_lim( kt )      ! Empty routine
      WRITE(*,*) 'bdy_ice_lim: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_ice_lim

   !!=================================================================================
END MODULE bdyice_lim
