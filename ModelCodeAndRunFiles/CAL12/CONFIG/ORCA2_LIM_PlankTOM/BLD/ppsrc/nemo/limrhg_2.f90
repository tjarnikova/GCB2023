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

MODULE limrhg_2
   !!======================================================================
   !!                     ***  MODULE  limrhg_2  ***
   !!   Ice rheology :  performs sea ice rheology
   !!======================================================================
   !! History :  0.0  !  1993-12  (M.A. Morales Maqueda.)  Original code
   !!            1.0  !  1994-12  (H. Goosse) 
   !!            2.0  !  2003-12  (C. Ethe, G. Madec)  F90, mpp
   !!             -   !  2006-08  (G. Madec)  surface module, ice-stress at I-point
   !!             -   !  2009-09  (G. Madec)  Huge verctor optimisation
   !!            3.3  !  2009-05  (G.Garric, C. Bricaud) addition of the lim2_evp case
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option        Dummy module      NO VP & LIM-2 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_rhg_2( k1 , k2 )       ! Dummy routine
      WRITE(*,*) 'lim_rhg_2: You should not have seen this print! error?', k1, k2
   END SUBROUTINE lim_rhg_2

   !!==============================================================================
END MODULE limrhg_2
