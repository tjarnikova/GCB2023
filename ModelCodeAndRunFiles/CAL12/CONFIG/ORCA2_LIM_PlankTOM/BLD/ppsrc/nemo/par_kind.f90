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

MODULE par_kind
   !!======================================================================
   !!                   ***  MODULE par_kind  ***
   !! Ocean :  define the kind of real for the whole model
   !!======================================================================
   !! History :  1.0  ! 2002-06  (G. Madec)  Original code
   !!            3.3  ! 2010-12  (G. Madec)  add a standard length of character strings
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE

   INTEGER, PUBLIC, PARAMETER ::   jpbyt   = 8    !: real size for mpp communications
   INTEGER, PUBLIC, PARAMETER ::   jpbytda = 4    !: real size in input data files 4 or 8

   ! Number model from which the SELECTED_*_KIND are requested:
   !             4 byte REAL       8 byte REAL
   ! CRAY:           -            precision = 13
   !                              exponent = 2465
   ! IEEE:      precision = 6     precision = 15
   !            exponent = 37     exponent = 307

   !                                                                !!** Floating point **
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision

   !                                                                !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   i4 = SELECTED_INT_KIND( 9)        !: single precision (integer 4)
   INTEGER, PUBLIC, PARAMETER ::   i8 = SELECTED_INT_KIND(14)        !: double precision (integer 8)
   
   !                                                                !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   lc = 256                          !: Lenght of Character strings

   !!----------------------------------------------------------------------
   !! NEMO 3.3 , NEMO Consortium (2010)
   !! $Id: par_kind.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
END MODULE par_kind
