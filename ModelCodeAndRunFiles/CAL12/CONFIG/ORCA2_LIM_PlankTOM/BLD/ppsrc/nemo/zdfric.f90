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

MODULE zdfric
   !!======================================================================
   !!                       ***  MODULE  zdfric  ***
   !! Ocean physics:  vertical mixing coefficient compute from the local
   !!                 Richardson number dependent formulation
   !!======================================================================
   !! History :  OPA  ! 1987-09  (P. Andrich)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!            7.0  ! 1996-01  (G. Madec)  complete rewriting of multitasking suppression of common work arrays
   !!            8.0  ! 1997-06  (G. Madec)  complete rewriting of zdfmix
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.3.1! 2011-09  (P. Oddo) Mixed layer depth parameterization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :              NO Richardson dependent vertical mixing
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .FALSE.   !: Richardson mixing flag
CONTAINS
   SUBROUTINE zdf_ric_init         ! Dummy routine
   END SUBROUTINE zdf_ric_init
   SUBROUTINE zdf_ric( kt )        ! Dummy routine
      WRITE(*,*) 'zdf_ric: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_ric

   !!======================================================================
END MODULE zdfric
