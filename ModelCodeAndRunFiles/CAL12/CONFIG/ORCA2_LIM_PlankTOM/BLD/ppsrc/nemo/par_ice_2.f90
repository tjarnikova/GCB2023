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

MODULE par_ice_2
   !!======================================================================
   !!                       ***  MODULE par_ice_2   ***
   !! Sea-Ice model : definition of the parameters
   !!======================================================================
   !!----------------------------------------------------------------------
   !!  'key_lim2'                                       LIM-2 sea-ice model
   !!----------------------------------------------------------------------
   USE par_oce

   IMPLICIT NONE
   PUBLIC               ! allows par_oce and par_kind to be known in ice modules

   INTEGER, PUBLIC, PARAMETER ::   jpl        = 1              !: number of ice categories (only 1 in LIM-2)

   INTEGER, PUBLIC, PARAMETER ::   jplayers   = 2              !: number of vertical ice layers
   INTEGER, PUBLIC, PARAMETER ::   jplayersp1 = jplayers + 1   !: ???

   !!----------------------------------------------------------------------
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: par_ice_2.F90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE par_ice_2
