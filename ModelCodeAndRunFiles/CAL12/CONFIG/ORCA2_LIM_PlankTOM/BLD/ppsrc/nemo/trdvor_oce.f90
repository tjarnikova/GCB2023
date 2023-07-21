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

MODULE trdvor_oce
   !!======================================================================
   !!                   ***  MODULE trdvor_oce  ***
   !! Ocean trends :   set vorticity trend variables
   !!======================================================================
   !! History :  1.0  !  04-2006  (L. Brunier, A-M. Treguier) Original code 
   !!----------------------------------------------------------------------
   
   USE par_oce      ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   !                                               !!* vorticity trends index
   INTEGER, PUBLIC, PARAMETER ::   jpltot_vor = 11  !: Number of vorticity trend terms
   !
   INTEGER, PUBLIC, PARAMETER ::   jpvor_prg =  1   !: Pressure Gradient Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_keg =  2   !: KE Gradient Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_rvo =  3   !: Relative Vorticity Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_pvo =  4   !: Planetary Vorticity Term Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_ldf =  5   !: Horizontal Diffusion Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_zad =  6   !: Vertical Advection Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_zdf =  7   !: Vertical Diffusion Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_spg =  8   !: Surface Pressure Grad. Trend
   INTEGER, PUBLIC, PARAMETER ::   jpvor_bev =  9   !: Beta V
   INTEGER, PUBLIC, PARAMETER ::   jpvor_swf = 10   !: wind stress forcing term
   INTEGER, PUBLIC, PARAMETER ::   jpvor_bfr = 11   !: bottom friction term

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: trdvor_oce.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trdvor_oce
