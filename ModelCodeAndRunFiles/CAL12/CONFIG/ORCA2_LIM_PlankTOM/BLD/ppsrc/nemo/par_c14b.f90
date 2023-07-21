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

MODULE par_c14b
   !!======================================================================
   !!                        ***  par_c14b ***
   !! TOP :   set the C14 bomb parameters
   !!======================================================================
   !! History :   2.0  !  2008-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   USE par_planktom , ONLY : jp_planktom       !: number of tracers in PISCES
   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

   USE par_cfc    , ONLY : jp_cfc          !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_2d       !: number of 2D diag in CFC
   USE par_cfc    , ONLY : jp_cfc_3d       !: number of 3D diag in CFC
   USE par_cfc    , ONLY : jp_cfc_trd      !: number of biological diag in CFC


   IMPLICIT NONE

   INTEGER, PARAMETER ::   jp_lb      =  jp_planktom !    + jp_cfc     !: cum. number of pass. tracers
   INTEGER, PARAMETER ::   jp_lb_2d   =  jp_pisces_2d  + jp_cfc_2d  !:
   INTEGER, PARAMETER ::   jp_lb_3d   =  jp_pisces_3d  + jp_cfc_3d  !:
   INTEGER, PARAMETER ::   jp_lb_trd  =  jp_pisces_trd + jp_cfc_trd !:
   
   !!---------------------------------------------------------------------
   !!   'key_c14b'   :                                   C14 bomb tracer
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_c14b     = .TRUE.      !: C14 bomb flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b     =  3          !: number of passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b_2d  =  2          !: additional 2d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b_3d  =  1          !: additional 3d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b_trd =  0          !: number of sms trends for C14
   INTEGER, PUBLIC, PARAMETER ::   jpc14       = jp_lb + 2   !: assign an index in trc arrays for C14 bomb 
   INTEGER, PUBLIC, PARAMETER ::   jpb14       = jp_lb + 1
   INTEGER, PUBLIC, PARAMETER ::   jpd14       = jp_lb + 3

   ! Starting/ending C14 do-loop indices (N.B. no C14 : jp_c14b0 > jp_c14b1 the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b0     = jp_lb     + 1            !: First index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b1     = jp_lb     + jp_c14b      !: Last  index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b0_2d  = jp_lb_2d  + 1            !: First index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b1_2d  = jp_lb_2d  + jp_c14b_2d   !: Last  index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b0_3d  = jp_lb_3d  + 1            !: First index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b1_3d  = jp_lb_3d  + jp_c14b_3d   !: Last  index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b0_trd = jp_lb_trd + 1            !: First index of C14 tracer
   INTEGER, PUBLIC, PARAMETER ::   jp_c14b1_trd = jp_lb_trd + jp_c14b_trd  !: Last  index of C14 tracer

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_c14b.F90 3680 2012-11-27 14:42:24Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE par_c14b
