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

MODULE par_my_trc
   !!======================================================================
   !!                        ***  par_my_trc  ***
   !! TOP :   set the MY_TRC parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_my_trc.F90 3680 2012-11-27 14:42:24Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   USE par_pisces , ONLY : jp_pisces       !: number of tracers in PISCES
   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

   USE par_cfc    , ONLY : jp_cfc          !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_2d       !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_3d       !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_trd      !: number of tracers in CFC

   USE par_c14b   , ONLY : jp_c14b         !: number of tracers in C14
   USE par_c14b   , ONLY : jp_c14b_2d      !: number of tracers in C14
   USE par_c14b   , ONLY : jp_c14b_3d      !: number of tracers in C14
   USE par_c14b   , ONLY : jp_c14b_trd     !: number of tracers in C14

   IMPLICIT NONE

   INTEGER, PARAMETER ::   jp_lm      =  jp_pisces     + jp_cfc     + jp_c14b     !: 
   INTEGER, PARAMETER ::   jp_lm_2d   =  jp_pisces_2d  + jp_cfc_2d  + jp_c14b_2d  !:
   INTEGER, PARAMETER ::   jp_lm_3d   =  jp_pisces_3d  + jp_cfc_3d  + jp_c14b_3d  !:
   INTEGER, PARAMETER ::   jp_lm_trd  =  jp_pisces_trd + jp_cfc_trd + jp_c14b_trd !:

   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (MY_TRC)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_my_trc     = .FALSE.  !: MY_TRC flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_my_trc     =  0       !: No MY_TRC tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_my_trc_2d  =  0       !: No MY_TRC additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_my_trc_3d  =  0       !: No MY_TRC additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_my_trc_trd =  0       !: number of sms trends for MY_TRC

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_myt0     = jp_lm     + 1              !: First index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt1     = jp_lm     + jp_my_trc      !: Last  index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt0_2d  = jp_lm_2d  + 1              !: First index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt1_2d  = jp_lm_2d  + jp_my_trc_2d   !: Last  index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt0_3d  = jp_lm_3d  + 1              !: First index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt1_3d  = jp_lm_3d  + jp_my_trc_3d   !: Last  index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt0_trd = jp_lm_trd + 1              !: First index of MY_TRC passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_myt1_trd = jp_lm_trd + jp_my_trc_trd  !: Last  index of MY_TRC passive tracers

   !!======================================================================
END MODULE par_my_trc
