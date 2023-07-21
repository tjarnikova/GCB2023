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

MODULE par_trc
   !!======================================================================
   !!                        ***  par_trc  ***
   !! TOP :   set the passive tracers parameters
   !!======================================================================
   !! History :    -   !  1996-01  (M. Levy)  original code
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!             1.0  !  2004-03  (C. Ethe) Free form and module
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters
   !
   USE par_pisces    ! PISCES  model
   USE par_c14b      ! C14 bomb tracer
   USE par_cfc       ! CFC 11 and 12 tracers
   USE par_my_trc    ! user defined passive tracers
   USE par_planktom  ! PlankTOM bgc

   IMPLICIT NONE

   ! Passive tracers : Maximum number of tracers. Needed to define data structures
   ! --------------- 
   INTEGER, PUBLIC,  PARAMETER ::   jpmaxtrc = 100

   ! Passive tracers : Total size
   ! ---------------               ! total number of passive tracers, of 2d and 3d output and trend arrays
   INTEGER, PUBLIC,  PARAMETER ::   jptra    =  jp_pisces     + jp_cfc     + jp_c14b    + jp_planktom
   INTEGER, PUBLIC,  PARAMETER ::   jpdia2d  =  jp_pisces_2d  + jp_cfc_2d  + jp_c14b_2d + jp_my_trc_2d
   INTEGER, PUBLIC,  PARAMETER ::   jpdia3d  =  jp_pisces_3d  + jp_cfc_3d  + jp_c14b_3d + jp_my_trc_3d
   !                     ! total number of sms diagnostic arrays
   INTEGER, PUBLIC,  PARAMETER ::   jpdiabio =  jp_pisces_trd + jp_cfc_trd + jp_c14b_trd + jp_my_trc_trd
   
   !  1D configuration ("key_c1d")
   ! -----------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trc_c1d   = .FALSE.  !: 1D pass. tracer configuration flag

   REAL(wp), PUBLIC  :: rtrn  = 0.5 * EPSILON( 1.e0 )    !: truncation value

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_trc.F90 4529 2014-03-15 11:00:04Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE par_trc
