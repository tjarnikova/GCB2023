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

MODULE limtab_2
   !!======================================================================
   !!                       ***  MODULE limtab_2   ***
   !!   LIM : transform 1D (2D) array to a 2D (1D) table
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   tab_2d_1d  : 2-D <==> 1-D
   !!   tab_1d_2d  : 1-D <==> 2-D
   !!----------------------------------------------------------------------
   USE par_kind

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tab_2d_1d_2   ! called by limthd
   PUBLIC   tab_1d_2d_2   ! called by limthd

   !!----------------------------------------------------------------------
   !! NEMO/LIM2 4.0 , UCL - NEMO Consortium (2010)
   !! $Id: limtab_2.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tab_2d_1d_2 ( ndim1d, tab1d, tab2d, ndim2d_x, ndim2d_y, tab_ind )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                               , INTENT(in   ) ::   ndim1d, ndim2d_x, ndim2d_y   ! 1D & 2D sizes
      REAL(wp), DIMENSION(ndim2d_x,ndim2d_y), INTENT(in   ) ::   tab2d                        ! input 2D field
      INTEGER , DIMENSION(ndim1d)           , INTENT(in   ) ::   tab_ind                      ! input index
      REAL(wp), DIMENSION(ndim1d)           , INTENT(  out) ::   tab1d                        ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , ndim2d_x ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / ndim2d_x + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO 
   END SUBROUTINE tab_2d_1d_2


   SUBROUTINE tab_1d_2d_2 ( ndim1d, tab2d, tab_ind, tab1d, ndim2d_x, ndim2d_y )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                               , INTENT(in   ) ::   ndim1d, ndim2d_x, ndim2d_y   ! 1d & 2D sizes
      REAL(wp), DIMENSION(ndim1d)           , INTENT(in   ) ::   tab1d                        ! input 1D field
      INTEGER , DIMENSION(ndim1d)           , INTENT(in   ) ::   tab_ind                      ! input index
      REAL(wp), DIMENSION(ndim2d_x,ndim2d_y), INTENT(  out) ::   tab2d                        ! output 2D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1 , ndim2d_x ) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / ndim2d_x  + 1
         tab2d(jid, jjd) = tab1d( jn)
      END DO
   END SUBROUTINE tab_1d_2d_2

   !!======================================================================
END MODULE limtab_2
