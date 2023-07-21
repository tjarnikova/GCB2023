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

MODULE diadct
  !!=====================================================================
  !!                       ***  MODULE  diadct  ***
  !! Ocean diagnostics: Compute the transport trough a sec.
  !!===============================================================
  !! History : 
  !!
  !!         original  : 02/99 (Y Drillet)
  !!         addition  : 10/01 (Y Drillet, R Bourdalle Badie)
  !!                   : 10/05 (M Laborie) F90
  !!         addition  : 04/07 (G Garric) Ice sections
  !!         bugfix    : 04/07 (C Bricaud) test on sec%nb_point
  !!                                      initialisation of ztransp1,ztransp2,...
  !!         nemo_v_3_4: 09/2011 (C Bricaud)
  !!
  !!
  !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                       Dummy module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diadct = .FALSE.    !: diamht flag
   PUBLIC 
   !! $Id: diadct.F90 5505 2015-06-29 12:51:57Z timgraham $
CONTAINS

   SUBROUTINE dia_dct_init          ! Dummy routine
      WRITE(*,*) 'dia_dct_init: You should not have seen this print! error?', kt
   END SUBROUTINE dia_dct_init

   SUBROUTINE dia_dct( kt )         ! Dummy routine
      INTEGER, INTENT( in ) :: kt   ! ocean time-step index
      WRITE(*,*) 'dia_dct: You should not have seen this print! error?', kt
   END SUBROUTINE dia_dct

END MODULE diadct
