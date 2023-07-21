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

MODULE nc4interface
!-
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
      !!--------------------------------------------------------------------
      !! NOT 'key_netcdf4' Defines dummy routines for netcdf4
      !!                   calls when compiling without netcdf4 libraries
      !!--------------------------------------------------------------------
  !- netcdf4 chunking control structure
  !- (optional on histbeg and histend calls)
!$AGRIF_DO_NOT_TREAT
  TYPE, PUBLIC :: snc4_ctl
     SEQUENCE
     INTEGER :: ni
     INTEGER :: nj
     INTEGER :: nk
     LOGICAL :: luse
  END TYPE snc4_ctl
!$AGRIF_END_DO_NOT_TREAT

CONTAINS
!===
   SUBROUTINE GET_NF90_SYMBOL(sym_name, ivalue)
      CHARACTER(len=*),      INTENT(in)  :: sym_name
      INTEGER,               INTENT(out) :: ivalue
      ivalue = -999
   END SUBROUTINE GET_NF90_SYMBOL
   INTEGER FUNCTION SET_NF90_DEF_VAR_CHUNKING(idum1, idum2, idum3, iarr1)
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE NF90_DEF_VAR_CHUNKING  ***
      !!
      !! ** Purpose :   Dummy NetCDF4 routine to enable compiling with NetCDF3 libraries
      !!--------------------------------------------------------------------
      INTEGER,               INTENT(in) :: idum1, idum2, idum3
      INTEGER, DIMENSION(4), INTENT(in) :: iarr1
      WRITE(*,*) 'Warning: Attempt to chunk output variable without NetCDF4 support'
      SET_NF90_DEF_VAR_CHUNKING = -1
   END FUNCTION SET_NF90_DEF_VAR_CHUNKING

   INTEGER FUNCTION SET_NF90_DEF_VAR_DEFLATE(idum1, idum2, idum3, idum4, idum5)
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE NF90_DEF_VAR_DEFLATE  ***
      !!
      !! ** Purpose :   Dummy NetCDF4 routine to enable compiling with NetCDF3 libraries
      !!--------------------------------------------------------------------
      INTEGER,               INTENT(in) :: idum1, idum2, idum3, idum4, idum5
      WRITE(*,*) 'Warning: Attempt to compress output variable without NetCDF4 support'
      SET_NF90_DEF_VAR_DEFLATE = -1
   END FUNCTION SET_NF90_DEF_VAR_DEFLATE

!------------------
END MODULE nc4interface
