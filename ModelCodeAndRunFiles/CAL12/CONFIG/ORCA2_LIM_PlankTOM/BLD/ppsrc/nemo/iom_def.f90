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

MODULE iom_def
   !!=====================================================================
   !!                    ***  MODULE  iom_def ***
   !! IOM variables definitions
   !!====================================================================
   !! History :  9.0  ! 06 09  (S. Masson) Original code
   !!             "   ! 07 07  (D. Storkey) Add uldname
   !!--------------------------------------------------------------------
   !!---------------------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: iom_def.F90 4205 2013-11-14 17:02:06Z clem $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!---------------------------------------------------------------------------------

   USE par_kind

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC ::   jpdom_data          = 1   !: ( 1  :jpidta, 1  :jpjdta)
   INTEGER, PARAMETER, PUBLIC ::   jpdom_global        = 2   !: ( 1  :jpiglo, 1  :jpjglo)
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local         = 3   !: One of the 3 following cases
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local_full    = 4   !: ( 1  :jpi   , 1  :jpi   )
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local_noextra = 5   !: ( 1  :nlci  , 1  :nlcj  )
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local_noovlap = 6   !: (nldi:nlei  ,nldj:nlej  )
   INTEGER, PARAMETER, PUBLIC ::   jpdom_unknown       = 7   !: No dimension checking
   INTEGER, PARAMETER, PUBLIC ::   jpdom_autoglo       = 8   !: 
   INTEGER, PARAMETER, PUBLIC ::   jpdom_autodta       = 9   !: 

   INTEGER, PARAMETER, PUBLIC ::   jpioipsl    = 100      !: Use ioipsl (fliocom only) library
   INTEGER, PARAMETER, PUBLIC ::   jpnf90      = 101      !: Use nf90 library
   INTEGER, PARAMETER, PUBLIC ::   jprstdimg   = 102      !: Use restart dimgs (fortran direct acces) library
   INTEGER, PARAMETER, PUBLIC ::   jprstlib  = jpnf90     !: restarts io library

   INTEGER, PARAMETER, PUBLIC ::   jp_r8    = 200      !: write REAL(8)
   INTEGER, PARAMETER, PUBLIC ::   jp_r4    = 201      !: write REAL(4)
   INTEGER, PARAMETER, PUBLIC ::   jp_i4    = 202      !: write INTEGER(4)
   INTEGER, PARAMETER, PUBLIC ::   jp_i2    = 203      !: write INTEGER(2)
   INTEGER, PARAMETER, PUBLIC ::   jp_i1    = 204      !: write INTEGER(1)

   INTEGER, PARAMETER, PUBLIC ::   jpmax_files  = 100   !: maximum number of simultaneously opened file
   INTEGER, PARAMETER, PUBLIC ::   jpmax_vars   = 600  !: maximum number of variables in one file
   INTEGER, PARAMETER, PUBLIC ::   jpmax_dims   =  4   !: maximum number of dimensions for one variable
   INTEGER, PARAMETER, PUBLIC ::   jpmax_digits =  5   !: maximum number of digits for the cpu number in the file name

!$AGRIF_DO_NOT_TREAT
   INTEGER, PUBLIC            ::   iom_open_init = 0   !: used to initialize iom_file(:)%nfid to 0

   TYPE, PUBLIC ::   file_descriptor
      CHARACTER(LEN=240)                        ::   name     !: name of the file
      INTEGER                                   ::   nfid     !: identifier of the file (0 if closed)
      INTEGER                                   ::   iolib    !: library used to read the file (jpioipsl, jpnf90 or jprstdimg)
      INTEGER                                   ::   nvars    !: number of identified varibles in the file
      INTEGER                                   ::   iduld    !: id of the unlimited dimension
      INTEGER                                   ::   irec     !: writing record position  
      CHARACTER(LEN=32)                         ::   uldname  !: name of the unlimited dimension
      CHARACTER(LEN=32), DIMENSION(jpmax_vars)  ::   cn_var   !: names of the variables
      INTEGER, DIMENSION(jpmax_vars)            ::   nvid     !: id of the variables
      INTEGER, DIMENSION(jpmax_vars)            ::   ndims    !: number of dimensions of the variables
      LOGICAL, DIMENSION(jpmax_vars)            ::   luld     !: variable using the unlimited dimension
      INTEGER, DIMENSION(jpmax_dims,jpmax_vars) ::   dimsz    !: size of variables dimensions 
      REAL(kind=wp), DIMENSION(jpmax_vars)      ::   scf      !: scale_factor of the variables
      REAL(kind=wp), DIMENSION(jpmax_vars)      ::   ofs      !: add_offset of the variables
   END TYPE file_descriptor
   TYPE(file_descriptor), DIMENSION(jpmax_files), PUBLIC ::   iom_file !: array containing the info for all opened files
!$AGRIF_END_DO_NOT_TREAT

   !!=====================================================================
END MODULE iom_def
