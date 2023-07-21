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

MODULE zdfkpp
   !!======================================================================
   !!                       ***  MODULE  zdfkpp  ***
   !! Ocean physics:  vertical mixing coefficient compute from the KPP 
   !!                 turbulent closure parameterization
   !!=====================================================================
   !! History :  OPA  ! 2000-03 (W.G. Large, J. Chanut) Original code
   !!            8.1  ! 2002-06 (J.M. Molines) for real case CLIPPER  
   !!            8.2  ! 2003-10 (Chanut J.) re-writting
   !!   NEMO     1.0  ! 2005-01 (C. Ethe, G. Madec) Free form, F90 + creation of tra_kpp routine
   !!            3.3  ! 2010-10 (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO KPP scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfkpp = .FALSE.   !: KPP flag
CONTAINS
   SUBROUTINE zdf_kpp_init           ! Dummy routine
      WRITE(*,*) 'zdf_kpp_init: You should not have seen this print! error?'
   END SUBROUTINE zdf_kpp_init
   SUBROUTINE zdf_kpp( kt )          ! Dummy routine
      WRITE(*,*) 'zdf_kpp: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_kpp
   SUBROUTINE tra_kpp( kt )          ! Dummy routine
      WRITE(*,*) 'tra_kpp: You should not have seen this print! error?', kt
   END SUBROUTINE tra_kpp
   SUBROUTINE trc_kpp( kt )          ! Dummy routine
      WRITE(*,*) 'trc_kpp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_kpp

   !!======================================================================
END MODULE zdfkpp
