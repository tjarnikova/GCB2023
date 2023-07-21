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

MODULE asmpar
   !!======================================================================
   !!                       ***  MODULE asmpar  ***
   !! Assimilation increment : Parameters for assimilation interface
   !!======================================================================

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   !! * Shared Modules variables
   CHARACTER (LEN=40), PUBLIC, PARAMETER :: &
      & c_asmbkg = 'assim_background_state_Jb',  & !: Filename for storing the 
                                                   !: background state for use 
                                                   !: in the Jb term
      & c_asmdin = 'assim_background_state_DI',  & !: Filename for storing the 
                                                   !: background state for direct 
                                                   !: initialization
      & c_asmtrj = 'assim_trj',                  & !: Filename for storing the 
                                                   !: reference trajectory
      & c_asminc = 'assim_background_increments'   !: Filename for storing the 
                                                   !: increments to the background
                                                   !: state

   INTEGER, PUBLIC :: nitbkg_r      !: Background time step referenced to nit000
   INTEGER, PUBLIC :: nitdin_r      !: Direct Initialization time step referenced to nit000
   INTEGER, PUBLIC :: nitiaustr_r   !: IAU starting time step referenced to nit000
   INTEGER, PUBLIC :: nitiaufin_r   !: IAU final time step referenced to nit000
   INTEGER, PUBLIC :: nittrjfrq     !: Frequency of trajectory output for 4D-VAR

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: asmpar.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

END MODULE asmpar
