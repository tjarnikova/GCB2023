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

MODULE obs_profiles
   !!=====================================================================
   !!                       ***  MODULE  obs_profiles  ***
   !! Observation diagnostics: Storage space for profile observations
   !!                          arrays and additional flags etc.
   !!=====================================================================
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_profiles.F90 2733 2011-04-08 15:55:31Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   
   !! * Modules used 
   USE obs_profiles_def ! Definition of profile data types and tools

   IMPLICIT NONE
   
   SAVE

   !! * Routine accessibility
   PRIVATE

   PUBLIC nprofsets, nprofvars, nprofextr, profdata, prodatqc
   PUBLIC nvelosets, nvelovars, nveloextr, velodata, veldatqc
   
   !! * Shared Module variables
   INTEGER :: nprofsets                    ! Total number of profile data sets
   INTEGER :: nprofvars                    ! Total number of variables for profiles
   INTEGER :: nprofextr                    ! Extra fields for each variable
   TYPE(obs_prof), POINTER ::  profdata(:) ! Initial profile data
   TYPE(obs_prof), POINTER ::  prodatqc(:) ! Profile data after quality control

   INTEGER :: nvelosets                     ! Total number of velocity profile data sets
   INTEGER :: nvelovars                     ! Total number of variables for profiles
   INTEGER :: nveloextr                     ! Extra fields for each variable
   TYPE(obs_prof), POINTER ::  velodata(:)  ! Initial velocity profile data
   TYPE(obs_prof), POINTER ::  veldatqc(:)  ! Velocity profile data after quality control
END MODULE obs_profiles
