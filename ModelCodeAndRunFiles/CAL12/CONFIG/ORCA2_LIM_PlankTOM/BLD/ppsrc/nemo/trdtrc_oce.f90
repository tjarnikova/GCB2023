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

MODULE trdtrc_oce
   !!======================================================================
   !!                   ***  MODULE trdtrc_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   USE par_oce       ! ocean parameters
   USE par_trc       ! passive tracers parameters

   IMPLICIT NONE
   PUBLIC

   !                                         !!* Namelist namtoptrd:  diagnostics on passive tracers trends
   INTEGER  ::    nn_trd_trc                  !: time step frequency dynamics and tracers trends
   INTEGER  ::    nn_ctls_trc                 !: control surface type for trends vertical integration
   REAL(wp) ::    rn_ucf_trc                  !: unit conversion factor (for netCDF trends outputs)
   LOGICAL  ::    ln_trdmxl_trc_instant       !: flag to diagnose inst./mean ML trc trends
   LOGICAL  ::    ln_trdmxl_trc_restart       !: flag to restart mixed-layer trc diagnostics
   CHARACTER(len=50) ::  cn_trdrst_trc_in     !: suffix of pass. tracer restart name (input)
   CHARACTER(len=50) ::  cn_trdrst_trc_out    !: suffix of pass. tracer restart name (output)
   LOGICAL, DIMENSION(jptra) ::   ln_trdtrc   !: large trends diagnostic to write or not (namelist)

   LOGICAL, PARAMETER ::   lk_trdtrc = .FALSE.   !: ML trend flag

   LOGICAL, PARAMETER ::   lk_trdmxl_trc = .FALSE.   !: ML trend flag

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trdtrc_oce.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_trc_oce_alloc()
      !!----------------------------------------------------------------------
      !!         *** ROUTINE trd_trc_oce_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_warn
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      !
      !
      trd_trc_oce_alloc = MAXVAL(ierr)
      !
      IF( trd_trc_oce_alloc /= 0 )   CALL ctl_warn('trd_trc_oce_alloc: failed to allocate arrays')
      !
      !
   END FUNCTION trd_trc_oce_alloc


   !!======================================================================
END MODULE trdtrc_oce
