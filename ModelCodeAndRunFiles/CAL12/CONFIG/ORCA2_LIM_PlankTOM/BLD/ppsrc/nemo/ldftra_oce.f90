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

MODULE ldftra_oce
   !!=====================================================================
   !!                      ***  MODULE  ldftra_oce  ***
   !! Ocean physics :  lateral tracer mixing coefficient defined in memory 
   !!=====================================================================
   !! History :  9.0  !  2002-11  (G. Madec)  Original code
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC ldftra_oce_alloc ! called by nemo_init->nemo_alloc, nemogcm.F90

   !!----------------------------------------------------------------------
   !! Lateral eddy diffusivity coefficients (tracers)
   !!----------------------------------------------------------------------
   !                                     !!* Namelist namtra_ldf : lateral mixing *
   LOGICAL , PUBLIC ::   ln_traldf_lap    !: laplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_bilap  !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_level  !: iso-level direction
   LOGICAL , PUBLIC ::   ln_traldf_hor    !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_traldf_iso    !: iso-neutral direction
   LOGICAL , PUBLIC ::   ln_traldf_grif   !: griffies skew flux
   LOGICAL , PUBLIC ::   ln_traldf_gdia   !: griffies skew flux streamfunction diagnostics
   REAL(wp), PUBLIC ::   rn_aht_0         !: lateral eddy diffusivity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahtb_0        !: lateral background eddy diffusivity (m2/s)
   REAL(wp), PUBLIC ::   rn_aeiv_0        !: eddy induced velocity coefficient (m2/s)
   REAL(wp), PUBLIC ::   rn_slpmax        !: slope limit
   REAL(wp), PUBLIC ::   rn_chsmag        !:  multiplicative factor in Smagorinsky diffusivity
   REAL(wp), PUBLIC ::   rn_smsh          !:  Smagorinsky diffusivity: = 0 - use only sheer
   REAL(wp), PUBLIC ::   rn_aht_m         !:  upper limit or stability criteria for lateral eddy diffusivity (m2/s)

   REAL(wp), PUBLIC ::   aht0, ahtb0, aeiv0         !!: OLD namelist names

   LOGICAL , PUBLIC ::   ln_triad_iso    !: calculate triads twice
   LOGICAL , PUBLIC ::   ln_botmix_grif  !: mixing on bottom
   LOGICAL , PUBLIC ::   l_grad_zps      = .FALSE.   !: special treatment for Horz Tgradients w partial steps 

   REAL(wp), PUBLIC ::   rldf                        !: multiplicative factor of diffusive coefficient
                                                     !: Needed to define the ratio between passive and active tracer diffusion coef. 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ahtt, ahtu, ahtv, ahtw   !: ** 2D coefficients ** at T-,U-,V-,W-points

   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'                              eddy induced velocity
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER               ::   lk_traldf_eiv   = .TRUE.   !: eddy induced velocity flag
   
   !                                                                              !!! eddy coefficients at U-, V-, W-points  [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   aeiu , aeiv , aeiw   !: ** 2D coefficients **
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   u_eiv, v_eiv, w_eiv   !: eddy induced velocity [m/s]


   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_oce.F90 4147 2013-11-04 11:51:55Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ldftra_oce_alloc()
     !!----------------------------------------------------------------------
      !!                 ***  FUNCTION ldftra_oce_alloc  ***
     !!----------------------------------------------------------------------
     INTEGER, DIMENSION(3) :: ierr
     !!----------------------------------------------------------------------
     ierr(:) = 0

      ALLOCATE( ahtt(jpi,jpj    ) , ahtu(jpi,jpj    ) , ahtv(jpi,jpj    ) , ahtw(jpi,jpj    ) , STAT=ierr(1) )
      !
      ALLOCATE( aeiu(jpi,jpj    ) , aeiv(jpi,jpj    ) , aeiw(jpi,jpj    ) , STAT=ierr(2) )
      ALLOCATE( u_eiv(jpi,jpj,jpk), v_eiv(jpi,jpj,jpk), w_eiv(jpi,jpj,jpk), STAT=ierr(3))
      ldftra_oce_alloc = MAXVAL( ierr )
      IF( ldftra_oce_alloc /= 0 )   CALL ctl_warn('ldftra_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION ldftra_oce_alloc

   !!=====================================================================
END MODULE ldftra_oce
