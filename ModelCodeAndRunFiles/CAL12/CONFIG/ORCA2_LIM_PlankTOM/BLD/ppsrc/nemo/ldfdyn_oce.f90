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

MODULE ldfdyn_oce
   !!======================================================================
   !!                  ***  MODULE  ldfdyn_oce  ***
   !! Ocean physics:  lateral momentum mixing coefficient defined in memory 
   !!======================================================================
   !! History :  1.0  ! 2002-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PUBLIC

   !                                       !!* Namelist namdyn_ldf : lateral mixing *
   LOGICAL , PUBLIC ::   ln_dynldf_lap      !: laplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_bilap    !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_level    !: iso-level direction
   LOGICAL , PUBLIC ::   ln_dynldf_hor      !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_dynldf_iso      !: iso-neutral direction
   REAL(wp), PUBLIC ::   rn_ahm_0_lap       !: lateral laplacian eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahmb_0          !: lateral laplacian background eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahm_0_blp       !: lateral bilaplacian eddy viscosity (m4/s)
   REAL(wp), PUBLIC ::   ahm0, ahmb0, ahm0_blp         !: OLD namelist names
   REAL(wp), PUBLIC ::   rn_cmsmag_1        !: constant in laplacian Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_cmsmag_2        !: constant in bilaplacian Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_cmsh            !: 1 or 0 , if 0 -use only shear for Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_ahm_m_blp       !: upper limit for bilap  abs(ahm) < min( dx^4/128rdt, rn_ahm_m_blp)
   REAL(wp), PUBLIC ::   rn_ahm_m_lap       !: upper limit for lap  ahm < min(dx^2/16rdt, rn_ahm_m_lap)

   INTEGER , PUBLIC ::   nkahm_smag      =  0          !: 

   !                                                                                  !!! eddy coeff. at U-,V-,W-pts [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahm1, ahm2, ahm3, ahm4   !: ** 3D coefficients **

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: ldfdyn_oce.F90 4147 2013-11-04 11:51:55Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ldfdyn_oce_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION ldfdyn_oce_alloc  ***
      !!----------------------------------------------------------------------
      ldfdyn_oce_alloc = 0
      ALLOCATE( ahm1(jpi,jpj,jpk) , ahm2(jpi,jpj,jpk) , ahm3(jpi,jpj,jpk) , ahm4(jpi,jpj,jpk) , STAT=ldfdyn_oce_alloc )
      IF( ldfdyn_oce_alloc /= 0 )   CALL ctl_warn('ldfdyn_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION ldfdyn_oce_alloc

   !!======================================================================
END MODULE ldfdyn_oce
