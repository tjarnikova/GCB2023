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

MODULE dynspg_oce
   !!======================================================================
   !!                       ***  MODULE dynspg_oce  ***
   !!       
   !! Ocean dynamics: Define in memory some surface pressure gradient variables
   !!======================================================================
   !! History :  1.0  ! 2005-12  (C. Talandier, G. Madec)  Original code
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PUBLIC           

   PUBLIC   dynspg_oce_alloc   ! called in dynspg.F90
   
   !                                                       !!! Surface pressure gradient logicals
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .FALSE.  !: Explicit free surface flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .FALSE.  !: Free surface with time splitting flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .TRUE.   !: Filtered free surface cst volume flag

  !                                                                         !!! Time splitting scheme (key_dynspg_ts) 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sshn_e, ssha_e   ! sea surface heigth (now, after)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ua_e  , va_e     ! barotropic velocities (after)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hu_e  , hv_e     ! now ocean depth ( = Ho+sshn_e )
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hur_e , hvr_e    ! inverse of hu_e and hv_e
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_adv, vn_adv   ! Advection vel. at "now" barocl. step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ub2_b,  vb2_b    ! Half step fluxes (ln_bt_fw=T)

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: dynspg_oce.F90 4486 2014-02-05 11:23:56Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dynspg_oce_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dynspg_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( sshn_e(jpi,jpj) , ua_e(jpi,jpj) , hu_e(jpi,jpj) , hur_e(jpi,jpj) ,      &
         &      ssha_e(jpi,jpj) , va_e(jpi,jpj) , hv_e(jpi,jpj) , hvr_e(jpi,jpj) ,      &
         &      ub2_b(jpi,jpj)  , vb2_b(jpi,jpj)                                 ,      &
         &      un_adv(jpi,jpj) , vn_adv(jpi,jpj)                                , STAT = dynspg_oce_alloc )
         !
      IF( lk_mpp                )   CALL mpp_sum ( dynspg_oce_alloc )
      IF( dynspg_oce_alloc /= 0 )   CALL ctl_warn('dynspg_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION dynspg_oce_alloc

   !!======================================================================
END MODULE dynspg_oce
