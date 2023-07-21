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

MODULE zdf_oce
   !!======================================================================
   !!              ***  MODULE  zdf_oce  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================
   !! history :  1.0  !  2002-06  (G. Madec) Original code
   !!            3.2  !  2009-07  (G.Madec) addition of avm
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC  zdf_oce_alloc    ! Called in nemogcm.F90

   LOGICAL, PARAMETER, PUBLIC ::   lk_zdfcst        = .FALSE.        !: constant vertical mixing flag

   !                                 !!* namelist namzdf: vertical diffusion *
   REAL(wp), PUBLIC ::   rn_avm0     !: vertical eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_avt0     !: vertical eddy diffusivity (m2/s)
   INTEGER , PUBLIC ::   nn_avb      !: constant or profile background on avt (=0/1)
   INTEGER , PUBLIC ::   nn_havtb    !: horizontal shape or not for avtb (=0/1)
   LOGICAL , PUBLIC ::   ln_zdfexp   !: explicit vertical diffusion scheme flag
   INTEGER , PUBLIC ::   nn_zdfexp   !: number of sub-time step (explicit time stepping)
   LOGICAL , PUBLIC ::   ln_zdfevd   !: convection: enhanced vertical diffusion flag
   INTEGER , PUBLIC ::   nn_evdm     !: =0/1 flag to apply enhanced avm or not
   REAL(wp), PUBLIC ::   rn_avevd    !: vertical eddy coeff. for enhanced vert. diff. (m2/s)
   LOGICAL , PUBLIC ::   ln_zdfnpc   !: convection: non-penetrative convection flag
   INTEGER , PUBLIC ::   nn_npc      !: non penetrative convective scheme call  frequency
   INTEGER , PUBLIC ::   nn_npcp     !: non penetrative convective scheme print frequency


   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:)     ::   avmb , avtb    !: background profile of avm and avt
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   avtb_2d        !: horizontal shape of background Kz profile
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   bfrua, bfrva   !: Bottom friction coefficients set in zdfbfr
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   tfrua, tfrva   !: top friction coefficients set in zdfbfr
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avmu , avmv    !: vertical viscosity coef at uw- & vw-pts       [m2/s]
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avm  , avt     !: vertical viscosity & diffusivity coef at w-pt [m2/s]
 
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdf_oce.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_oce_alloc()
      !!----------------------------------------------------------------------
      !!            *** FUNCTION zdf_oce_alloc ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE(avmb(jpk) , bfrua(jpi,jpj) ,                         &
         &     avtb(jpk) , bfrva(jpi,jpj) , avtb_2d(jpi,jpj) ,      &
         &     tfrua(jpi, jpj), tfrva(jpi, jpj)              ,      &
         &     avmu(jpi,jpj,jpk), avm(jpi,jpj,jpk)           ,      &
         &     avmv(jpi,jpj,jpk), avt(jpi,jpj,jpk)           , STAT = zdf_oce_alloc )
         !
      IF( zdf_oce_alloc /= 0 )   CALL ctl_warn('zdf_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION zdf_oce_alloc

   !!======================================================================
END MODULE zdf_oce
