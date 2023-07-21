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

MODULE dom_ice_2
   !!======================================================================
   !!                   ***  MODULE  dom_ice  ***
   !! LIM 2.0 Sea Ice :   Domain  variables
   !!======================================================================
   !! History :   2.0  !  03-08  (C. Ethe)  Free form and module
   !!             3.3  !  2009-05 (G.Garric, C. Bricaud) addition of lim2_evp case
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2'                                       LIM2 sea-ice model
   !!----------------------------------------------------------------------
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: dom_ice_2.F90 3764 2013-01-23 14:33:04Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   USE par_ice_2   ! LIM parameters

   IMPLICIT NONE
   PRIVATE

   PUBLIC    dom_ice_alloc_2    ! Called from nemogcm.F90

   LOGICAL, PUBLIC ::   l_jeq     = .TRUE.     !: Equator inside the domain flag

   INTEGER, PUBLIC ::   njeq , njeqm1          !: j-index of the equator if it is inside the domain
      !                                        !  (otherwise = jpj+10 (SH) or -10 (SH) )

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   fs2cor , fcor     !: coriolis factor and coeficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   covrai            !: sine of geographic latitude
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   area              !: surface of grid cell 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   tms    , tmu      !: temperature and velocity points masks
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)     ::   wght              !: weight of the 4 neighbours to compute averages
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   tmv               !: y-velocity mask used for evp rheology 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)     ::   akappa , bkappa   !: first and third group of metric coefficients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:,:) ::   alambd            !: second group of metric coefficients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   tmf               !: F-points masks
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)         ::   tmi               !: ice mask: =1 if ice thick > 0
   !!----------------------------------------------------------------------
   CONTAINS

   INTEGER FUNCTION dom_ice_alloc_2()
      !!----------------------------------------------------------------------
      USE lib_mpp, ONLY:   ctl_warn   ! MPP library
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( fs2cor(jpi,jpj)     , fcor(jpi,jpj) ,                                   &
         &      covrai(jpi,jpj)     , area(jpi,jpj) , tms(jpi,jpj) , tmu(jpi,jpj) ,     &
         &      wght  (jpi,jpj,2,2)                                               , STAT=ierr(1) )
         !
      ALLOCATE(                                                    &
         &        tmv(jpi,jpj) , tmf(jpi,jpj) , tmi(jpi,jpj) ,     &
         &        STAT=ierr(2) )
         !
      dom_ice_alloc_2 = MAXVAL(ierr)
      IF( dom_ice_alloc_2 /= 0 )   CALL ctl_warn('dom_ice_alloc_2: failed to allocate arrays')
      !
   END FUNCTION dom_ice_alloc_2

   !!======================================================================
END MODULE dom_ice_2
