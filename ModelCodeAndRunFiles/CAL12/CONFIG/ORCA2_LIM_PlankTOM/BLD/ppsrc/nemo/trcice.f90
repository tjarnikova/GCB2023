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

MODULE trcice
   !!======================================================================
   !!                         ***  MODULE trcice  ***
   !! TOP :   Manage the communication between TOP and sea ice
   !!======================================================================
   !! History :  3.5  ! 2013    (M. Vancoppenolle, O. Aumont, G. Madec), original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_ice   :  Call the appropriate sea ice tracer subroutine
   !!----------------------------------------------------------------------

   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables
   USE trcice_cfc      ! CFC      initialisation
   USE trcice_pisces   ! PISCES   initialisation
   USE trcice_my_trc   ! MY_TRC   initialisation
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_ice_ini ! called by trc_nam

CONTAINS
   
   SUBROUTINE trc_ice_ini
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ice_ini ***
      !!
      !! ** Purpose :   Initialization of the ice module for tracers
      !!
      !! ** Method  : - 
      !!            
      !!---------------------------------------------------------------------
      ! --- Variable declarations --- !

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_ice_ini : Initialize sea ice tracer boundary condition'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF

      IF( nn_timing == 1 )  CALL timing_start('trc_ice_ini')

      !
      trc_i(:,:,:) = 0.0d0 ! by default
      trc_o(:,:,:) = 0.0d0 ! by default

      IF ( nn_ice_tr == 1 ) THEN
         IF( lk_pisces  )    CALL trc_ice_ini_pisces       ! PISCES  bio-model
         IF( lk_cfc     )    CALL trc_ice_ini_cfc          ! CFC     tracers
         IF( lk_my_trc  )    CALL trc_ice_ini_my_trc       ! MY_TRC  tracers
      ENDIF

      IF( nn_timing == 1 )   CALL timing_stop('trc_ice_ini')
      !
   END SUBROUTINE trc_ice_ini


   !!======================================================================
END MODULE trcice
