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

MODULE icbstp

   !!======================================================================
   !!                       ***  MODULE  icbstp  ***
   !! Icebergs:  initialise variables for iceberg tracking
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!            -    !                            Move budgets to icbdia routine
   !!            -    !  2011-05  (Alderson)       Add call to copy forcing arrays
   !!            -    !                            into icb copies with haloes
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icb_stp       : start iceberg tracking
   !!   icb_end       : end   iceberg tracking
   !!----------------------------------------------------------------------
   USE par_oce        ! nemo parameters
   USE dom_oce        ! ocean domain
   USE sbc_oce        ! ocean surface forcing
   USE phycst
   USE in_out_manager ! nemo IO
   USE lib_mpp
   USE iom
   USE fldread
   USE timing         ! timing

   USE icb_oce        ! define iceberg arrays
   USE icbini         ! iceberg initialisation routines
   USE icbutl         ! iceberg utility routines
   USE icbrst         ! iceberg restart routines
   USE icbdyn         ! iceberg dynamics (ie advection) routines
   USE icbclv         ! iceberg calving routines
   USE icbthm         ! iceberg thermodynamics routines
   USE icblbc         ! iceberg lateral boundary routines (including mpp)
   USE icbtrj         ! iceberg trajectory I/O routines
   USE icbdia         ! iceberg budget

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_stp        ! routine called in sbcmod.F90 module
   PUBLIC   icb_end        ! routine called in nemogcm.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2011)
   !! $Id: icbstp.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_stp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_stp  ***
      !!
      !! ** Purpose :   iceberg time stepping.
      !!
      !! ** Method  : - top level routine to do things in the correct order
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      !
      LOGICAL ::   ll_sample_traj, ll_budget, ll_verbose   ! local logical
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('icb_stp')

      !! start of timestep housekeeping

      nktberg = kt

      IF( nn_test_icebergs < 0 ) THEN      ! read calving data
         !
         CALL fld_read ( kt, 1, sf_icb )
         src_calving(:,:)      = sf_icb(1)%fnow(:,:,1)    ! calving in km^3/year (water equivalent)
         src_calving_hflx(:,:) = 0._wp                    ! NO heat flux for now
         !
      ENDIF

      berg_grid%floating_melt(:,:) = 0._wp

      ! anything that needs to be reset to zero each timestep for budgets is dealt with here
      CALL icb_dia_step()

      ll_verbose = .FALSE.
      IF( nn_verbose_write > 0 .AND. &
          MOD(kt-1,nn_verbose_write ) == 0 )   ll_verbose = nn_verbose_level >= 0

      ! write out time
      IF( ll_verbose ) WRITE(numicb,9100) nktberg, ndastp, nsec_day
 9100 FORMAT('kt= ',i8, ' day= ',i8,' secs=',i8)

      ! copy nemo forcing arrays into iceberg versions with extra halo
      ! only necessary for variables not on T points
      CALL icb_utl_copy()

      !!----------------------------------------------------------------------
      !! process icebergs

                                     CALL icb_clv_flx( kt )   ! Accumulate ice from calving

                                     CALL icb_clv()           ! Calve excess stored ice into icebergs


!                               !==  For each berg, evolve  ==!
      !
      IF( ASSOCIATED(first_berg) )   CALL icb_dyn( kt )       ! ice berg dynamics

      IF( lk_mpp ) THEN          ;   CALL icb_lbc_mpp()       ! Send bergs to other PEs
      ELSE                       ;   CALL icb_lbc()           ! Deal with any cyclic boundaries in non-mpp case
      ENDIF

      IF( ASSOCIATED(first_berg) )   CALL icb_thm( kt )       ! Ice berg thermodynamics (melting) + rolling

      !!----------------------------------------------------------------------
      !! end of timestep housekeeping

      ll_sample_traj = .FALSE.
      IF( nn_sample_rate > 0 .AND. MOD(kt-1,nn_sample_rate) == 0 )   ll_sample_traj = .TRUE.
      IF( ll_sample_traj .AND.   &
          ASSOCIATED(first_berg) )   CALL icb_trj_write( kt )  ! For each berg, record trajectory

      ! Gridded diagnostics
      ! To get these iom_put's and those preceding to actually do something
      ! use 1 in cpp file and create content for XML file

      CALL iom_put( "calving"           , berg_grid%calving      (:,:)   )  ! 'calving mass input'
      CALL iom_put( "berg_floating_melt", berg_grid%floating_melt(:,:)   )  ! 'Melt rate of icebergs + bits' , 'kg/m2/s'
      CALL iom_put( "berg_stored_ice"   , berg_grid%stored_ice   (:,:,:) )  ! 'Accumulated ice mass by class', 'kg'

      ! store mean budgets
      CALL icb_dia_put()

      ! Dump icebergs to screen
      if ( nn_verbose_level >= 2 )   CALL icb_utl_print( 'icb_stp, status', kt )

      ! Diagnose budgets
      ll_budget = .FALSE.
      IF( nn_verbose_write > 0 .AND. MOD(kt-1,nn_verbose_write) == 0 )   ll_budget = ln_bergdia
      CALL icb_dia( ll_budget )

      IF( MOD(kt,nn_stock) == 0 ) THEN
         CALL icb_rst_write( kt )
         IF( nn_sample_rate > 0 )   CALL icb_trj_sync()
      ENDIF

      IF( nn_timing == 1 ) CALL timing_stop('icb_stp')
      !
   END SUBROUTINE icb_stp


   SUBROUTINE icb_end( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_end  ***
      !!
      !! ** Purpose :   close iceberg files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  :: kt
      !!----------------------------------------------------------------------

      ! only write a restart if not done in icb_stp
      IF( MOD(kt,nn_stock) .NE. 0 ) CALL icb_rst_write( kt )

      ! finish with trajectories if they were written
      IF( nn_sample_rate .GT. 0 ) CALL icb_trj_end()

      IF(lwp)   WRITE(numout,'(a,i6)') 'icebergs: icb_end complete', narea
      CALL flush( numicb )
      CLOSE( numicb )
      !
   END SUBROUTINE icb_end

   !!-------------------------------------------------------------------------

END MODULE icbstp
