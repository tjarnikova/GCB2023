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

MODULE obs_rot_vel
   !!======================================================================
   !!                       ***  MODULE obs_rot_vel  ***
   !! Observation diagnostics: Read the velocity profile observations
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_rotvel : Rotate velocity data into N-S,E-W directorions
   !!----------------------------------------------------------------------
   !! * Modules used   
   USE wrk_nemo                 ! Memory Allocation
   USE par_kind                 ! Precision variables
   USE par_oce                  ! Ocean parameters
   USE in_out_manager           ! I/O manager
   USE dom_oce                  ! Ocean space and time domain variables
   USE obs_grid                 ! Grid search
   USE obs_utils                ! For error handling
   USE obs_profiles_def         ! Profile definitions
   USE obs_inter_h2d            ! Horizontal interpolation
   USE obs_inter_sup            ! MPP support routines for interpolation
   USE geo2ocean                ! Rotation of vectors
   USE obs_fbm                  ! Feedback definitions

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_rotvel            ! Rotate the observations

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_rot_vel.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_rotvel( profdata, k2dint, pu, pv )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_pro_dri ***
      !!
      !! ** Purpose : Rotate velocity data into N-S,E-W directorions
      !!
      !! ** Method  : Interpolation of geo2ocean coefficients on U,V grid
      !!              to observation point followed by a similar computations
      !!              as in geo2ocean.
      !!
      !! ** Action  : Review if there is a better way to do this.
      !!
      !! References : 
      !!
      !! History :  
      !!      ! :  2009-02 (K. Mogensen) : New routine
      !!----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: profdata    ! Profile data to be read
      INTEGER, INTENT(IN) :: k2dint     ! Horizontal interpolation methed
      REAL(wp), DIMENSION(*) :: &
         & pu, &
         & pv
      !! * Local declarations
      REAL(wp), DIMENSION(2,2,1) :: zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmasku, &
         & zmaskv, &
         & zcoslu, &
         & zsinlu, &
         & zcoslv, &
         & zsinlv, &
         & zglamu, &
         & zgphiu, &
         & zglamv, &
         & zgphiv
      REAL(wp), DIMENSION(1) :: &
         & zsinu, &
         & zcosu, &
         & zsinv, &
         & zcosv
      REAL(wp) :: zsin
      REAL(wp) :: zcos
      REAL(wp), DIMENSION(1) :: zobsmask
      REAL(wp), POINTER, DIMENSION(:,:) :: zsingu,zcosgu,zsingv,zcosgv
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdiu, &
         & igrdju, &
         & igrdiv, &
         & igrdjv
      INTEGER :: ji
      INTEGER :: jk

      CALL wrk_alloc(jpi,jpj,zsingu,zcosgu,zsingv,zcosgv) 

      !-----------------------------------------------------------------------
      ! Allocate data for message parsing and interpolation
      !-----------------------------------------------------------------------

      ALLOCATE( &
         & igrdiu(2,2,profdata%nprof), &
         & igrdju(2,2,profdata%nprof), &
         & zglamu(2,2,profdata%nprof), &
         & zgphiu(2,2,profdata%nprof), &
         & zmasku(2,2,profdata%nprof), &
         & zcoslu(2,2,profdata%nprof), &
         & zsinlu(2,2,profdata%nprof), &
         & igrdiv(2,2,profdata%nprof), &
         & igrdjv(2,2,profdata%nprof), &
         & zglamv(2,2,profdata%nprof), &
         & zgphiv(2,2,profdata%nprof), &
         & zmaskv(2,2,profdata%nprof), &
         & zcoslv(2,2,profdata%nprof), &
         & zsinlv(2,2,profdata%nprof)  &
         & )

      !-----------------------------------------------------------------------
      ! Receive the angles on the U and V grids.
      !-----------------------------------------------------------------------

      CALL obs_rot( zsingu, zcosgu, zsingv, zcosgv )

      DO ji = 1, profdata%nprof
         igrdiu(1,1,ji) = profdata%mi(ji,1)-1
         igrdju(1,1,ji) = profdata%mj(ji,1)-1
         igrdiu(1,2,ji) = profdata%mi(ji,1)-1
         igrdju(1,2,ji) = profdata%mj(ji,1)
         igrdiu(2,1,ji) = profdata%mi(ji,1)
         igrdju(2,1,ji) = profdata%mj(ji,1)-1
         igrdiu(2,2,ji) = profdata%mi(ji,1)
         igrdju(2,2,ji) = profdata%mj(ji,1)
         igrdiv(1,1,ji) = profdata%mi(ji,2)-1
         igrdjv(1,1,ji) = profdata%mj(ji,2)-1
         igrdiv(1,2,ji) = profdata%mi(ji,2)-1
         igrdjv(1,2,ji) = profdata%mj(ji,2)
         igrdiv(2,1,ji) = profdata%mi(ji,2)
         igrdjv(2,1,ji) = profdata%mj(ji,2)-1
         igrdiv(2,2,ji) = profdata%mi(ji,2)
         igrdjv(2,2,ji) = profdata%mj(ji,2)
      END DO

      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiu, igrdju, &
         &                  glamu, zglamu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiu, igrdju, &
         &                  gphiu, zgphiu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiu, igrdju, &
         &                  umask(:,:,1), zmasku )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiu, igrdju, &
         &                  zsingu, zsinlu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiu, igrdju, &
         &                  zcosgu, zcoslu )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiv, igrdjv, &
         &                  glamv, zglamv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiv, igrdjv, &
         &                  gphiv, zgphiv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiv, igrdjv, &
         &                  vmask(:,:,1), zmaskv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiv, igrdjv, &
         &                  zsingv, zsinlv )
      CALL obs_int_comm_2d( 2, 2, profdata%nprof, igrdiv, igrdjv, &
         &                  zcosgv, zcoslv )

      DO ji = 1, profdata%nprof
            
         CALL obs_int_h2d_init( 1, 1, k2dint, &
            &                   profdata%rlam(ji), profdata%rphi(ji), &
            &                   zglamu(:,:,ji), zgphiu(:,:,ji), &
            &                   zmasku(:,:,ji), zweig, zobsmask )
         
         CALL obs_int_h2d( 1, 1, zweig, zsinlu(:,:,ji),  zsinu )

         CALL obs_int_h2d( 1, 1, zweig, zcoslu(:,:,ji),  zcosu )

         CALL obs_int_h2d_init( 1, 1, k2dint, &
            &                   profdata%rlam(ji), profdata%rphi(ji), &
            &                   zglamv(:,:,ji), zgphiv(:,:,ji), &
            &                   zmaskv(:,:,ji), zweig, zobsmask )
         
         CALL obs_int_h2d( 1, 1, zweig, zsinlv(:,:,ji),  zsinv )

         CALL obs_int_h2d( 1, 1, zweig, zcoslv(:,:,ji),  zcosv )

         ! Assume that the angle at observation point is the 
         ! mean of u and v cosines/sines

         zcos = 0.5_wp * ( zcosu(1) + zcosv(1) )
         zsin = 0.5_wp * ( zsinu(1) + zsinv(1) )
         
         IF ( ( profdata%npvsta(ji,1) /= profdata%npvsta(ji,2) ) .OR. &
            & ( profdata%npvend(ji,1) /= profdata%npvend(ji,2) ) ) THEN
            CALL fatal_error( 'Different number of U and V observations '// &
               'in a profile in obs_rotvel', 192 )
         ENDIF

         DO jk = profdata%npvsta(ji,1), profdata%npvend(ji,1)
            IF ( ( profdata%var(1)%vmod(jk) /= fbrmdi ) .AND. &
               & ( profdata%var(1)%vmod(jk) /= fbrmdi ) ) THEN
               pu(jk) = profdata%var(1)%vmod(jk) * zcos - &
                  &     profdata%var(2)%vmod(jk) * zsin 
               pv(jk) = profdata%var(2)%vmod(jk) * zcos + &
                  &     profdata%var(1)%vmod(jk) * zsin
            ELSE
               pu(jk) = fbrmdi
               pv(jk) = fbrmdi
            ENDIF
         END DO

      END DO
      
      DEALLOCATE( &
         & igrdiu, &
         & igrdju, &
         & zglamu, &
         & zgphiu, &
         & zmasku, &
         & zcoslu, &
         & zsinlu, &
         & igrdiv, &
         & igrdjv, &
         & zglamv, &
         & zgphiv, &
         & zmaskv, &
         & zcoslv, &
         & zsinlv  &
         & )

      CALL wrk_dealloc(jpi,jpj,zsingu,zcosgu,zsingv,zcosgv) 

   END SUBROUTINE obs_rotvel

END MODULE obs_rot_vel
