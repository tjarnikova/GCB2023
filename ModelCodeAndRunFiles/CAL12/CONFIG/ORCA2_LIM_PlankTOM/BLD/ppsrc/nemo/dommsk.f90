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

MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask and sup-
   !!                 !                      pression of the double computation of bmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of mbathy and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk        : compute land/ocean mask
   !!   dom_msk_nsa    : update land/ocean mask when no-slip accurate option is used.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk         ! routine called by inidom.F90
   PUBLIC   dom_msk_alloc   ! routine called by nemogcm.F90

   !                            !!* Namelist namlbc : lateral boundary condition *
   REAL(wp)        :: rn_shlat   ! type of lateral boundary condition on velocity
   LOGICAL, PUBLIC :: ln_vorlat  !  consistency of vorticity boundary condition 
   !                                            with analytical eqs.


   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  icoord ! Workspace for dom_msk_nsa()

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: vectopt_loop_substitute.h90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: dommsk.F90 5551 2015-07-06 16:38:51Z acc $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   INTEGER FUNCTION dom_msk_alloc()
      !!---------------------------------------------------------------------
      !!                 ***  FUNCTION dom_msk_alloc  ***
      !!---------------------------------------------------------------------
      dom_msk_alloc = 0
      IF( dom_msk_alloc /= 0 )   CALL ctl_warn('dom_msk_alloc: failed to allocate icoord array')
      !
   END FUNCTION dom_msk_alloc


   SUBROUTINE dom_msk
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) and baro-
      !!      tropic stream function  points (b).
      !!
      !! ** Method  :   The ocean/land mask is computed from the basin bathy-
      !!      metry in level (mbathy) which is defined or read in dommba.
      !!      mbathy equals 0 over continental T-point 
      !!      and the number of ocean level over the ocean.
      !!
      !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
      !!      t-point : 0. IF mbathy( ji ,jj) =< 0
      !!                1. IF mbathy( ji ,jj) >= jk
      !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
      !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
      !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
      !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
      !!                  and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
      !!      b-point : the same definition as for f-point of the first ocean
      !!                level (surface level) but with 0 along coastlines.
      !!      tmask_i : interior ocean mask at t-point, i.e. excluding duplicated
      !!                rows/lines due to cyclic or North Fold boundaries as well
      !!                as MPP halos.
      !!
      !!        The lateral friction is set through the value of fmask along
      !!      the coast and topography. This value is defined by rn_shlat, a
      !!      namelist parameter:
      !!         rn_shlat = 0, free slip  (no shear along the coast)
      !!         rn_shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_shlat, strong slip        | in the lateral boundary layer
      !!
      !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
      !!      are defined with the proper value at lateral domain boundaries,
      !!      but bmask. indeed, bmask defined the domain over which the
      !!      barotropic stream function is computed. this domain cannot
      !!      contain identical columns because the matrix associated with
      !!      the barotropic stream function equation is then no more inverti-
      !!      ble. therefore bmask is set to 0 along lateral domain boundaries
      !!      even IF nperio is not zero.
      !!
      !!      In case of open boundaries (lk_bdy=T):
      !!        - tmask is set to 1 on the points to be computed bay the open
      !!          boundaries routines.
      !!        - bmask is  set to 0 on the open boundaries.
      !!
      !! ** Action :   tmask    : land/ocean mask at t-point (=0. or 1.)
      !!               umask    : land/ocean mask at u-point (=0. or 1.)
      !!               vmask    : land/ocean mask at v-point (=0. or 1.)
      !!               fmask    : land/ocean mask at f-point (=0. or 1.)
      !!                          =rn_shlat along lateral boundaries
      !!               bmask    : land/ocean mask at barotropic stream
      !!                          function point (=0. or 1.) and set to 0 along lateral boundaries
      !!               tmask_i  : interior ocean mask
      !!----------------------------------------------------------------------
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iif, iil, ii0, ii1, ii   ! local integers
      INTEGER  ::   ijf, ijl, ij0, ij1       !   -       -
      INTEGER  ::   ios
      INTEGER  ::   isrow                    ! index for ORCA1 starting row
      INTEGER , POINTER, DIMENSION(:,:) ::  imsk
      REAL(wp), POINTER, DIMENSION(:,:) ::  zwf
      !!
      NAMELIST/namlbc/ rn_shlat, ln_vorlat
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_msk')
      !
      CALL wrk_alloc( jpi, jpj, imsk )
      CALL wrk_alloc( jpi, jpj, zwf  )
      !
      REWIND( numnam_ref )              ! Namelist namlbc in reference namelist : Lateral momentum boundary condition
      READ  ( numnam_ref, namlbc, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlbc in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namlbc in configuration namelist : Lateral momentum boundary condition
      READ  ( numnam_cfg, namlbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namlbc )
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '   Namelist namlbc'
         WRITE(numout,*) '      lateral momentum boundary cond.    rn_shlat  = ',rn_shlat
         WRITE(numout,*) '      consistency with analytical form   ln_vorlat = ',ln_vorlat 
      ENDIF

      IF     (      rn_shlat == 0.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  free-slip '
      ELSEIF (      rn_shlat == 2.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  no-slip '
      ELSEIF ( 0. < rn_shlat .AND. rn_shlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  partial-slip '
      ELSEIF ( 2. < rn_shlat                     ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  strong-slip '
      ELSE
         WRITE(ctmp1,*) ' rn_shlat is negative = ', rn_shlat
         CALL ctl_stop( ctmp1 )
      ENDIF

      ! 1. Ocean/land mask at t-point (computed from mbathy)
      ! -----------------------------
      ! N.B. tmask has already the right boundary conditions since mbathy is ok
      !
      tmask(:,:,:) = 0._wp
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( REAL( mbathy(ji,jj) - jk, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
            END DO  
         END DO  
      END DO  
      
      ! (ISF) define barotropic mask and mask the ice shelf point
      ssmask(:,:)=tmask(:,:,1) ! at this stage ice shelf is not masked
      
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( REAL( misfdep(ji,jj) - jk, wp ) - 0.1_wp >= 0._wp )   THEN
                  tmask(ji,jj,jk) = 0._wp
               END IF
            END DO  
         END DO  
      END DO  

!!gm  ????
!!gm end

      ! Interior domain mask (used for global sum)
      ! --------------------
      tmask_i(:,:) = ssmask(:,:)            ! (ISH) tmask_i = 1 even on the ice shelf
      iif = jpreci                         ! ???
      iil = nlci - jpreci + 1
      ijf = jprecj                         ! ???
      ijl = nlcj - jprecj + 1

      tmask_i( 1 :iif,   :   ) = 0._wp      ! first columns
      tmask_i(iil:jpi,   :   ) = 0._wp      ! last  columns (including mpp extra columns)
      tmask_i(   :   , 1 :ijf) = 0._wp      ! first rows
      tmask_i(   :   ,ijl:jpj) = 0._wp      ! last  rows (including mpp extra rows)

      ! north fold mask
      ! ---------------
      tpol(1:jpiglo) = 1._wp 
      fpol(1:jpiglo) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0._wp
         fpol(     1    :jpiglo) = 0._wp
         IF( mjg(nlej) == jpjglo ) THEN                  ! only half of the nlcj-1 row
            DO ji = iif+1, iil-1
               tmask_i(ji,nlej-1) = tmask_i(ji,nlej-1) * tpol(mig(ji))
            END DO
         ENDIF
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0._wp
         fpol(jpiglo/2+1:jpiglo) = 0._wp
      ENDIF

      ! 2. Ocean/land mask at u-,  v-, and z-points (computed from tmask)
      ! -------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector loop
               umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
               vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
            END DO
            DO ji = 1, jpim1      ! NO vector opt.
               fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                  &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO
      ! (ISF) MIN(1,SUM(umask)) is here to check if you have effectively at least 1 wet u point
      DO jj = 1, jpjm1
         DO ji = 1, jpim1   ! vector loop
            umask_i(ji,jj)  = ssmask(ji,jj) * ssmask(ji+1,jj  )  * MIN(1._wp,SUM(umask(ji,jj,:)))
            vmask_i(ji,jj)  = ssmask(ji,jj) * ssmask(ji  ,jj+1)  * MIN(1._wp,SUM(vmask(ji,jj,:)))
         END DO
         DO ji = 1, jpim1      ! NO vector opt.
            fmask_i(ji,jj) =  ssmask(ji,jj  ) * ssmask(ji+1,jj  )   &
               &            * ssmask(ji,jj+1) * ssmask(ji+1,jj+1) * MIN(1._wp,SUM(fmask(ji,jj,:)))
         END DO
      END DO
      CALL lbc_lnk( umask, 'U', 1._wp )      ! Lateral boundary conditions
      CALL lbc_lnk( vmask, 'V', 1._wp )
      CALL lbc_lnk( fmask, 'F', 1._wp )
      CALL lbc_lnk( umask_i, 'U', 1._wp )      ! Lateral boundary conditions
      CALL lbc_lnk( vmask_i, 'V', 1._wp )
      CALL lbc_lnk( fmask_i, 'F', 1._wp )

      ! 3. Ocean/land mask at wu-, wv- and w points 
      !----------------------------------------------
      wmask (:,:,1) = tmask(:,:,1) ! ????????
      wumask(:,:,1) = umask(:,:,1) ! ????????
      wvmask(:,:,1) = vmask(:,:,1) ! ????????
      DO jk=2,jpk
         wmask (:,:,jk)=tmask(:,:,jk) * tmask(:,:,jk-1)
         wumask(:,:,jk)=umask(:,:,jk) * umask(:,:,jk-1)   
         wvmask(:,:,jk)=vmask(:,:,jk) * vmask(:,:,jk-1)
      END DO

      ! 4. ocean/land mask for the elliptic equation
      ! --------------------------------------------
      bmask(:,:) = ssmask(:,:)       ! elliptic equation is written at t-point
      !
      !                               ! Boundary conditions
      !                                    ! cyclic east-west : bmask must be set to 0. on rows 1 and jpi
      IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
         bmask( 1 ,:) = 0._wp
         bmask(jpi,:) = 0._wp
      ENDIF
      IF( nperio == 2 ) THEN               ! south symmetric :  bmask must be set to 0. on row 1
         bmask(:, 1 ) = 0._wp
      ENDIF
      !                                    ! north fold : 
      IF( nperio == 3 .OR. nperio == 4 ) THEN   ! T-pt pivot : bmask set to 0. on row jpj and on half jpjglo-1 row
         DO ji = 1, jpi                      
            ii = ji + nimpp - 1
            bmask(ji,jpj-1) = bmask(ji,jpj-1) * tpol(ii)
            bmask(ji,jpj  ) = 0._wp
         END DO
      ENDIF
      IF( nperio == 5 .OR. nperio == 6 ) THEN   ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
         bmask(:,jpj) = 0._wp
      ENDIF
      !
      IF( lk_mpp ) THEN                    ! mpp specificities
         !                                      ! bmask is set to zero on the overlap region
         IF( nbondi /= -1 .AND. nbondi /= 2 )   bmask(  1 :jpreci,:) = 0._wp
         IF( nbondi /=  1 .AND. nbondi /= 2 )   bmask(nlci:jpi   ,:) = 0._wp
         IF( nbondj /= -1 .AND. nbondj /= 2 )   bmask(:,  1 :jprecj) = 0._wp
         IF( nbondj /=  1 .AND. nbondj /= 2 )   bmask(:,nlcj:jpj   ) = 0._wp
         !
         IF( npolj == 3 .OR. npolj == 4 ) THEN  ! north fold : bmask must be set to 0. on rows jpj-1 and jpj
            DO ji = 1, nlci
               ii = ji + nimpp - 1
               bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * tpol(ii)
               bmask(ji,nlcj  ) = 0._wp
            END DO
         ENDIF
         IF( npolj == 5 .OR. npolj == 6 ) THEN  ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
            DO ji = 1, nlci
               bmask(ji,nlcj  ) = 0._wp
            END DO
         ENDIF
      ENDIF


      ! mask for second order calculation of vorticity
      ! ----------------------------------------------
      CALL dom_msk_nsa

      
      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------     
      DO jk = 1, jpk
         zwf(:,:) = fmask(:,:,jk)         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               IF( fmask(ji,jj,jk) == 0._wp ) THEN
                  fmask(ji,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                     &                                           zwf(ji-1,jj), zwf(ji,jj-1)  )  )
               ENDIF
            END DO
         END DO
         DO jj = 2, jpjm1
            IF( fmask(1,jj,jk) == 0._wp ) THEN
               fmask(1  ,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
            ENDIF
            IF( fmask(jpi,jj,jk) == 0._wp ) THEN
               fmask(jpi,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
            ENDIF
         END DO         
         DO ji = 2, jpim1
            IF( fmask(ji,1,jk) == 0._wp ) THEN
               fmask(ji, 1 ,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
            ENDIF
            IF( fmask(ji,jpj,jk) == 0._wp ) THEN
               fmask(ji,jpj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
            ENDIF
         END DO
      END DO
      !
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN   ! ORCA_R2 configuration
         !                                                 ! Increased lateral friction near of some straits
         IF( nn_cla == 0 ) THEN
            !                                ! Gibraltar strait  : partial slip (fmask=0.5)
            ij0 = 101   ;   ij1 = 101
            ii0 = 139   ;   ii1 = 140   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            ij0 = 102   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 140   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
            !
            !                                ! Bab el Mandeb : partial slip (fmask=1)
            ij0 =  87   ;   ij1 =  88
            ii0 = 160   ;   ii1 = 160   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            ij0 =  88   ;   ij1 =  88
            ii0 = 159   ;   ii1 = 159   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
            !
         ENDIF
         !                                ! Danish straits  : strong slip (fmask > 2)
! We keep this as an example but it is instable in this case 
!         ij0 = 115   ;   ij1 = 115
!         ii0 = 145   ;   ii1 = 146   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
!         ij0 = 116   ;   ij1 = 116
!         ii0 = 145   ;   ii1 = 146   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
         !
      ENDIF
      !
      IF( cp_cfg == "orca" .AND. jp_cfg == 1 ) THEN   ! ORCA R1 configuration
         !                                                 ! Increased lateral friction near of some straits
         ! This dirty section will be suppressed by simplification process:
         ! all this will come back in input files
         ! Currently these hard-wired indices relate to configuration with
         ! extend grid (jpjglo=332)
         !
         isrow = 332 - jpjglo
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   orca_r1: increase friction near the following straits : '
         IF(lwp) WRITE(numout,*) '      Gibraltar '
         ii0 = 282           ;   ii1 = 283        ! Gibraltar Strait 
         ij0 = 241 - isrow   ;   ij1 = 241 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Bhosporus '
         ii0 = 314           ;   ii1 = 315        ! Bhosporus Strait 
         ij0 = 248 - isrow   ;   ij1 = 248 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Makassar (Top) '
         ii0 =  48           ;   ii1 =  48        ! Makassar Strait (Top) 
         ij0 = 189 - isrow   ;   ij1 = 190 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  

         IF(lwp) WRITE(numout,*) '      Lombok '
         ii0 =  44           ;   ii1 =  44        ! Lombok Strait 
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Ombai '
         ii0 =  53           ;   ii1 =  53        ! Ombai Strait 
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Timor Passage '
         ii0 =  56           ;   ii1 =  56        ! Timor Passage 
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      West Halmahera '
         ii0 =  58           ;   ii1 =  58        ! West Halmahera Strait 
         ij0 = 181 - isrow   ;   ij1 = 182 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  

         IF(lwp) WRITE(numout,*) '      East Halmahera '
         ii0 =  55           ;   ii1 =  55        ! East Halmahera Strait 
         ij0 = 181 - isrow   ;   ij1 = 182 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
         !
      ENDIF
      !
      CALL lbc_lnk( fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask

      ! CAUTION : The fmask may be further modified in dyn_vor_init ( dynvor.F90 )
            
      IF( nprint == 1 .AND. lwp ) THEN      ! Control print
         imsk(:,:) = INT( tmask_i(:,:) )
         WRITE(numout,*) ' tmask_i : '
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                           1, jpj, 1, 1, numout)
         WRITE (numout,*)
         WRITE (numout,*) ' dommsk: tmask for each level'
         WRITE (numout,*) ' ----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( tmask(:,:,jk) )

            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: vmask for each level'
         WRITE(numout,*) ' -----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( vmask(:,:,jk) )
            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: fmask for each level'
         WRITE(numout,*) ' -----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( fmask(:,:,jk) )
            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout )
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: bmask '
         WRITE(numout,*) ' ---------------'
         WRITE(numout,*)
         imsk(:,:) = INT( bmask(:,:) )
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
            &                              1, jpj, 1, 1, numout )
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, imsk )
      CALL wrk_dealloc( jpi, jpj, zwf  )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_msk')
      !
   END SUBROUTINE dom_msk

   !!----------------------------------------------------------------------
   !!   Default option :                                      Empty routine
   !!----------------------------------------------------------------------
   SUBROUTINE dom_msk_nsa       
   END SUBROUTINE dom_msk_nsa
   
   !!======================================================================
END MODULE dommsk
