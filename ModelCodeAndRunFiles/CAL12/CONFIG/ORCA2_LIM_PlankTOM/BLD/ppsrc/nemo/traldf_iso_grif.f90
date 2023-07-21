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

MODULE traldf_iso_grif
   !!======================================================================
   !!                   ***  MODULE  traldf_iso_grif  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History : 3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)
   !!                !          Griffies operator version adapted from traldf_iso.F90
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'               slope of the lateral diffusive direction
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso_grif  : update the tracer trend with the horizontal component
   !!                       of the Griffies iso-neutral laplacian operator
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE trc_oce         ! share passive tracers/Ocean variables
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfslp          ! iso-neutral slopes
   USE diaptr          ! poleward transport diagnostics
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing


   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_iso_grif   ! routine called by traldf.F90

   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   psix_eiv, psiy_eiv   !: eiv stream function (diag only)
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   ah_wslp2             !: aeiv*w-slope^2
   REAL(wp),         DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   zdkt3d               !: vertical tracer gradient at 2 levels

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c2d' :                 aht: 2D coefficient
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
   !!                   ***  ldfeiv_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaei. the eddy induced velocity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfeiv_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'traldf_c2d' :                           eiv: 2D coefficient
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traldf_iso_grif.F90 5147 2015-03-13 10:01:32Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

  SUBROUTINE tra_ldf_iso_grif( kt, kit000, cdtype, pgu, pgv,              &
       &                                   ptb, pta, kjpt, pahtb0 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_iso_grif  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend for a laplacian tensor (ezxcept the dz[ dz[.] ] term) and
      !!      add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The horizontal component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      1st part :  masked horizontal derivative of T  ( di[ t ] )
      !!      ========    with partial cell update if ln_zps=T.
      !!
      !!      2nd part :  horizontal fluxes of the lateral mixing operator
      !!      ========
      !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
      !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
      !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
      !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !!      3rd part: vertical trends of the lateral mixing operator
      !!      ========  (excluding the vertical flux proportional to dk[t] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ]
      !!      Add this trend to the general trend (ta,sa):
      !!         pta = pta + difft
      !!
      !! ** Action :   Update pta arrays with the before rotated diffusion
      !!----------------------------------------------------------------------
      USE oce     , ONLY:   zftu => ua       , zftv => va            ! (ua,va) used as 3D workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend
      REAL(wp)                             , INTENT(in   ) ::   pahtb0     ! background diffusion coef
      !
      INTEGER  ::  ji, jj, jk,jn   ! dummy loop indices
      INTEGER  ::  ip,jp,kp        ! dummy loop indices
      INTEGER  ::  ierr            ! temporary integer
      REAL(wp) ::  zmsku, zabe1, zcof1, zcoef3   ! local scalars
      REAL(wp) ::  zmskv, zabe2, zcof2, zcoef4   !   -      -
      REAL(wp) ::  zcoef0, zbtr                  !   -      -
      !
      REAL(wp) ::   zslope_skew, zslope_iso, zslope2, zbu, zbv
      REAL(wp) ::   ze1ur, zdxt, ze2vr, ze3wr, zdyt, zdzt
      REAL(wp) ::   zah, zah_slp, zaei_slp
      REAL(wp), POINTER, DIMENSION(:,:  ) :: z2d
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdit, zdjt, ztfw 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zw3d       ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf_iso_grif')
      !
      CALL wrk_alloc( jpi, jpj,      z2d ) 
      CALL wrk_alloc( jpi, jpj, jpk, zdit, zdjt, ztfw  ) 
      !

      IF( kt == kit000 .AND. .NOT.ALLOCATED(ah_wslp2) )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_iso_grif : rotated laplacian diffusion operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ALLOCATE( ah_wslp2(jpi,jpj,jpk) , zdkt3d(jpi,jpj,0:1), STAT=ierr )
         IF( lk_mpp   )   CALL mpp_sum ( ierr )
         IF( ierr > 0 )   CALL ctl_stop('STOP', 'tra_ldf_iso_grif: unable to allocate arrays')
         IF( ln_traldf_gdia ) THEN
            IF (.NOT. ALLOCATED(psix_eiv))THEN
                ALLOCATE( psix_eiv(jpi,jpj,jpk) , psiy_eiv(jpi,jpj,jpk) , STAT=ierr )
                IF( lk_mpp   )   CALL mpp_sum ( ierr )
                IF( ierr > 0 )   CALL ctl_stop('STOP', 'tra_ldf_iso_grif: unable to allocate diagnostics')
            ENDIF
         ENDIF
      ENDIF

      !!----------------------------------------------------------------------
      !!   0 - calculate  ah_wslp2, psix_eiv, psiy_eiv
      !!----------------------------------------------------------------------

      !!gm Future development: consider using Ah defined at T-points and attached to the 4 t-point triads

      ah_wslp2(:,:,:) = 0._wp
      IF( ln_traldf_gdia ) THEN
         psix_eiv(:,:,:) = 0._wp
         psiy_eiv(:,:,:) = 0._wp
      ENDIF

      DO ip = 0, 1
         DO kp = 0, 1
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     ze1ur = 1._wp / e1u(ji,jj)
                     ze3wr = 1._wp / e3w_0(ji+ip,jj,jk+kp)
                     zbu   = 0.25_wp * e1u(ji,jj) * e2u(ji,jj) * e3u_0(ji,jj,jk)
                     zah   = rldf * ahtu(ji,jj)                                  ! fsaht(ji+ip,jj,jk)
                     zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                     ! Subtract s-coordinate slope at t-points to give slope rel to s surfaces
                     ! (do this by *adding* gradient of depth)
                     zslope2 = zslope_skew + ( gdept_0(ji+1,jj,jk) - gdept_0(ji ,jj ,jk) ) * ze1ur * umask(ji,jj,jk+kp)
                     zslope2 = zslope2 *zslope2
                     ah_wslp2(ji+ip,jj,jk+kp) = ah_wslp2(ji+ip,jj,jk+kp)    &
                        &                     + zah * ( zbu * ze3wr / ( e1t(ji+ip,jj) * e2t(ji+ip,jj) ) ) * zslope2
                     IF( ln_traldf_gdia ) THEN
                        zaei_slp = aeiw(ji+ip,jj) * zslope_skew           ! fsaeit(ji+ip,jj,jk)*zslope_skew
                        psix_eiv(ji,jj,jk+kp) = psix_eiv(ji,jj,jk+kp) + 0.25_wp * zaei_slp
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      DO jp = 0, 1
         DO kp = 0, 1
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji=1,jpim1
                     ze2vr = 1._wp / e2v(ji,jj)
                     ze3wr = 1.0_wp / e3w_0(ji,jj+jp,jk+kp)
                     zbv   = 0.25_wp * e1v(ji,jj) * e2v(ji,jj) * e3v_0(ji,jj,jk)
                     zah   = rldf * ahtv(ji,jj)                                  ! fsaht(ji,jj+jp,jk)
                     zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                     ! Subtract s-coordinate slope at t-points to give slope rel to s surfaces
                     !    (do this by *adding* gradient of depth)
                     zslope2 = zslope_skew + ( gdept_0(ji,jj+1,jk) - gdept_0(ji,jj,jk) ) * ze2vr * vmask(ji,jj,jk+kp)
                     zslope2 = zslope2 * zslope2
                     ah_wslp2(ji,jj+jp,jk+kp) = ah_wslp2(ji,jj+jp,jk+kp)   &
                        &                     + zah * ( zbv * ze3wr / ( e1t(ji,jj+jp) * e2t(ji,jj+jp) ) ) * zslope2
                     IF( ln_traldf_gdia ) THEN
                        zaei_slp = aeiw(ji,jj+jp) * zslope_skew           ! fsaeit(ji,jj+jp,jk)*zslope_skew
                        psiy_eiv(ji,jj,jk+kp) = psiy_eiv(ji,jj,jk+kp) + 0.25_wp * zaei_slp
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      IF( iom_use("uoce_eiv") .OR. iom_use("voce_eiv") .OR. iom_use("woce_eiv") )  THEN
         !
         IF( ln_traldf_gdia .AND. cdtype == 'TRA' ) THEN
            CALL wrk_alloc( jpi , jpj , jpk  , zw3d )
            DO jk=1,jpkm1
               zw3d(:,:,jk) = (psix_eiv(:,:,jk+1) - psix_eiv(:,:,jk))/e3u_0(:,:,jk)  ! u_eiv = -dpsix/dz
            END DO
            zw3d(:,:,jpk) = 0._wp
            CALL iom_put( "uoce_eiv", zw3d )    ! i-eiv current

            DO jk=1,jpk-1
               zw3d(:,:,jk) = (psiy_eiv(:,:,jk+1) - psiy_eiv(:,:,jk))/e3v_0(:,:,jk)  ! v_eiv = -dpsiy/dz
            END DO
            zw3d(:,:,jpk) = 0._wp
            CALL iom_put( "voce_eiv", zw3d )    ! j-eiv current

            DO jk=1,jpk-1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1  ! vector opt.
                     zw3d(ji,jj,jk) = (psiy_eiv(ji,jj,jk) - psiy_eiv(ji,jj-1,jk))/e2t(ji,jj) + &
                          &    (psix_eiv(ji,jj,jk) - psix_eiv(ji-1,jj,jk))/e1t(ji,jj) ! w_eiv = dpsiy/dy + dpsiy/dx
                  END DO
               END DO
            END DO
            zw3d(:,:,jpk) = 0._wp
            CALL iom_put( "woce_eiv", zw3d )    ! vert. eiv current
            CALL wrk_dealloc( jpi , jpj , jpk  , zw3d )
         ENDIF
         !
      ENDIF
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         ! Zero fluxes for each tracer
         ztfw(:,:,:) = 0._wp
         zftu(:,:,:) = 0._wp
         zftv(:,:,:) = 0._wp
         !
         DO jk = 1, jpkm1                          !==  before lateral T & S gradients at T-level jk  ==!
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zdit(ji,jj,jk) = ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                  zdjt(ji,jj,jk) = ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_zps.and.l_grad_zps ) THEN              ! partial steps: correction at the last level
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zdit(ji,jj,mbku(ji,jj)) = pgu(ji,jj,jn)
                  zdjt(ji,jj,mbkv(ji,jj)) = pgv(ji,jj,jn)
               END DO
            END DO
         ENDIF

         !!----------------------------------------------------------------------
         !!   II - horizontal trend  (full)
         !!----------------------------------------------------------------------
         !
         DO jk = 1, jpkm1
            !
            !                    !==  Vertical tracer gradient at level jk and jk+1
            zdkt3d(:,:,1) = ( ptb(:,:,jk,jn) - ptb(:,:,jk+1,jn) ) * tmask(:,:,jk+1)
            !
            !                    ! surface boundary condition: zdkt3d(jk=0)=zdkt3d(jk=1)
            IF( jk == 1 ) THEN   ;   zdkt3d(:,:,0) = zdkt3d(:,:,1)
            ELSE                 ;   zdkt3d(:,:,0) = ( ptb(:,:,jk-1,jn) - ptb(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF


            IF (ln_botmix_grif) THEN
               DO ip = 0, 1              !==  Horizontal & vertical fluxes
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, jpim1
                           ze1ur = 1._wp / e1u(ji,jj)
                           zdxt  = zdit(ji,jj,jk) * ze1ur
                           ze3wr = 1._wp / e3w_0(ji+ip,jj,jk+kp)
                           zdzt  = zdkt3d(ji+ip,jj,kp) * ze3wr
                           zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                           zslope_iso  = triadi(ji+ip,jj,jk,1-ip,kp)

                           zbu = 0.25_wp * e1u(ji,jj) * e2u(ji,jj) * e3u_0(ji,jj,jk)
                           ! ln_botmix_grif is .T. don't mask zah for bottom half cells
                           zah = rldf * ahtu(ji,jj)   !*umask(ji,jj,jk+kp)         !fsaht(ji+ip,jj,jk)           ===>>  ????
                           zah_slp  = zah * zslope_iso
                           zaei_slp = aeiw(ji+ip,jj) * zslope_skew    !fsaeit(ji+ip,jj,jk)*zslope_skew
                           zftu(ji,jj,jk) = zftu(ji,jj,jk) - ( zah * zdxt + (zah_slp - zaei_slp) * zdzt ) * zbu * ze1ur
                           ztfw(ji+ip,jj,jk+kp) = ztfw(ji+ip,jj,jk+kp) - (zah_slp + zaei_slp) * zdxt * zbu * ze3wr
                        END DO
                     END DO
                  END DO
               END DO

               DO jp = 0, 1
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, jpim1
                           ze2vr = 1._wp / e2v(ji,jj)
                           zdyt  = zdjt(ji,jj,jk) * ze2vr
                           ze3wr = 1._wp / e3w_0(ji,jj+jp,jk+kp)
                           zdzt  = zdkt3d(ji,jj+jp,kp) * ze3wr
                           zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                           zslope_iso  = triadj(ji,jj+jp,jk,1-jp,kp)
                           zbv = 0.25_wp * e1v(ji,jj) * e2v(ji,jj) * e3v_0(ji,jj,jk)
                           ! ln_botmix_grif is .T. don't mask zah for bottom half cells
                           zah = rldf * ahtv(ji,jj)        !*vmask(ji,jj,jk+kp)  ! fsaht(ji,jj+jp,jk)
                           zah_slp = zah * zslope_iso
                           zaei_slp = aeiw(ji,jj+jp) * zslope_skew        ! fsaeit(ji,jj+jp,jk)*zslope_skew
                           zftv(ji,jj,jk) = zftv(ji,jj,jk) - ( zah * zdyt + (zah_slp - zaei_slp) * zdzt ) * zbv * ze2vr
                           ztfw(ji,jj+jp,jk+kp) = ztfw(ji,jj+jp,jk+kp) - (zah_slp + zaei_slp) * zdyt * zbv * ze3wr
                        END DO
                     END DO
                  END DO
               END DO
            ELSE
               DO ip = 0, 1              !==  Horizontal & vertical fluxes
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, jpim1
                           ze1ur = 1._wp / e1u(ji,jj)
                           zdxt  = zdit(ji,jj,jk) * ze1ur
                           ze3wr = 1._wp / e3w_0(ji+ip,jj,jk+kp)
                           zdzt  = zdkt3d(ji+ip,jj,kp) * ze3wr
                           zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                           zslope_iso  = triadi(ji+ip,jj,jk,1-ip,kp)

                           zbu = 0.25_wp * e1u(ji,jj) * e2u(ji,jj) * e3u_0(ji,jj,jk)
                           ! ln_botmix_grif is .F. mask zah for bottom half cells
                           zah = rldf * ahtu(ji,jj) * umask(ji,jj,jk+kp)         ! fsaht(ji+ip,jj,jk)   ===>>  ????
                           zah_slp  = zah * zslope_iso
                           zaei_slp = aeiw(ji+ip,jj) * zslope_skew        ! fsaeit(ji+ip,jj,jk)*zslope_skew
                           zftu(ji,jj,jk) = zftu(ji,jj,jk) - ( zah * zdxt + (zah_slp - zaei_slp) * zdzt ) * zbu * ze1ur
                           ztfw(ji+ip,jj,jk+kp) = ztfw(ji+ip,jj,jk+kp) - (zah_slp + zaei_slp) * zdxt * zbu * ze3wr
                        END DO
                     END DO
                  END DO
               END DO

               DO jp = 0, 1
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, jpim1
                           ze2vr = 1._wp / e2v(ji,jj)
                           zdyt  = zdjt(ji,jj,jk) * ze2vr
                           ze3wr = 1._wp / e3w_0(ji,jj+jp,jk+kp)
                           zdzt  = zdkt3d(ji,jj+jp,kp) * ze3wr
                           zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                           zslope_iso  = triadj(ji,jj+jp,jk,1-jp,kp)
                           zbv = 0.25_wp * e1v(ji,jj) * e2v(ji,jj) * e3v_0(ji,jj,jk)
                           ! ln_botmix_grif is .F. mask zah for bottom half cells
                           zah = rldf * ahtv(ji,jj) * vmask(ji,jj,jk+kp)         ! fsaht(ji,jj+jp,jk)
                           zah_slp = zah * zslope_iso
                           zaei_slp = aeiw(ji,jj+jp) * zslope_skew        ! fsaeit(ji,jj+jp,jk)*zslope_skew
                           zftv(ji,jj,jk) = zftv(ji,jj,jk) - ( zah * zdyt + (zah_slp - zaei_slp) * zdzt ) * zbv * ze2vr
                           ztfw(ji,jj+jp,jk+kp) = ztfw(ji,jj+jp,jk+kp) - (zah_slp + zaei_slp) * zdyt * zbv * ze3wr
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            !                          !==  divergence and add to the general trend  ==!
            DO jj = 2 , jpjm1
               DO ji = 2, jpim1  ! vector opt.
                  zbtr = 1._wp / ( e1t(ji,jj) * e2t(ji,jj) * e3t_0(ji,jj,jk) )
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zbtr * (   zftu(ji-1,jj,jk) - zftu(ji,jj,jk)   &
                     &                                           + zftv(ji,jj-1,jk) - zftv(ji,jj,jk)   )
               END DO
            END DO
            !
         END DO
         !
         DO jk = 1, jpkm1              !== Divergence of vertical fluxes added to the general tracer trend
            DO jj = 2, jpjm1
               DO ji = 2, jpim1  ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (  ztfw(ji,jj,jk+1) - ztfw(ji,jj,jk)  )   &
                     &                                / ( e1t(ji,jj) * e2t(ji,jj) * e3t_0(ji,jj,jk) )
               END DO
            END DO
         END DO
         !
         !                             ! "Poleward" diffusive heat or salt transports (T-S case only)
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN
            IF( jn == jp_tem)   htr_ldf(:) = ptr_sj( zftv(:,:,:) )        ! 3.3  names
            IF( jn == jp_sal)   str_ldf(:) = ptr_sj( zftv(:,:,:) )
         ENDIF

         IF( iom_use("udiff_heattr") .OR. iom_use("vdiff_heattr") ) THEN
           !
           IF( cdtype == 'TRA' .AND. jn == jp_tem  ) THEN
               z2d(:,:) = 0._wp 
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        z2d(ji,jj) = z2d(ji,jj) + zftu(ji,jj,jk) 
                     END DO
                  END DO
               END DO
               z2d(:,:) = rau0_rcp * z2d(:,:) 
               CALL lbc_lnk( z2d, 'U', -1. )
               CALL iom_put( "udiff_heattr", z2d )                  ! heat transport in i-direction
               !
               z2d(:,:) = 0._wp 
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        z2d(ji,jj) = z2d(ji,jj) + zftv(ji,jj,jk) 
                     END DO
                  END DO
               END DO
               z2d(:,:) = rau0_rcp * z2d(:,:)     
               CALL lbc_lnk( z2d, 'V', -1. )
               CALL iom_put( "vdiff_heattr", z2d )                  !  heat transport in i-direction
            END IF
            !
         ENDIF
         !
      END DO
      !
      CALL wrk_dealloc( jpi, jpj,      z2d ) 
      CALL wrk_dealloc( jpi, jpj, jpk, zdit, zdjt, ztfw  ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf_iso_grif')
      !
  END SUBROUTINE tra_ldf_iso_grif


   !!==============================================================================
END MODULE traldf_iso_grif
