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

MODULE traadv_cen2
   !!======================================================================
   !!                     ***  MODULE  traadv_cen2  ***
   !! Ocean  tracers:  horizontal & vertical advective trend
   !!======================================================================
   !! History :  OPA  ! 2001-08  (G. Madec, E. Durand) v8.2 trahad+trazad=traadv 
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!             -   ! 2004-08  (C. Talandier) New trends organization
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.0  ! 2006-04  (R. Benshila, G. Madec) Step reorganization
   !!             -   ! 2006-07  (G. madec)  add ups_orca_set routine
   !!            3.2  ! 2009-07  (G. Madec) add avmb, avtb in restart for cen2 advection
   !!            3.3  ! 2010-05  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_cen2 : update the tracer trend with the advection trends using a 2nd order centered scheme
   !!   ups_orca_set : allow mixed upstream/centered scheme in specific area (set for orca 2 and 4 only)
   !!----------------------------------------------------------------------
   USE oce, ONLY: tsn  ! now ocean temperature and salinity
   USE dom_oce         ! ocean space and time domain
   USE eosbn2          ! equation of state
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE closea          ! closed sea
   USE sbcrnf          ! river runoffs
   USE in_out_manager  ! I/O manager
   USE iom             ! IOM library
   USE diaptr          ! poleward transport diagnostics
   USE zdf_oce         ! ocean vertical physics
   USE trc_oce         ! share passive tracers/Ocean variables
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE phycst

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_cen2   ! routine called by traadv.F90

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: upsmsk !: mixed upstream/centered scheme near some straits 
   !                                                             !  and in closed seas (orca 2 and 4 configurations)
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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traadv_cen2.F90 5540 2015-07-02 15:11:23Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_cen2( kt, kit000, cdtype, pun, pvn, pwn,     &
      &                                         ptb, ptn, pta, kjpt   ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_cen2  ***
      !!                 
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The advection is evaluated by a second order centered
      !!      scheme using now fields (leap-frog scheme). In specific areas
      !!      (vicinity of major river mouths, some straits, or where tn is
      !!      approaching the freezing point) it is mixed with an upstream
      !!      scheme for stability reasons.
      !!         Part 0 : compute the upstream / centered flag
      !!                  (3D array, zind, defined at T-point (0<zind<1))
      !!         Part I : horizontal advection
      !!       * centered flux:
      !!               zcenu = e2u*e3u  un  mi(ptn)
      !!               zcenv = e1v*e3v  vn  mj(ptn)
      !!       * upstream flux:
      !!               zupsu = e2u*e3u  un  (ptb(i) or ptb(i-1) ) [un>0 or <0]
      !!               zupsv = e1v*e3v  vn  (ptb(j) or ptb(j-1) ) [vn>0 or <0]
      !!       * mixed upstream / centered horizontal advection scheme
      !!               zcofi = max(zind(i+1), zind(i))
      !!               zcofj = max(zind(j+1), zind(j))
      !!               zwx = zcofi * zupsu + (1-zcofi) * zcenu
      !!               zwy = zcofj * zupsv + (1-zcofj) * zcenv
      !!       * horizontal advective trend (divergence of the fluxes)
      !!               ztra = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
      !!       * Add this trend now to the general trend of tracer (ta,sa):
      !!               pta = pta + ztra
      !!       * trend diagnostic (l_trdtra=T or l_trctra=T): the trend is
      !!      saved for diagnostics. The trends saved is expressed as
      !!      Uh.gradh(T), i.e.  save trend = ztra + ptn divn
      !!
      !!         Part II : vertical advection
      !!      For temperature (idem for salinity) the advective trend is com-
      !!      puted as follows :
      !!            ztra = 1/e3t dk+1[ zwz ]
      !!      where the vertical advective flux, zwz, is given by :
      !!            zwz = zcofk * zupst + (1-zcofk) * zcent
      !!      with
      !!        zupsv = upstream flux = wn * (ptb(k) or ptb(k-1) ) [wn>0 or <0]
      !!        zcenu = centered flux = wn * mk(tn)
      !!         The surface boundary condition is :
      !!      variable volume (lk_vvl = T) : zero advective flux
      !!      lin. free-surf  (lk_vvl = F) : wn(:,:,1) * ptn(:,:,1)
      !!         Add this trend now to the general trend of tracer (ta,sa):
      !!             pta = pta + ztra
      !!         Trend diagnostic (l_trdtra=T or l_trctra=T): the trend is
      !!      saved for diagnostics. The trends saved is expressed as :
      !!             save trend =  w.gradz(T) = ztra - ptn divn.
      !!
      !! ** Action :  - update pta  with the now advective tracer trends
      !!              - save trends if needed
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn, ikt   ! dummy loop indices
      INTEGER  ::   ierr                 ! local integer
      REAL(wp) ::   zbtr, ztra                            ! local scalars
      REAL(wp) ::   zfp_ui, zfp_vj, zfp_w, zcofi          !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj, zfm_w, zcofj, zcofk   !   -      -
      REAL(wp) ::   zupsut, zcenut, zupst                 !   -      -
      REAL(wp) ::   zupsvt, zcenvt, zcent, zice           !   -      -
      REAL(wp), POINTER, DIMENSION(:,:)   :: zfzp, zpres   ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwx, zwy     ! 3D     -
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwz, zind    !  -     -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_cen2')
      !
      CALL wrk_alloc( jpi, jpj, zpres, zfzp )
      CALL wrk_alloc( jpi, jpj, jpk, zwx, zwy, zwz, zind )
      !

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_cen2 : 2nd order centered advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~ '
         IF(lwp) WRITE(numout,*)
         !
         IF( .NOT. ALLOCATED( upsmsk ) )  THEN
             ALLOCATE( upsmsk(jpi,jpj), STAT=ierr )
             IF( ierr /= 0 )   CALL ctl_stop('STOP', 'tra_adv_cen2: unable to allocate array')
         ENDIF

         !
         upsmsk(:,:) = 0._wp                             ! not upstream by default
         ! 
         IF( cp_cfg == "orca" )   CALL ups_orca_set      ! set mixed Upstream/centered scheme near some straits
         !                                               ! and in closed seas (orca2 and orca4 only)
         IF( jp_cfg == 2 .AND. .NOT. ln_rstart ) THEN    ! Increase the background in the surface layers
            avmb(1) = 10.  * avmb(1)      ;      avtb(1) = 10.  * avtb(1)
            avmb(2) = 10.  * avmb(2)      ;      avtb(2) = 10.  * avtb(2)
            avmb(3) =  5.  * avmb(3)      ;      avtb(3) =  5.  * avtb(3)
            avmb(4) =  2.5 * avmb(4)      ;      avtb(4) =  2.5 * avtb(4)
         ENDIF
      ENDIF
      !
      ! Upstream / centered scheme indicator
      ! ------------------------------------
!!gm  not strickly exact : the freezing point should be computed at each ocean levels...
!!gm  not a big deal since cen2 is no more used in global ice-ocean simulations
!!ch  changes for ice shelf to retain standard behaviour elsewhere, even if not optimal 
      DO jj = 1, jpj 
         DO ji = 1, jpi 
            ikt = mikt(ji,jj) 
            IF (ikt > 1 ) THEN 
               zpres(ji,jj) = grav * rau0 * gdept_0(ji,jj,ikt) * 1.e-04  
            ELSE 
               zpres(ji,jj) = 0.0 
            ENDIF  
         END DO 
      END DO 
      CALL eos_fzp( tsn(:,:,1,jp_sal), zfzp(:,:), zpres(:,:) )
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               !                                        ! below ice covered area (if tn < "freezing"+0.1 )
               IF( tsn(ji,jj,jk,jp_tem) <= zfzp(ji,jj) + 0.1 ) THEN   ;   zice = 1._wp
               ELSE                                                   ;   zice = 0._wp
               ENDIF
               zind(ji,jj,jk) = MAX (   &
                  rnfmsk(ji,jj) * rnfmsk_z(jk),      &  ! near runoff mouths (& closed sea outflows)
                  upsmsk(ji,jj)               ,      &  ! some of some straits
                  zice                               &  ! below ice covered area (if tn < "freezing"+0.1 )
                  &                  ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      DO jn = 1, kjpt
         !
         ! I. Horizontal advection
         !    ====================
         !
         DO jk = 1, jpkm1
            !                        ! Second order centered tracer flux at u- and v-points
            DO jj = 1, jpjm1
               !
               DO ji = 1, jpim1   ! vector opt.
                  ! upstream indicator
                  zcofi = MAX( zind(ji+1,jj,jk), zind(ji,jj,jk) )
                  zcofj = MAX( zind(ji,jj+1,jk), zind(ji,jj,jk) )
                  !
                  ! upstream scheme
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  zupsut = zfp_ui * ptb(ji,jj,jk,jn) + zfm_ui * ptb(ji+1,jj  ,jk,jn)
                  zupsvt = zfp_vj * ptb(ji,jj,jk,jn) + zfm_vj * ptb(ji  ,jj+1,jk,jn)
                  ! centered scheme
                  zcenut = pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn) )
                  zcenvt = pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn) )
                  ! mixed centered / upstream scheme
                  zwx(ji,jj,jk) = 0.5 * ( zcofi * zupsut + (1.-zcofi) * zcenut )
                  zwy(ji,jj,jk) = 0.5 * ( zcofj * zupsvt + (1.-zcofj) * zcenvt )
               END DO
            END DO
         END DO

         ! II. Vertical advection
         !     ==================
         !
         !                                                ! Vertical advective fluxes
         zwz(:,:,jpk) = 0.e0                                   ! Bottom  value : flux set to zero
         !                                                     ! Surface value : 
         IF( lk_vvl ) THEN   ;   zwz(:,:, 1 ) = 0.e0                         ! volume variable
         ELSE
            DO jj = 1, jpj   ! vector opt.
               DO ji = 1, jpi   ! vector opt.
                  ikt = mikt(ji,jj)                
                  zwz(ji,jj,ikt ) = pwn(ji,jj,ikt) * ptn(ji,jj,ikt,jn)   ! linear free surface 
                  zwz(ji,jj,1:ikt-1) = 0.e0
               END DO
            END DO
         ENDIF
         !
         DO jk = 2, jpk              ! Second order centered tracer flux at w-point
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ! upstream indicator
                  zcofk = MAX( zind(ji,jj,jk-1), zind(ji,jj,jk) ) 
                  ! mixed centered / upstream scheme
                  zfp_w = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                  zfm_w = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                  zupst = zfp_w * ptb(ji,jj,jk,jn) + zfm_w * ptb(ji,jj,jk-1,jn)
                  ! centered scheme
                  zcent = pwn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) )
                  ! mixed centered / upstream scheme
                  zwz(ji,jj,jk) = 0.5 * ( zcofk * zupst + (1.-zcofk) * zcent )
               END DO
            END DO
         END DO

         ! II. Divergence of advective fluxes
         ! ----------------------------------
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) *  e3t_0(ji,jj,jk) )
                  ! advective trends
                  ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &                + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
                  &                + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1)  )
                  ! advective trends added to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO

         !                                 ! trend diagnostics
         IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR.    &
            &( cdtype == 'TRC' .AND. l_trdtrc ) )   THEN
            CALL trd_tra( kt, cdtype, jn, jptra_xad, zwx, pun, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, zwy, pvn, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_zad, zwz, pwn, ptn(:,:,:,jn) )
         END IF
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
           IF( jn == jp_tem )   htr_adv(:) = ptr_sj( zwy(:,:,:) )
           IF( jn == jp_sal )   str_adv(:) = ptr_sj( zwy(:,:,:) )
         ENDIF
         !
      END DO

      ! ---------------------------  required in restart file to ensure restartability)
      ! avmb, avtb will be read in zdfini in restart case as they are used in zdftke, kpp etc...
      IF( lrst_oce .AND. cdtype == 'TRA' ) THEN
         CALL iom_rstput( kt, nitrst, numrow, 'avmb', avmb )
         CALL iom_rstput( kt, nitrst, numrow, 'avtb', avtb )
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zpres, zfzp )
      CALL wrk_dealloc( jpi, jpj, jpk, zwx, zwy, zwz, zind )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_cen2')
      !
   END SUBROUTINE tra_adv_cen2
   
   
   SUBROUTINE ups_orca_set
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ups_orca_set  ***
      !!       
      !! ** Purpose :   add a portion of upstream scheme in area where the
      !!                centered scheme generates too strong overshoot
      !!
      !! ** Method  :   orca (R4 and R2) confiiguration setting. Set upsmsk
      !!                array to nozero value in some straith. 
      !!
      !! ** Action : - upsmsk set to 1 at some strait, 0 elsewhere for orca
      !!----------------------------------------------------------------------
      INTEGER  ::   ii0, ii1, ij0, ij1      ! temporary integers
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ups_orca_set')
      !
      ! mixed upstream/centered scheme near river mouths
      ! ------------------------------------------------
      SELECT CASE ( jp_cfg )
      !                                        ! =======================
      CASE ( 4 )                               !  ORCA_R4 configuration 
         !                                     ! =======================
         !                                          ! Gibraltar Strait
         ii0 =  70   ;   ii1 =  71
         ij0 =  52   ;   ij1 =  53   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         !
         !                                     ! =======================
      CASE ( 2 )                               !  ORCA_R2 configuration 
         !                                     ! =======================
         !                                          ! Gibraltar Strait
         ij0 = 102   ;   ij1 = 102
         ii0 = 138   ;   ii1 = 138   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.20
         ii0 = 139   ;   ii1 = 139   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40
         ii0 = 140   ;   ii1 = 140   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         ij0 = 101   ;   ij1 = 102
         ii0 = 141   ;   ii1 = 141   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         !                                          ! Bab el Mandeb Strait
         ij0 =  87   ;   ij1 =  88
         ii0 = 164   ;   ii1 = 164   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.10
         ij0 =  88   ;   ij1 =  88
         ii0 = 163   ;   ii1 = 163   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
         ii0 = 162   ;   ii1 = 162   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40
         ii0 = 160   ;   ii1 = 161   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         ij0 =  89   ;   ij1 =  89
         ii0 = 158   ;   ii1 = 160   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
         ij0 =  90   ;   ij1 =  90
         ii0 = 160   ;   ii1 = 160   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
         !                                          ! Sound Strait
         ij0 = 116   ;   ij1 = 116
         ii0 = 144   ;   ii1 = 144   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
         ii0 = 145   ;   ii1 = 147   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         ii0 = 148   ;   ii1 = 148   ;   upsmsk( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
         !
      END SELECT 
      
      ! mixed upstream/centered scheme over closed seas
      ! -----------------------------------------------
      CALL clo_ups( upsmsk(:,:) )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ups_orca_set')
      !
   END SUBROUTINE ups_orca_set

   !!======================================================================
END MODULE traadv_cen2
