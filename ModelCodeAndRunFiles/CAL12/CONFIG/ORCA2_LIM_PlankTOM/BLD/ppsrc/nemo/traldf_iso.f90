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

MODULE traldf_iso
   !!======================================================================
   !!                   ***  MODULE  traldf_iso  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History :  OPA  !  1994-08  (G. Madec, M. Imbard)
   !!            8.0  !  1997-05  (G. Madec)  split into traldf and trazdf
   !!            NEMO !  2002-08  (G. Madec)  Free form, F90
   !!            1.0  !  2005-11  (G. Madec)  merge traldf and trazdf :-)
   !!            3.3  !  2010-09  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'               slope of the lateral diffusive direction
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso  : update the tracer trend with the horizontal 
   !!                  component of a iso-neutral laplacian operator
   !!                  and with the vertical part of
   !!                  the isopycnal or geopotential s-coord. operator 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE trc_oce         ! share passive tracers/Ocean variables
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfslp          ! iso-neutral slopes
   USE diaptr          ! poleward transport diagnostics
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_iso   ! routine called by step.F90

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: traldf_iso.F90 5149 2015-03-13 21:44:31Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf_iso( kt, kit000, cdtype, pgu, pgv,              &
      &                                pgui, pgvi,                    &
      &                                ptb, pta, kjpt, pahtb0 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_iso  ***
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
      USE oce     , ONLY:   zftu => ua       , zftv  => va         ! (ua,va) used as workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu , pgv    ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgui, pgvi   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      REAL(wp)                             , INTENT(in   ) ::   pahtb0     ! background diffusion coef
      !
      INTEGER  ::  ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::  ikt
      REAL(wp) ::  zmsku, zabe1, zcof1, zcoef3   ! local scalars
      REAL(wp) ::  zmskv, zabe2, zcof2, zcoef4   !   -      -
      REAL(wp) ::  zcoef0, zbtr, ztra            !   -      -
      REAL(wp), POINTER, DIMENSION(:,:  ) ::  z2d
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zdkt, zdk1t, zdit, zdjt, ztfw 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf_iso')
      !
      CALL wrk_alloc( jpi, jpj,      z2d ) 
      CALL wrk_alloc( jpi, jpj, jpk, zdit, zdjt, ztfw, zdkt, zdk1t ) 
      !

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_iso : rotated laplacian diffusion operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         !                                               
         !!----------------------------------------------------------------------
         !!   I - masked horizontal derivative 
         !!----------------------------------------------------------------------
         !!bug ajout.... why?   ( 1,jpj,:) and (jpi,1,:) should be sufficient....
         zdit (1,:,:) = 0.e0     ;     zdit (jpi,:,:) = 0.e0
         zdjt (1,:,:) = 0.e0     ;     zdjt (jpi,:,:) = 0.e0
         !!end

         ! Horizontal tracer gradient 
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zdit(ji,jj,jk) = ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                  zdjt(ji,jj,jk) = ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO

         ! partial cell correction
         IF( ln_zps ) THEN      ! partial steps correction at the last ocean level 
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
! IF useless if zpshde defines pgu everywhere
                  zdit(ji,jj,mbku(ji,jj)) = pgu(ji,jj,jn)          
                  zdjt(ji,jj,mbkv(ji,jj)) = pgv(ji,jj,jn)
               END DO
            END DO
         ENDIF
         IF( ln_zps .AND. ln_isfcav ) THEN      ! partial steps correction at the first wet level beneath a cavity
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  IF (miku(ji,jj) > 1) zdit(ji,jj,miku(ji,jj)) = pgui(ji,jj,jn)          
                  IF (mikv(ji,jj) > 1) zdjt(ji,jj,mikv(ji,jj)) = pgvi(ji,jj,jn)     
               END DO
            END DO
         END IF

         !!----------------------------------------------------------------------
         !!   II - horizontal trend  (full)
         !!----------------------------------------------------------------------
!!!!!!!!!!CDIR PARALLEL DO PRIVATE( zdk1t ) 
            ! 1. Vertical tracer gradient at level jk and jk+1
            ! ------------------------------------------------
         ! 
         ! interior value 
         DO jk = 2, jpkm1               
            DO jj = 1, jpj
               DO ji = 1, jpi   ! vector opt.
                  zdk1t(ji,jj,jk) = ( ptb(ji,jj,jk,jn  ) - ptb(ji,jj,jk+1,jn) ) * wmask(ji,jj,jk+1)
                  !
                  zdkt(ji,jj,jk)  = ( ptb(ji,jj,jk-1,jn) - ptb(ji,jj,jk,jn  ) ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
         zdk1t(:,:,1) = ( ptb(:,:,1,jn  ) - ptb(:,:,2,jn) ) * wmask(:,:,2)
         zdkt (:,:,1) = zdk1t(:,:,1)
         IF ( ln_isfcav ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi   ! vector opt.
                  ikt = mikt(ji,jj) ! surface level
                  zdk1t(ji,jj,ikt) = ( ptb(ji,jj,ikt,jn  ) - ptb(ji,jj,ikt+1,jn) ) * wmask(ji,jj,ikt+1)
                  zdkt (ji,jj,ikt) = zdk1t(ji,jj,ikt)
               END DO
            END DO
         END IF

         ! 2. Horizontal fluxes
         ! --------------------   
         DO jk = 1, jpkm1
            DO jj = 1 , jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zabe1 = ( rldf * ahtu(ji,jj) + pahtb0 ) * re2u_e1u(ji,jj) * e3u_0(ji,jj,jk)
                  zabe2 = ( rldf * ahtv(ji,jj) + pahtb0 ) * re1v_e2v(ji,jj) * e3v_0(ji,jj,jk)
                  !
                  zmsku = 1. / MAX(  tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
                     &             + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
                  !
                  zmskv = 1. / MAX(  tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
                     &             + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
                  !
                  zcof1 = - rldf * ahtu(ji,jj) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
                  zcof2 = - rldf * ahtv(ji,jj) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv
                  !
                  zftu(ji,jj,jk ) = (  zabe1 * zdit(ji,jj,jk)   &
                     &              + zcof1 * (  zdkt (ji+1,jj,jk) + zdk1t(ji,jj,jk)      &
                     &                         + zdk1t(ji+1,jj,jk) + zdkt (ji,jj,jk)  )  ) * umask(ji,jj,jk)
                  zftv(ji,jj,jk) = (  zabe2 * zdjt(ji,jj,jk)   &
                     &              + zcof2 * (  zdkt (ji,jj+1,jk) + zdk1t(ji,jj,jk)      &
                     &                         + zdk1t(ji,jj+1,jk) + zdkt (ji,jj,jk)  )  ) * vmask(ji,jj,jk)                  
               END DO
            END DO

            ! II.4 Second derivative (divergence) and add to the general trend
            ! ----------------------------------------------------------------
            DO jj = 2 , jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr = 1.0 / ( e12t(ji,jj) * e3t_0(ji,jj,jk) )
                  ztra = zbtr * ( zftu(ji,jj,jk) - zftu(ji-1,jj,jk) + zftv(ji,jj,jk) - zftv(ji,jj-1,jk)  )
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
            !                                          ! ===============
         END DO                                        !   End of slab  
         !                                             ! ===============
         !
         ! "Poleward" diffusive heat or salt transports (T-S case only)
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN
            ! note sign is reversed to give down-gradient diffusive transports (#1043)
            IF( jn == jp_tem)   htr_ldf(:) = ptr_sj( -zftv(:,:,:) )
            IF( jn == jp_sal)   str_ldf(:) = ptr_sj( -zftv(:,:,:) )
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
               z2d(:,:) = - rau0_rcp * z2d(:,:)     ! note sign is reversed to give down-gradient diffusive transports (#1043)
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
               z2d(:,:) = - rau0_rcp * z2d(:,:)     ! note sign is reversed to give down-gradient diffusive transports (#1043)
               CALL lbc_lnk( z2d, 'V', -1. )
               CALL iom_put( "vdiff_heattr", z2d )                  !  heat transport in i-direction
            END IF
            !
         ENDIF

         !!----------------------------------------------------------------------
         !!   III - vertical trend of T & S (extra diagonal terms only)
         !!----------------------------------------------------------------------
         
         ! Local constant initialization
         ! -----------------------------
         ztfw(1,:,:) = 0.e0     ;     ztfw(jpi,:,:) = 0.e0
         
         ! Vertical fluxes
         ! ---------------
         
         ! Surface and bottom vertical fluxes set to zero
         ztfw(:,:, 1 ) = 0.e0      ;      ztfw(:,:,jpk) = 0.e0
         
         ! interior (2=<jk=<jpk-1)
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zcoef0 = - rldf * ahtw(ji,jj) * wmask(ji,jj,jk)
                  !
                  zmsku = 1./MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
                     &            + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.  )
                  zmskv = 1./MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
                     &            + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.  )
                  !
                  zcoef3 = zcoef0 * e2t(ji,jj) * zmsku * wslpi (ji,jj,jk)
                  zcoef4 = zcoef0 * e1t(ji,jj) * zmskv * wslpj (ji,jj,jk)
                  !
                  ztfw(ji,jj,jk) = zcoef3 * (   zdit(ji  ,jj  ,jk-1) + zdit(ji-1,jj  ,jk)      &
                     &                        + zdit(ji-1,jj  ,jk-1) + zdit(ji  ,jj  ,jk)  )   &
                     &           + zcoef4 * (   zdjt(ji  ,jj  ,jk-1) + zdjt(ji  ,jj-1,jk)      &
                     &                        + zdjt(ji  ,jj-1,jk-1) + zdjt(ji  ,jj  ,jk)  )
               END DO
            END DO
         END DO
         
         
         ! I.5 Divergence of vertical fluxes added to the general tracer trend
         ! -------------------------------------------------------------------
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr = 1.0 / ( e12t(ji,jj) * e3t_0(ji,jj,jk) )
                  ztra = (  ztfw(ji,jj,jk) - ztfw(ji,jj,jk+1)  ) * zbtr
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO
         !
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, z2d ) 
      CALL wrk_dealloc( jpi, jpj, jpk, zdit, zdjt, ztfw, zdkt, zdk1t ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf_iso')
      !
   END SUBROUTINE tra_ldf_iso


   !!==============================================================================
END MODULE traldf_iso
