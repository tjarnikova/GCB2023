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

MODULE traldf_bilap
   !!==============================================================================
   !!                   ***  MODULE  traldf_bilap  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
   !! History :  OPA  !  1991-11  (G. Madec)  Original code
   !!                 !  1993-03  (M. Guyon)  symetrical conditions
   !!                 !  1995-11  (G. Madec)  suppress volumetric scale factors
   !!                 !  1996-01  (G. Madec)  statement function for e3
   !!                 !  1996-01  (M. Imbard)  mpp exchange
   !!                 !  1997-07  (G. Madec)  optimization, and ahtt
   !!            8.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!   NEMO     1.0  !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-11  (G. Madec)  zps or sco as default option
   !!            3.3  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_ldf_bilap : update the tracer trend with the horizontal diffusion
   !!                   using a iso-level biharmonic operator
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE in_out_manager  ! I/O manager
   USE ldfslp          ! iso-neutral slopes 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE diaptr          ! poleward transport diagnostics
   USE trc_oce         ! share passive tracers/Ocean variables
   USE lib_mpp         ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_bilap   ! routine called by step.F90

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
   !! $Id: traldf_bilap.F90 5147 2015-03-13 10:01:32Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE tra_ldf_bilap( kt, kit000, cdtype, pgu, pgv,            &
      &                                          pgui, pgvi,          &
      &                                  ptb, pta, kjpt )  
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   4th order diffusive operator along model level surfaces 
      !!      evaluated using before fields (forward time scheme). The hor.
      !!      diffusive trends  is given by:
      !!      Laplacian of tb:
      !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
      !!
      !!      Add this trend to the general trend
      !!         (pta) = (pta) + ( difft )
      !!
      !! ** Action : - Update pta arrays with the before iso-level
      !!               biharmonic mixing trend.
      !!----------------------------------------------------------------------
      USE oce     , ONLY:   ztu  => ua       , ztv  => va                           ! (ua,va) used as workspace
      !!
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgu , pgv  ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      !!
      INTEGER  ::  ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::  zbtr, ztra       ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::  zeeu, zeev, zlt
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_ldf_bilap')
      !
      CALL wrk_alloc( jpi, jpj, zeeu, zeev, zlt ) 
      !

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_bilap : iso-level biharmonic operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         !                                               
         DO jk = 1, jpkm1                                        ! Horizontal slab
            !                                             
            !                          !==  Initialization of metric arrays (for z- or s-coordinates)  ==!
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zeeu(ji,jj) = re2u_e1u(ji,jj) * e3u_0(ji,jj,jk) * umask(ji,jj,jk)
                  zeev(ji,jj) = re1v_e2v(ji,jj) * e3v_0(ji,jj,jk) * vmask(ji,jj,jk)
               END DO
            END DO
            !                          !==  Laplacian  ==!
            !
            DO jj = 1, jpjm1                 ! First derivative (gradient)
               DO ji = 1, jpim1   ! vector opt.
                  ztu(ji,jj,jk) = zeeu(ji,jj) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zeev(ji,jj) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
            !
            IF( ln_zps ) THEN                ! set gradient at partial step level (last ocean level)
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( mbku(ji,jj) == jk )  ztu(ji,jj,jk) = zeeu(ji,jj) * pgu(ji,jj,jn)
                     IF( mbkv(ji,jj) == jk )  ztv(ji,jj,jk) = zeev(ji,jj) * pgv(ji,jj,jn)
                  END DO
               END DO
            ENDIF
            ! (ISH)
            IF( ln_zps .AND. ln_isfcav ) THEN ! set gradient at partial step level (first ocean level in a cavity)
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( miku(ji,jj) == MAX(jk,2) )  ztu(ji,jj,jk) = zeeu(ji,jj) * pgui(ji,jj,jn)
                     IF( mikv(ji,jj) == MAX(jk,2) )  ztu(ji,jj,jk) = zeev(ji,jj) * pgvi(ji,jj,jn)
                  END DO
               END DO
            ENDIF
            !
            DO jj = 2, jpjm1                 ! Second derivative (divergence) time the eddy diffusivity coefficient
               DO ji = 2, jpim1   ! vector opt.
                  zbtr = 1.0 / ( e12t(ji,jj) * e3t_0(ji,jj,jk) )
                  zlt(ji,jj) = rldf * ahtt(ji,jj) * zbtr * (   ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                     &                                     + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)   )
               END DO
            END DO
            CALL lbc_lnk( zlt, 'T', 1. )     ! Lateral boundary conditions (unchanged sgn)

            !                          !==  Bilaplacian  ==!
            !
            DO jj = 1, jpjm1                 ! third derivative (gradient)
               DO ji = 1, jpim1   ! vector opt.
                  ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
                  ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
               END DO
            END DO
            DO jj = 2, jpjm1                 ! fourth derivative (divergence) and add to the general tracer trend
               DO ji = 2, jpim1   ! vector opt.
                  ! horizontal diffusive trends
                  zbtr = 1.0 / ( e12t(ji,jj) * e3t_0(ji,jj,jk) )
                  ztra = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
                  ! add it to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
            !                                             
         END DO                                           ! Horizontal slab
         !                                                
         ! "zonal" mean lateral diffusive heat and salt transport
         IF( cdtype == 'TRA' .AND. ln_diaptr ) THEN  
           IF( jn == jp_tem )  htr_ldf(:) = ptr_sj( ztv(:,:,:) )
           IF( jn == jp_sal )  str_ldf(:) = ptr_sj( ztv(:,:,:) )
         ENDIF
         !                                                ! ===========
      END DO                                              ! tracer loop
      !                                                   ! ===========
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_ldf_bilap')
      !
      CALL wrk_dealloc( jpi, jpj, zeeu, zeev, zlt ) 
      !
   END SUBROUTINE tra_ldf_bilap

   !!==============================================================================
END MODULE traldf_bilap
