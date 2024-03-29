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

MODULE dynspg_flt
   !!======================================================================
   !!                   ***  MODULE  dynspg_flt  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
   !! History    OPA  !  1998-05  (G. Roullet)  free surface
   !!                 !  1998-10  (G. Madec, M. Imbard)  release 8.2
   !!   NEMO     O.1  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!            1.0  !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!             -   !  2006-08  (J.Chanut, A.Sellar) Calls to BDY routines. 
   !!            3.2  !  2009-03  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !!            3.7  !  2014-04  (F. Roquet, G. Madec)  add some trends diag
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_dynspg_flt'                              filtered free surface
   !!----------------------------------------------------------------------
   !!   dyn_spg_flt  : update the momentum trend with the surface pressure gradient in the filtered free surface case 
   !!   flt_rst      : read/write the time-splitting restart fields in the ocean restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 
   USE zdf_oce         ! ocean vertical physics
   USE sbc_oce         ! surface boundary condition: ocean
   USE bdy_oce         ! Lateral open boundary condition
   USE sol_oce         ! ocean elliptic solver
   USE phycst          ! physical constants
   USE domvvl          ! variable volume
   USE dynadv          ! advection 
   USE solmat          ! matrix construction for elliptic solvers
   USE solpcg          ! preconditionned conjugate gradient solver
   USE solsor          ! Successive Over-relaxation solver
   USE bdydyn          ! ocean open boundary condition on dynamics
   USE bdyvol          ! ocean open boundary condition (bdy_vol routine)
   USE cla             ! cross land advection
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE wrk_nemo        ! Memory Allocation
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom
   USE lib_fortran
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_spg_flt  ! routine called by step.F90
   PUBLIC   flt_rst      ! routine called by istate.F90

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
   !! $Id: dynspg_flt.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_spg_flt( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_flt  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure 
      !!      gradient in case of filtered free surface formulation  and add
      !!      it to the general trend of momentum equation.
      !!
      !! ** Method  :   Filtered free surface formulation. The surface
      !!      pressure gradient is given by:
      !!         spgu = 1/rau0 d/dx(ps) =  1/e1u di( sshn + btda )
      !!         spgv = 1/rau0 d/dy(ps) =  1/e2v dj( sshn + btda )
      !!      where sshn is the free surface elevation and btda is the after
      !!      time derivative of the free surface elevation
      !!       -1- evaluate the surface presure trend (including the addi-
      !!      tional force) in three steps:
      !!        a- compute the right hand side of the elliptic equation:
      !!            gcb = 1/(e1t e2t) [ di(e2u spgu) + dj(e1v spgv) ]
      !!         where (spgu,spgv) are given by:
      !!            spgu = vertical sum[ e3u (ub+ 2 rdt ua ) ]
      !!                 - grav 2 rdt hu /e1u di[sshn + (emp-rnf)]
      !!            spgv = vertical sum[ e3v (vb+ 2 rdt va) ]
      !!                 - grav 2 rdt hv /e2v dj[sshn + (emp-rnf)]
      !!         and define the first guess from previous computation :
      !!            zbtd = btda
      !!            btda = 2 zbtd - btdb
      !!            btdb = zbtd
      !!        b- compute the relative accuracy to be reached by the
      !!         iterative solver
      !!        c- apply the solver by a call to sol... routine
      !!       -2- compute and add the free surface pressure gradient inclu-
      !!      ding the additional force used to stabilize the equation.
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References : Roullet and Madec, JGR, 2000.
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(  out) ::   kindic   ! solver convergence flag (<0 if not converge)
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   z2dt, z2dtg, zgcb, zbtd, ztdgu, ztdgv   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      REAL(wp), POINTER, DIMENSION(:,:)   ::  zpw
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_flt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_flt : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (free surface constant volume case)'
       
         ! set to zero free surface specific arrays
         spgu(:,:) = 0._wp                     ! surface pressure gradient (i-direction)
         spgv(:,:) = 0._wp                     ! surface pressure gradient (j-direction)

         ! read filtered free surface arrays in restart file
         ! when using agrif, sshn, gcx have to be read in istate
         IF(.NOT. lk_agrif)   CALL flt_rst( nit000, 'READ' )      ! read or initialize the following fields:
         !                                                        ! gcx, gcxb
      ENDIF

      ! Local constant initialization
      z2dt = 2. * rdt                                             ! time step: leap-frog
      IF( neuler == 0 .AND. kt == nit000   )   z2dt = rdt         ! time step: Euler if restart from rest
      IF( neuler == 0 .AND. kt == nit000+1 )   CALL sol_mat( kt )
      z2dtg  = grav * z2dt

      ! Evaluate the masked next velocity (effect of the additional force not included)
      ! ---------------------------------  
      IF( lk_vvl ) THEN          ! variable volume  (surface pressure gradient already included in dyn_hpg)
         !
         IF( ln_dynadv_vec ) THEN      ! vector form : applied on velocity
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ua(ji,jj,jk) = (  ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)  ) * umask(ji,jj,jk)
                     va(ji,jj,jk) = (  vb(ji,jj,jk) + z2dt * va(ji,jj,jk)  ) * vmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         ELSE                          ! flux form : applied on thickness weighted velocity
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ua(ji,jj,jk) = (        ub(ji,jj,jk) * e3u_0(ji,jj,jk)      &
                        &           + z2dt * ua(ji,jj,jk) * e3u_0(ji,jj,jk)  )   &
                        &         / e3u_0(ji,jj,jk) * umask(ji,jj,jk)
                     va(ji,jj,jk) = (        vb(ji,jj,jk) * e3v_0(ji,jj,jk)      &
                        &           + z2dt * va(ji,jj,jk) * e3v_0(ji,jj,jk)  )   &
                        &         / e3v_0(ji,jj,jk) * vmask(ji,jj,jk)
                 END DO
               END DO
            END DO
            !
         ENDIF
         !
      ELSE                       ! fixed volume  (add the surface pressure gradient + unweighted time stepping)
         !
         DO jj = 2, jpjm1              ! Surface pressure gradient (now)
            DO ji = 2, jpim1   ! vector opt.
               spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
               spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
            END DO 
         END DO 
         DO jk = 1, jpkm1              ! unweighted time stepping 
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ua(ji,jj,jk) = (  ub(ji,jj,jk) + z2dt * ( ua(ji,jj,jk) + spgu(ji,jj) )  ) * umask(ji,jj,jk)
                  va(ji,jj,jk) = (  vb(ji,jj,jk) + z2dt * ( va(ji,jj,jk) + spgv(ji,jj) )  ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         IF( l_trddyn )   THEN                      ! temporary save of spg trends
            CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv )
            DO jk = 1, jpkm1              ! unweighted time stepping 
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ztrdu(ji,jj,jk) = spgu(ji,jj) * umask(ji,jj,jk)
                     ztrdv(ji,jj,jk) = spgv(ji,jj) * vmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_spgexp, kt )
         ENDIF
         !
      ENDIF

      IF( nn_cla == 1 .AND. cp_cfg == 'orca' .AND. jp_cfg == 2 )   CALL cla_dynspg( kt )      ! Cross Land Advection (update (ua,va))

      ! compute the next vertically averaged velocity (effect of the additional force not included)
      ! ---------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            spgu(ji,jj) = e3u_0(ji,jj,1) * ua(ji,jj,1)
            spgv(ji,jj) = e3v_0(ji,jj,1) * va(ji,jj,1)
         END DO
      END DO
      DO jk = 2, jpkm1                     ! vertical sum
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               spgu(ji,jj) = spgu(ji,jj) + e3u_0(ji,jj,jk) * ua(ji,jj,jk)
               spgv(ji,jj) = spgv(ji,jj) + e3v_0(ji,jj,jk) * va(ji,jj,jk)
            END DO
         END DO
      END DO

      DO jj = 2, jpjm1                     ! transport: multiplied by the horizontal scale factor
         DO ji = 2, jpim1   ! vector opt.
            spgu(ji,jj) = spgu(ji,jj) * e2u(ji,jj)
            spgv(ji,jj) = spgv(ji,jj) * e1v(ji,jj)
         END DO
      END DO
      CALL lbc_lnk( spgu, 'U', -1. )       ! lateral boundary conditions 
      CALL lbc_lnk( spgv, 'V', -1. )

      IF( lk_vvl ) CALL sol_mat( kt )      ! build the matrix at kt (vvl case only)

      ! Right hand side of the elliptic equation and first guess
      ! --------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            ! Divergence of the after vertically averaged velocity
            zgcb =  spgu(ji,jj) - spgu(ji-1,jj)   &
                  + spgv(ji,jj) - spgv(ji,jj-1)
            gcb(ji,jj) = gcdprc(ji,jj) * zgcb
            ! First guess of the after barotropic transport divergence
            zbtd = gcx(ji,jj)
            gcx (ji,jj) = 2. * zbtd   - gcxb(ji,jj)
            gcxb(ji,jj) =      zbtd
         END DO
      END DO
      ! applied the lateral boundary conditions
      IF( nn_solv == 2 .AND. MAX( jpr2di, jpr2dj ) > 0 )   CALL lbc_lnk_e( gcb, c_solver_pt, 1., jpr2di, jpr2dj )   



      ! Relative precision (computation on one processor)
      ! ------------------
      rnorme =0.e0
      rnorme = GLOB_SUM( gcb(1:jpi,1:jpj) * gcdmat(1:jpi,1:jpj) * gcb(1:jpi,1:jpj) * bmask(:,:) )

      epsr = eps * eps * rnorme
      ncut = 0
      ! if rnorme is 0, the solution is 0, the solver is not called
      IF( rnorme == 0._wp ) THEN
         gcx(:,:) = 0._wp
         res   = 0._wp
         niter = 0
         ncut  = 999
      ENDIF

      ! Evaluate the next transport divergence
      ! --------------------------------------
      !    Iterarive solver for the elliptic equation (except IF sol.=0)
      !    (output in gcx with boundary conditions applied)
      kindic = 0
      IF( ncut == 0 ) THEN
         IF    ( nn_solv == 1 ) THEN   ;   CALL sol_pcg( kindic )      ! diagonal preconditioned conjuguate gradient
         ELSEIF( nn_solv == 2 ) THEN   ;   CALL sol_sor( kindic )      ! successive-over-relaxation
         ENDIF
      ENDIF

      ! Transport divergence gradient multiplied by z2dt
      ! --------------------------------------------====
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            ! trend of Transport divergence gradient
            ztdgu = z2dtg * (gcx(ji+1,jj  ) - gcx(ji,jj) ) / e1u(ji,jj)
            ztdgv = z2dtg * (gcx(ji  ,jj+1) - gcx(ji,jj) ) / e2v(ji,jj)
            ! multiplied by z2dt
            spgu(ji,jj) = z2dt * ztdgu
            spgv(ji,jj) = z2dt * ztdgv
         END DO
      END DO


      IF( l_trddyn )   THEN                     
         ztrdu(:,:,:) = ua(:,:,:)                 ! save the after velocity before the filtered SPG
         ztrdv(:,:,:) = va(:,:,:)
         !
         CALL wrk_alloc( jpi, jpj, zpw )
         !
         zpw(:,:) = - z2dt * gcx(:,:)
         CALL iom_put( "ssh_flt" , zpw )          ! output equivalent ssh modification due to implicit filter
         !
         !                                        ! save surface pressure flux: -pw at z=0
         zpw(:,:) = - rau0 * grav * sshn(:,:) * wn(:,:,1) * tmask(:,:,1)
         CALL iom_put( "pw0_exp" , zpw )
         zpw(:,:) = wn(:,:,1)
         CALL iom_put( "w0" , zpw )
         zpw(:,:) =  rau0 * z2dtg * gcx(:,:) * wn(:,:,1) * tmask(:,:,1)
         CALL iom_put( "pw0_flt" , zpw )
         !
         CALL wrk_dealloc( jpi, jpj, zpw ) 
         !                                   
      ENDIF
      
      ! Add the trends multiplied by z2dt to the after velocity
      ! -------------------------------------------------------
      !     ( c a u t i o n : (ua,va) here are the after velocity not the
      !                       trend, the leap-frog time stepping will not
      !                       be done in dynnxt.F90 routine)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua(ji,jj,jk) = ( ua(ji,jj,jk) + spgu(ji,jj) ) * umask(ji,jj,jk)
               va(ji,jj,jk) = ( va(ji,jj,jk) + spgv(ji,jj) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( l_trddyn )   THEN                      ! save the explicit SPG trends for further diagnostics
         ztrdu(:,:,:) = ( ua(:,:,:) - ztrdu(:,:,:) ) / z2dt
         ztrdv(:,:,:) = ( va(:,:,:) - ztrdv(:,:,:) ) / z2dt
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_spgflt, kt )
         !
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
      ENDIF

      IF( lrst_oce )   CALL flt_rst( kt, 'WRITE' )      ! write filtered free surface arrays in restart file
      !
      IF( nn_timing == 1 )   CALL timing_stop('dyn_spg_flt')
      !
   END SUBROUTINE dyn_spg_flt


   SUBROUTINE flt_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write filtered free surface arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN
         IF( iom_varid( numror, 'gcx', ldstop = .FALSE. ) > 0 ) THEN
! Caution : extra-hallow
! gcx and gcxb are defined as: DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)
            CALL iom_get( numror, jpdom_autoglo, 'gcx' , gcx (1:jpi,1:jpj) )
            CALL iom_get( numror, jpdom_autoglo, 'gcxb', gcxb(1:jpi,1:jpj) )
            IF( neuler == 0 )   gcxb(:,:) = gcx (:,:)
         ELSE
            gcx (:,:) = 0.e0
            gcxb(:,:) = 0.e0
         ENDIF
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
! Caution : extra-hallow
! gcx and gcxb are defined as: DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)
         CALL iom_rstput( kt, nitrst, numrow, 'gcx' , gcx (1:jpi,1:jpj) )
         CALL iom_rstput( kt, nitrst, numrow, 'gcxb', gcxb(1:jpi,1:jpj) )
      ENDIF
      !
   END SUBROUTINE flt_rst

   
   !!======================================================================
END MODULE dynspg_flt
