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

MODULE solmat
   !!======================================================================
   !!                       ***  MODULE  solmat  ***
   !! solver       : construction of the matrix 
   !!======================================================================
   !! History :   1.0  ! 1988-04  (G. Madec)  Original code
   !!                  ! 1993-03  (M. Guyon)  symetrical conditions
   !!                  ! 1993-06  (M. Guyon)  suppress pointers
   !!                  ! 1996-05  (G. Madec)  merge sor and pcg formulations
   !!                  ! 1996-11  (A. Weaver)  correction to preconditioning
   !!   NEMO      1.0  ! 1902-08  (G. Madec)  F90: Free form
   !!              -   ! 1902-11  (C. Talandier, A-M. Treguier) Free surface & Open boundaries
   !!             2.0  ! 2005-09  (R. Benshila)  add sol_exd for extra outer halo
   !!              -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!             3.2  ! 2009-06  (S. Masson)  distributed restart using iom
   !!              -   ! 2009-07  (R. Benshila)  suppression of rigid-lid option
   !!             3.3  ! 2010-09  (D. Storkey) update for BDY module.
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sol_mat : Construction of the matrix of used by the elliptic solvers
   !!   sol_exd :
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE sol_oce         ! ocean solver
   USE phycst          ! physical constants
   USE bdy_oce         ! unstructured open boundary conditions
   USE lbclnk          ! lateral boudary conditions
   USE lib_mpp         ! distributed memory computing
   USE c1d               ! 1D vertical configuration
   USE in_out_manager  ! I/O manager
   USE timing          ! timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sol_mat    ! routine called by inisol.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: solmat.F90 4328 2013-12-06 10:25:13Z davestorkey $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sol_mat( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_mat  ***
      !!
      !! ** Purpose :   Construction of the matrix of used by the elliptic 
      !!              solvers (either sor or pcg methods).
      !!
      !! ** Method  :   The matrix is built for the divergence of the transport 
      !!              system. a diagonal preconditioning matrix is also defined.
      !! 
      !! ** Action  : - gcp    : extra-diagonal elements of the matrix
      !!              - gcdmat : preconditioning matrix (diagonal elements)
      !!              - gcdprc : inverse of the preconditioning matrix
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt
      !!
      INTEGER ::   ji, jj                    ! dummy loop indices
      REAL(wp) ::   zcoefs, zcoefw, zcoefe, zcoefn  ! temporary scalars
      REAL(wp) ::   z2dt, zcoef
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sol_mat')
      !
      
      ! 1. Construction of the matrix
      ! -----------------------------
      zcoef = 0.e0                          ! initialize to zero
      gcp(:,:,1) = 0.e0
      gcp(:,:,2) = 0.e0
      gcp(:,:,3) = 0.e0
      gcp(:,:,4) = 0.e0
      !
      gcdprc(:,:) = 0.e0
      gcdmat(:,:) = 0.e0
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   z2dt = rdt
      ELSE                                        ;   z2dt = 2. * rdt
      ENDIF


      DO jj = 2, jpjm1                      ! matrix of free surface elliptic system
         DO ji = 2, jpim1
            zcoef = z2dt * z2dt * grav * bmask(ji,jj)
            zcoefs = -zcoef * hv(ji  ,jj-1) * e1v(ji  ,jj-1) / e2v(ji  ,jj-1)    ! south coefficient
            zcoefw = -zcoef * hu(ji-1,jj  ) * e2u(ji-1,jj  ) / e1u(ji-1,jj  )    ! west coefficient
            zcoefe = -zcoef * hu(ji  ,jj  ) * e2u(ji  ,jj  ) / e1u(ji  ,jj  )    ! east coefficient
            zcoefn = -zcoef * hv(ji  ,jj  ) * e1v(ji  ,jj  ) / e2v(ji  ,jj  )    ! north coefficient
            gcp(ji,jj,1) = zcoefs
            gcp(ji,jj,2) = zcoefw
            gcp(ji,jj,3) = zcoefe
            gcp(ji,jj,4) = zcoefn
            gcdmat(ji,jj) = e1t(ji,jj) * e2t(ji,jj) * bmask(ji,jj)    &          ! diagonal coefficient
               &          - zcoefs -zcoefw -zcoefe -zcoefn
         END DO
      END DO


      IF( .NOT. Agrif_Root() ) THEN
         !
         IF( nbondi == -1 .OR. nbondi == 2 )   bmask(2     ,:     ) = 0.e0
         IF( nbondi ==  1 .OR. nbondi == 2 )   bmask(nlci-1,:     ) = 0.e0
         IF( nbondj == -1 .OR. nbondj == 2 )   bmask(:     ,2     ) = 0.e0
         IF( nbondj ==  1 .OR. nbondj == 2 )   bmask(:     ,nlcj-1) = 0.e0
         !
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zcoef = z2dt * z2dt * grav * bmask(ji,jj)
               !  south coefficient
               IF( ( nbondj == -1 .OR. nbondj == 2 ) .AND. ( jj == 3 ) ) THEN
                  zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)*(1.-vmask(ji,jj-1,1))
               ELSE
                  zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)
               END IF
               gcp(ji,jj,1) = zcoefs
               ! 
               !  west coefficient
               IF( ( nbondi == -1 .OR. nbondi == 2 ) .AND. ( ji == 3 )  ) THEN
                  zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)*(1.-umask(ji-1,jj,1))
               ELSE
                  zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)
               END IF
               gcp(ji,jj,2) = zcoefw
               !
               !   east coefficient
               IF( ( nbondi == 1 .OR. nbondi == 2 ) .AND. ( ji == nlci-2 ) ) THEN
                  zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)*(1.-umask(ji,jj,1))
               ELSE
                  zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)
               END IF
               gcp(ji,jj,3) = zcoefe
               !
               !   north coefficient
               IF( ( nbondj == 1 .OR. nbondj == 2 ) .AND. ( jj == nlcj-2 ) ) THEN
                  zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)*(1.-vmask(ji,jj,1))
               ELSE
                  zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)
               END IF
               gcp(ji,jj,4) = zcoefn
               !
               ! diagonal coefficient
               gcdmat(ji,jj) = e1t(ji,jj)*e2t(ji,jj)*bmask(ji,jj)   &
                  &            - zcoefs -zcoefw -zcoefe -zcoefn
            END DO
         END DO
         ! 
      ENDIF

      ! 2. Boundary conditions 
      ! ----------------------
      
      ! Cyclic east-west boundary conditions
      !     ji=2 is the column east of ji=jpim1 and reciprocally,
      !     ji=jpim1 is the column west of ji=2
      !     all the coef are already set to zero as bmask is initialized to
      !     zero for ji=1 and ji=jpj in dommsk.
      
      ! Symetrical conditions
      ! free surface: no specific action
      ! bsf system: n-s gradient of bsf = 0 along j=2 (perhaps a bug !!!!!!)
      ! the diagonal coefficient of the southern grid points must be modify to
      ! account for the existence of the south symmetric bassin.
      
      ! North fold boundary condition
      !     all the coef are already set to zero as bmask is initialized to
      !     zero on duplicated lignes and portion of lignes
      
      ! 3. Preconditioned matrix
      ! ------------------------
      
      ! SOR and PCG solvers
      IF( lk_c1d ) CALL lbc_lnk( gcdmat, 'T', 1._wp ) ! 1D case bmask =/0  but gcdmat not define everywhere 
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( bmask(ji,jj) /= 0.e0 )   gcdprc(ji,jj) = 1.e0 / gcdmat(ji,jj)
         END DO
      END DO
         
      gcp(:,:,1) = gcp(:,:,1) * gcdprc(:,:)
      gcp(:,:,2) = gcp(:,:,2) * gcdprc(:,:)
      gcp(:,:,3) = gcp(:,:,3) * gcdprc(:,:)
      gcp(:,:,4) = gcp(:,:,4) * gcdprc(:,:)
      IF( nn_solv == 2 )  gccd(:,:) = rn_sor * gcp(:,:,2)

      IF( nn_solv == 2 .AND. MAX( jpr2di, jpr2dj ) > 0) THEN
         CALL lbc_lnk_e( gcp   (:,:,1), c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary conditions
         CALL lbc_lnk_e( gcp   (:,:,2), c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary conditions
         CALL lbc_lnk_e( gcp   (:,:,3), c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary conditions
         CALL lbc_lnk_e( gcp   (:,:,4), c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary conditions
         CALL lbc_lnk_e( gcdprc(:,:)  , c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary conditions
         CALL lbc_lnk_e( gcdmat(:,:)  , c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary conditions         
         IF( npolj /= 0 ) CALL sol_exd( gcp , c_solver_pt ) ! switch northernelements
      END IF
   
      ! 4. Initialization the arrays used in pcg
      ! ----------------------------------------
      gcb  (:,:) = 0.e0
      gcr  (:,:) = 0.e0
      gcdes(:,:) = 0.e0
      gccd (:,:) = 0.e0
      ! 
      IF( nn_timing == 1 )  CALL timing_stop('sol_mat')
      !
   END SUBROUTINE sol_mat


   SUBROUTINE sol_exd( pt3d, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  routine sol_exd  ***
      !!                  
      !! ** Purpose :   Reorder gcb coefficient on the extra outer  halo 
      !!                at north fold in case of T or F pivot
      !!
      !! ** Method  :   Perform a circular permutation of the coefficients on 
      !!                the total area strictly above the pivot point,
      !!                and on the semi-row of the pivot point   
      !!----------------------------------------------------------------------
      CHARACTER(len=1) , INTENT( in ) ::   cd_type   ! define the nature of pt2d array grid-points
         !                                           !  = T , U , V , F , W 
         !                                           !  = S : T-point, north fold treatment
         !                                           !  = G : F-point, north fold treatment
         !                                           !  = I : sea-ice velocity at F-point with index shift
      REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4), INTENT(inout) ::   pt3d   ! 2D field to be treated
      !!
      INTEGER  ::   ji, jk   ! dummy loop indices
      INTEGER  ::   iloc     ! local integers
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ztab   ! workspace allocated one for all
      !!----------------------------------------------------------------------

      IF( .NOT. ALLOCATED( ztab ) ) THEN
         ALLOCATE( ztab(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4), STAT=iloc )
         IF( lk_mpp    )   CALL mpp_sum ( iloc )
         IF( iloc /= 0 )   CALL ctl_stop('STOP', 'sol_exd: failed to allocate array')
      ENDIF
      
      ztab = pt3d

      SELECT CASE ( npolj )            ! north fold type
      ! 
      CASE ( 3 , 4 )                        !==  T pivot  ==!
         iloc = jpiglo/2 +1 
         !   
         SELECT CASE ( cd_type )
         ! 
         CASE ( 'T' , 'U', 'W' )
            DO jk = 1, 4
               DO ji = 1-jpr2di, nlci+jpr2di
                  pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
               END DO
            END DO
            DO jk =1, 4
               DO ji = nlci+jpr2di, 1-jpr2di,  -1
                  IF( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
                     & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
                     pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
               END DO
            END DO
            !
         CASE ( 'F' , 'I', 'V' )
            DO jk =1, 4
               DO ji = 1-jpr2di, nlci+jpr2di
                  pt3d(ji,nlcj-1:nlcj+jpr2dj,jk) = ztab(ji,nlcj-1:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
               END DO
            END DO
            !
         END SELECT   ! cd_type
          ! 
      CASE ( 5 , 6 )                        !==  F pivot  ==!
         iloc=jpiglo/2
         !
         SELECT CASE (cd_type )
         !
         CASE ( 'T' , 'U', 'W')
            DO jk =1, 4
               DO ji = 1-jpr2di, nlci+jpr2di
                  pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
               END DO
            END DO
            !
         CASE ( 'F' , 'I', 'V' )
            DO jk =1, 4
               DO ji = 1-jpr2di, nlci+jpr2di
                  pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
               END DO
            END DO
            DO jk =1, 4
               DO ji = nlci+jpr2di, 1-jpr2di,  -1
                  IF( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) .OR. ( mi0(iloc) == jpi+1 ) )   EXIT
                    pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
               END DO
            END DO
            !
         END SELECT   ! cd_type
         !
      END SELECT   ! npolj
      !   
   END SUBROUTINE sol_exd

   !!======================================================================
END MODULE solmat
