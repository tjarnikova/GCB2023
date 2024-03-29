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

MODULE limhdf_2
   !!======================================================================
   !!                    ***  MODULE limhdf_2   ***
   !! LIM 2.0 ice model : horizontal diffusion of sea-ice quantities
   !!======================================================================
   !! History :  LIM  !  2000-01 (LIM) Original code
   !!             -   !  2001-05 (G. Madec, R. Hordoir) opa norm
   !!            1.0  !  2002-08 (C. Ethe)  F90, free form
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2'                                    LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_hdf_2  : diffusion trend on sea-ice variable
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean domain
   USE ice_2            ! LIM-2: ice variables
   USE lbclnk           ! lateral boundary condition - MPP exchanges
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE prtctl           ! Print control
   USE in_out_manager   ! I/O manager
   USE lib_fortran      ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_hdf_2         ! called by limtrp_2.F90

   LOGICAL  ::   linit = .TRUE.   ! ! initialization flag (set to flase after the 1st call)
   REAL(wp) ::   epsi04 = 1e-04   ! constant
   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   efact   ! metric coefficient

   !! * Substitution 
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
   !! NEMO/LIM2 4.0 , UCL - NEMO Consortium (2010)
   !! $Id: limhdf_2.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_hdf_2( ptab )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_hdf_2  ***
      !!
      !! ** purpose :   Compute and add the diffusive trend on sea-ice variables
      !!
      !! ** method  :   Second order diffusive operator evaluated using a
      !!              Cranck-Nicholson time Scheme.
      !!
      !! ** Action  :    update ptab with the diffusive contribution
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   ptab   ! Field on which the diffusion is applied  
      !
      INTEGER  ::   ji, jj            ! dummy loop indices
      INTEGER  ::   its, iter, ierr   ! local integers
      REAL(wp) ::   zalfa, zrlxint, zconv, zeps   ! local scalars
      REAL(wp), DIMENSION(:,:), POINTER :: zrlx, zflu, zflv, zdiv0, zdiv, ztab0 
      CHARACTER (len=55) :: charout
      !!-------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, zrlx, zflu, zflv, zdiv0, zdiv, ztab0 )

      !                       !==  Initialisation  ==!
      !
      IF( linit ) THEN              ! Metric coefficient (compute at the first call and saved in efact)
         ALLOCATE( efact(jpi,jpj) , STAT=ierr )
         IF( lk_mpp    )   CALL mpp_sum( ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'lim_hdf_2 : unable to allocate standard arrays' )
         DO jj = 2, jpjm1  
            DO ji = 2 , jpim1   ! vector opt.
               efact(ji,jj) = ( e2u(ji,jj) + e2u(ji-1,jj) + e1v(ji,jj) + e1v(ji,jj-1) ) / ( e1t(ji,jj) * e2t(ji,jj) )
            END DO
         END DO
         linit = .FALSE.
      ENDIF
      !
      !                             ! Time integration parameters
      zalfa = 0.5_wp                      ! =1.0/0.5/0.0 = implicit/Cranck-Nicholson/explicit
      its   = 100                         ! Maximum number of iteration
      zeps  =  2._wp * epsi04
      !
      ztab0(:, : ) = ptab(:,:)      ! Arrays initialization
      zdiv0(:, 1 ) = 0._wp
      zdiv0(:,jpj) = 0._wp
      zflu (jpi,:) = 0._wp   
      zflv (jpi,:) = 0._wp
      zdiv0(1,  :) = 0._wp
      zdiv0(jpi,:) = 0._wp

      zconv = 1._wp           !==  horizontal diffusion using a Crant-Nicholson scheme  ==!
      iter  = 0
      !
      DO WHILE (  zconv > zeps  .AND.  iter <= its  )    ! Sub-time step loop
         !
         iter = iter + 1                                       ! incrementation of the sub-time step number
         !
         DO jj = 1, jpjm1                                      ! diffusive fluxes in U- and V- direction
            DO ji = 1 , jpim1   ! vector opt.
               zflu(ji,jj) = pahu(ji,jj) * e2u(ji,jj) / e1u(ji,jj) * ( ptab(ji+1,jj) - ptab(ji,jj) )
               zflv(ji,jj) = pahv(ji,jj) * e1v(ji,jj) / e2v(ji,jj) * ( ptab(ji,jj+1) - ptab(ji,jj) )
            END DO
         END DO
         !
         DO jj= 2, jpjm1                                       ! diffusive trend : divergence of the fluxes
            DO ji = 2 , jpim1   ! vector opt. 
               zdiv (ji,jj) = (  zflu(ji,jj) - zflu(ji-1,jj  )   &
                  &            + zflv(ji,jj) - zflv(ji  ,jj-1)  ) / ( e1t (ji,jj) * e2t (ji,jj) )
            END DO
         END DO
         !
         IF( iter == 1 )   zdiv0(:,:) = zdiv(:,:)              ! save the 1st evaluation of the diffusive trend in zdiv0
         !
         DO jj = 2, jpjm1                                      ! iterative evaluation
            DO ji = 2 , jpim1   ! vector opt.
               zrlxint = (   ztab0(ji,jj)    &
                  &       +  rdt_ice * (           zalfa   * ( zdiv(ji,jj) + efact(ji,jj) * ptab(ji,jj) )   &
                  &                      + ( 1.0 - zalfa ) *   zdiv0(ji,jj) )  )                             & 
                  &    / ( 1.0 + zalfa * rdt_ice * efact(ji,jj) )
               zrlx(ji,jj) = ptab(ji,jj) + om * ( zrlxint - ptab(ji,jj) )
            END DO
         END DO
         CALL lbc_lnk( zrlx, 'T', 1. )                         ! lateral boundary condition

         zconv = 0._wp                                         ! convergence test

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zconv = MAX( zconv, ABS( zrlx(ji,jj) - ptab(ji,jj) )  )
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_max( zconv )                  ! max over the global domain

         ptab(:,:) = zrlx(:,:)
         !
      END DO                                             ! end of sub-time step loop

      IF(ln_ctl)   THEN
         zrlx(:,:) = ptab(:,:) - ztab0(:,:)
         WRITE(charout,FMT="(' lim_hdf  : zconv =',D23.16, ' iter =',I4,2X)") zconv, iter
         CALL prt_ctl( tab2d_1=zrlx, clinfo1=charout )
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zrlx, zflu, zflv, zdiv0, zdiv, ztab0 )
      !
   END SUBROUTINE lim_hdf_2


   !!======================================================================
END MODULE limhdf_2
