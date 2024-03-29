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

MODULE stopts
   !!==============================================================================
   !!                       ***  MODULE  stopts  ***
   !! Stochastic parameterization: compute stochastic tracer fluctuations
   !!==============================================================================
   !! History :  3.3  ! 2011-12 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sto_pts        : compute current stochastic tracer fluctuations
   !!   sto_pts_init   : initialisation for stochastic tracer fluctuations
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE lbclnk          ! lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE stopar          ! stochastic parameterization

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sto_pts         ! called by step.F90
   PUBLIC   sto_pts_init    ! called by nemogcm.F90

   ! Public array with random tracer fluctuations
   REAL(wp), PUBLIC, DIMENSION(:,:,:,:,:), ALLOCATABLE :: pts_ran

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: eosbn2.F90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sto_pts( pts )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_pts  ***
      !!
      !! ** Purpose :   Compute current stochastic tracer fluctuations
      !!
      !! ** Method  :   Compute tracer differences from a random walk
      !!                around every model grid point
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pts   ! 1 : potential temperature  [Celsius]
      !                                                               ! 2 : salinity               [psu]
      INTEGER  ::   ji, jj, jk, jts, jdof ! dummy loop indices
      INTEGER  ::   jim1, jjm1, jkm1  ! incremented indices
      INTEGER  ::   jip1, jjp1, jkp1  !     -          -
      REAL(wp) ::   zdtsim, zdtsjm, zdtskm         ! temporary scalars
      REAL(wp) ::   zdtsip, zdtsjp, zdtskp, zdts   !     -        -
      !!----------------------------------------------------------------------

      DO jts = 1, jpts
        CALL lbc_lnk( pts(:,:,:,jts), 'T' , 1._wp )
      ENDDO

      DO jdof = 1, nn_sto_eos
        DO jts = 1, jpts
           DO jk = 1, jpkm1
              jkm1 = MAX(jk-1,1) ; jkp1 = MIN(jk+1,jpkm1)
              DO jj = 1, jpj
                 jjm1 = MAX(jj-1,1) ; jjp1 = MIN(jj+1,jpj)
                 DO ji = 1, jpi
                    jim1 = MAX(ji-1,1) ; jip1 = MIN(ji+1,jpi)
                    !
                    ! compute tracer gradient
                    zdtsip = ( pts(jip1,jj,jk,jts) - pts(ji,jj,jk,jts) ) * tmask(jip1,jj,jk)
                    zdtsim = ( pts(ji,jj,jk,jts) - pts(jim1,jj,jk,jts) ) * tmask(jim1,jj,jk)
                    zdtsjp = ( pts(ji,jjp1,jk,jts) - pts(ji,jj,jk,jts) ) * tmask(ji,jjp1,jk)
                    zdtsjm = ( pts(ji,jj,jk,jts) - pts(ji,jjm1,jk,jts) ) * tmask(ji,jjm1,jk)
                    zdtskp = ( pts(ji,jj,jkp1,jts) - pts(ji,jj,jk,jts) ) * tmask(ji,jj,jkp1)
                    zdtskm = ( pts(ji,jj,jk,jts) - pts(ji,jj,jkm1,jts) ) * tmask(ji,jj,jkm1)
                    !
                    ! compute random tracer fluctuation (zdts)
                    zdts   = ( zdtsip + zdtsim ) * sto2d(ji,jj,jsto_eosi(jdof)) + &
                           & ( zdtsjp + zdtsjm ) * sto2d(ji,jj,jsto_eosj(jdof)) + &
                           & ( zdtskp + zdtskm ) * sto2d(ji,jj,jsto_eosk(jdof))
!                   zdts   = zdtsip * MAX(sto2d(ji,jj,jsto_eosi),0._wp) + &
!                          & zdtsim * MIN(sto2d(ji,jj,jsto_eosi),0._wp) + &
!                          & zdtsjp * MAX(sto2d(ji,jj,jsto_eosj),0._wp) + &
!                          & zdtsjm * MIN(sto2d(ji,jj,jsto_eosj),0._wp) + &
!                          & zdtskp * MAX(sto2d(ji,jj,jsto_eosk),0._wp) + &
!                          & zdtskm * MIN(sto2d(ji,jj,jsto_eosk),0._wp)
                    zdts   = zdts  * tmask(ji,jj,jk) *SIN( gphit(ji,jj) * rad )
                    pts_ran(ji,jj,jk,jts,jdof) = zdts * 0.5_wp
                    !
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! Eliminate any possible negative salinity
      DO jdof = 1, nn_sto_eos
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  pts_ran(ji,jj,jk,jp_sal,jdof) = MIN( ABS(pts_ran(ji,jj,jk,jp_sal,jdof)) ,  &
                                                &      MAX(pts(ji,jj,jk,jp_sal),0._wp) )     &
                                                &  * SIGN(1._wp,pts_ran(ji,jj,jk,jp_sal,jdof))
               END DO
            END DO
         END DO
      END DO

      ! Eliminate any temperature lower than -2 degC
!     DO jdof = 1, nn_sto_eos
!        DO jk = 1, jpkm1
!           DO jj = 1, jpj
!              DO ji = 1, jpi
!                 pts_ran(ji,jj,jk,jp_tem,jdof) = MIN( ABS(pts_ran(ji,jj,jk,jp_tem,jdof)) ,    &
!                                               &      MAX(pts(ji,jj,jk,jp_tem)+2._wp,0._wp) ) &
!                                               &  * SIGN(1._wp,pts_ran(ji,jj,jk,jp_tem,jdof))
!              END DO
!           END DO
!        END DO
!     END DO


      ! Lateral boundary conditions on pts_ran
      DO jdof = 1, nn_sto_eos
         DO jts = 1, jpts
            CALL lbc_lnk( pts_ran(:,:,:,jts,jdof), 'T' , 1._wp )
         END DO
      END DO

   END SUBROUTINE sto_pts


   SUBROUTINE sto_pts_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sto_pts_init  ***
      !!
      !! ** Purpose :   Initialisation for stochastic tracer fluctuations
      !!
      !! ** Method  :   Allocate required array
      !!
      !!----------------------------------------------------------------------

      ALLOCATE(pts_ran(jpi,jpj,jpk,jpts,nn_sto_eos))

   END SUBROUTINE sto_pts_init

END MODULE stopts
