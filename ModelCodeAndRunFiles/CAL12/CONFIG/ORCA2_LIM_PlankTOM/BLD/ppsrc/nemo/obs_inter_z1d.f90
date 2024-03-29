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

MODULE obs_inter_z1d
   !!======================================================================
   !!                       ***  MODULE obs_inter_z1d  ***
   !! Observation diagnostics: Perform the vertical interpolation
   !!                          from model grid to observation location
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_int_z1d     : Vertical interpolation to the observation point
   !!   obs_int_z1d_spl : Compute the vertical 2nd derivative of the
   !!                     interpolating function for a cubic spline (n1dint=1)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : &  ! Precision variables
      & wp

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_int_z1d,    &  ! Vertical interpolation to the observation pt.
      &   obs_int_z1d_spl    ! Compute the vertical 2nd derivative of the
                             ! interpolating function used with a cubic spline

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_inter_z1d.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obsinter_z1d.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE obs_int_z1d( kpk, kkco, k1dint, kdep, &
      &                    pobsdep, pobsk, pobs2k,  &
      &                    pobs, pdep, pobsmask )
      !!---------------------------------------------------------------------
      !!
      !!                   ***  ROUTINE obs_int_z1d ***
      !!
      !! ** Purpose : Vertical interpolation to the observation point.
      !!  
      !! ** Method  : If k1dint = 0 then use linear interpolation.
      !!              If k1dint = 1 then use cubic spline interpolation.
      !! 
      !! ** Action  :
      !!
      !! References :
      !!
      !! History
      !!      ! 97-11 (A. Weaver, S. Ricci, N. Daget)
      !!      ! 06-03 (G. Smith) Conversion to F90 for use with NEMOVAR
      !!      ! 06-10 (A. Weaver) Cleanup
      !!      ! 07-01 (K. Mogensen) Use profile rather than single level
      !!---------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT(IN) :: kpk        ! Number of vertical levels
      INTEGER, INTENT(IN) :: k1dint     ! 0 = linear; 1 = cubic spline interpolation 
      INTEGER, INTENT(IN) :: kdep       ! Number of levels in profile
      INTEGER, INTENT(IN), DIMENSION(kdep) :: &
         & kkco                 ! Array indicies for interpolation
      REAL(KIND=wp), INTENT(IN), DIMENSION(kdep) :: &
         & pobsdep              ! Depth of the observation
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpk) :: &
         & pobsk,  &            ! Model profile at a given (lon,lat)
         & pobs2k, &            ! 2nd derivative of the interpolating function
         & pdep,   &            ! Model depth array
         & pobsmask             ! Vertical mask
      REAL(KIND=wp), INTENT(OUT), DIMENSION(kdep) :: &
         & pobs                 ! Model equivalent at observation point
  
      !! * Local declarations
      REAL(KIND=wp) :: z1dm       ! Distance above and below obs to model grid points
      REAL(KIND=wp) :: z1dp         
      REAL(KIND=wp) :: zsum       ! Dummy variables for computation
      REAL(KIND=wp) :: zsum2
      INTEGER :: jdep             ! Observation depths loop variable
    
      !------------------------------------------------------------------------
      ! Loop over all observation depths
      !------------------------------------------------------------------------

      DO jdep = 1, kdep

         !---------------------------------------------------------------------
         ! Initialization
         !---------------------------------------------------------------------
         z1dm = ( pdep(kkco(jdep)) - pobsdep(jdep)      )
         z1dp = ( pobsdep(jdep)    - pdep(kkco(jdep)-1) )
         IF ( pobsmask(kkco(jdep)) == 0.0_wp ) z1dp = 0.0_wp

         zsum = z1dm + z1dp
         
         IF ( k1dint == 0 ) THEN

            !-----------------------------------------------------------------
            !  Linear interpolation
            !-----------------------------------------------------------------
            pobs(jdep) = (   z1dm * pobsk(kkco(jdep)-1) &
               &           + z1dp * pobsk(kkco(jdep)  ) ) / zsum

         ELSEIF ( k1dint == 1 ) THEN

            !-----------------------------------------------------------------
            ! Cubic spline interpolation
            !-----------------------------------------------------------------
            zsum2 = zsum * zsum
            pobs(jdep)  = (  z1dm                             * pobsk (kkco(jdep)-1) &
               &           + z1dp                             * pobsk (kkco(jdep)  ) &
               &           + ( z1dm * ( z1dm * z1dm - zsum2 ) * pobs2k(kkco(jdep)-1) &
               &           +   z1dp * ( z1dp * z1dp - zsum2 ) * pobs2k(kkco(jdep)  ) &
               &             ) / 6.0_wp                                              &
               &          ) / zsum

         ENDIF
      END DO

   END SUBROUTINE obs_int_z1d

   SUBROUTINE obs_int_z1d_spl( kpk, pobsk, pobs2k, &
      &                        pdep, pobsmask )
      !!--------------------------------------------------------------------
      !!
      !!                  *** ROUTINE obs_int_z1d_spl ***
      !!
      !! ** Purpose : Compute the local vector of vertical second-derivatives 
      !!              of the interpolating function used with a cubic spline.
      !!  
      !! ** Method  : 
      !! 
      !!    Top and bottom boundary conditions on the 2nd derivative are
      !!    set to zero.
      !!
      !! ** Action  :
      !!
      !! References : 
      !!
      !! History
      !!      ! 01-11 (A. Weaver, S. Ricci, N. Daget)
      !!      ! 06-03 (G. Smith) Conversion to F90 for use with NEMOVAR
      !!      ! 06-10 (A. Weaver) Cleanup
      !!----------------------------------------------------------------------
     
      !! * Arguments
      INTEGER, INTENT(IN) :: kpk               ! Number of vertical levels
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpk) :: &
         & pobsk, &          ! Model profile at a given (lon,lat)
         & pdep,  &          ! Model depth array
         & pobsmask          ! Vertical mask
      REAL(KIND=wp), INTENT(OUT), DIMENSION(kpk) :: &
         & pobs2k            ! 2nd derivative of the interpolating function
  
      !! * Local declarations
      INTEGER :: jk
      REAL(KIND=wp) :: za
      REAL(KIND=wp) :: zb
      REAL(KIND=wp) :: zc
      REAL(KIND=wp) :: zpa
      REAL(KIND=wp) :: zkm
      REAL(KIND=wp) :: zkp
      REAL(KIND=wp) :: zk
      REAL(KIND=wp), DIMENSION(kpk-1) :: &
         & zs, &
         & zp, &
         & zu, &
         & zv

      !-----------------------------------------------------------------------
      ! Matrix initialisation
      !-----------------------------------------------------------------------
      zs(1) =  0.0_wp
      zp(1) =  0.0_wp
      zv(1) = -0.5_wp
      DO jk = 2, kpk-1
         zs(jk) =  ( pdep(jk  ) - pdep(jk-1) ) &
            &    / ( pdep(jk+1) - pdep(jk-1) )
         zp(jk) = zs(jk) * zv(jk-1) + 2.0_wp
         zv(jk) = ( zs(jk) - 1.0_wp ) / zp(jk)
      END DO
 
      !-----------------------------------------------------------------------
      ! Solution of the tridiagonal system
      !-----------------------------------------------------------------------
 
      ! Top boundary condition
      zu(1) = 0.0_wp
 
      DO jk = 2, kpk-1
         za = pdep(jk+1) - pdep(jk-1)
         zb = pdep(jk+1) - pdep(jk  )
         zc = pdep(jk  ) - pdep(jk-1)
 
         zpa = 6.0_wp / ( zp(jk) * za )
         zkm = zpa / zc
         zkp = zpa / zb
         zk  = - ( zkm + zkp )
  
         zu(jk) =  pobsk(jk+1) * zkp  &
            &    + pobsk(jk  ) * zk   &
            &    + pobsk(jk-1) * zkm  &
            &    + zu(jk-1) * ( -zs(jk) / zp(jk) )
      END DO
 
      !-----------------------------------------------------------------------
      ! Second derivative
      !-----------------------------------------------------------------------
      pobs2k(kpk) = 0.0_wp
 
      ! Bottom boundary condition
      DO jk = kpk-1, 1, -1
         pobs2k(jk) = zv(jk) * pobs2k(jk+1) + zu(jk)
         IF ( pobsmask(jk+1) == 0.0_wp ) pobs2k(jk) = 0.0_wp
      END DO
 
  END SUBROUTINE obs_int_z1d_spl



END MODULE obs_inter_z1d

