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

MODULE obs_conv
   !!=====================================================================
   !!                       ***  MODULE  obs_conv  ***
   !! Observation diagnostics: Various conversion functions
   !!=====================================================================
   !!
   !!   potemp   : Compute potential temperature from insitu temperature,
   !!              salinity and pressure
   !!   fspott   : Compute potential temperature from insitu temperature,
   !!              salinity and pressure
   !!   atg      : Compute adiabatic temperature gradient deg c per decibar
   !!   theta    : Compute potential temperature from insitu temperature,
   !!              salinity and pressure
   !!   depth    : Compute depth from pressure and latitude.
   !!   p_to_dep : Compute depth from pressure and latitude 
   !!              (approximate version)
   !!   dep_to_p : Compute pressure from depth and latitude 
   !!              (approximate version)
   !!---------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & wp   
   IMPLICIT NONE
 
   !! * Function accessibility
   PRIVATE
   PUBLIC &
      & potemp,   &
      & fspott,   &
      & atg,      &
      & theta,    &
      & depth,    &
      & p_to_dep, &
      & dep_to_p
   
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_conv.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

    !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_conv_functions.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

  REAL(KIND=wp) FUNCTION potemp( ps, pt, pp, ppr )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION potemp  ***
      !!          
      !! ** Purpose : Compute potential temperature
      !!
      !! ** Method  : A regression formula is used. 
      !!
      !! ** Action  : The code is kept as close to the F77 code as possible
      !!              Check value: potemp(35,20,2000,0) = 19.621967
      !!
      !! References : T. J. Mcdougall, D. R. Jackett, D. G. Wright 
      !!              and R. Feistel
      !!              Accurate and computationally efficient algoritms for
      !!              potential temperatures and density of seawater
      !!              Journal of atmospheric and oceanic technology
      !!              Vol 20, 2003, pp 730-741
      !!              
      !!
      !! History :
      !!        !  07-05 (K. Mogensen) Original code
      !!----------------------------------------------------------------------

      !! * Arguments

      REAL(KIND=wp), INTENT(IN) :: ps
      REAL(KIND=wp), INTENT(IN) :: pt
      REAL(KIND=wp), INTENT(IN) :: pp
      REAL(KIND=wp), INTENT(IN) :: ppr

      !! * Local declarations
      REAL(KIND=wp) :: zpol
      REAL(KIND=wp), PARAMETER :: a1 =  1.067610e-05
      REAL(KIND=wp), PARAMETER :: a2 = -1.434297e-06
      REAL(KIND=wp), PARAMETER :: a3 = -7.566349e-09
      REAL(KIND=wp), PARAMETER :: a4 = -8.535585e-06
      REAL(KIND=wp), PARAMETER :: a5 =  3.074672e-08
      REAL(KIND=wp), PARAMETER :: a6 =  1.918639e-08
      REAL(KIND=wp), PARAMETER :: a7 =  1.788718e-10

      zpol = a1 + a2 * ps + a3 * ( pp + ppr ) + a4 * pt &
         & + a5 * ps * pt + a6 * pt * pt + a7 * pt * ( pp + ppr )
      
      potemp = pt + ( pp - ppr ) * zpol
      
   END FUNCTION potemp

   REAL(KIND=wp) FUNCTION fspott( pft, pfs, pfp )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION fspott  ***
      !!          
      !! ** Purpose : Compute potential temperature
      !!
      !! ** Method  : A regression formula is used. 
      !!
      !! ** Action  : Check value: fspott(10,25,1000) = 8.4678516
      !!
      !! References : A. E. Gill
      !!              Atmosphere-Ocean Dynamics
      !!              Volume 30 (International Geophysics)
      !!
      !! History :
      !!        !  07-05 (K. Mogensen) NEMO adopting of OPAVAR code.
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=wp) :: pft   ! in situ temperature in degrees celcius
      REAL(KIND=wp) :: pfs   ! salinity in psu
      REAL(KIND=wp) :: pfp   ! pressure in bars
      
      fspott = &
         &  pft - pfp * (   (            3.6504e-4                     &
         &                    + pft * (  8.3198e-5                     &
         &                    + pft * ( -5.4065e-7                     &
         &                    + pft *    4.0274e-9  ) ) )              &
         &                + ( pfs - 35.0 ) * (         1.7439e-5       &
         &                                     - pft * 2.9778e-7 )     &
         &                + pfp * (            8.9309e-7               &
         &                          + pft * ( -3.1628e-8               &
         &                          + pft *    2.1987e-10 )            &
         &                          - ( pfs - 35.0 ) *  4.1057e-9      &
         &                          + pfp * (          -1.6056e-10     &
         &                                    + pft * 5.0484e-12 ) ) )

   END FUNCTION fspott

   REAL(KIND=wp) FUNCTION atg( p_s, p_t, p_p )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION atg  ***
      !!          
      !! ** Purpose : Compute adiabatic temperature gradient deg c per decibar
      !!
      !! ** Method  : A regression formula is used
      !!
      !! ** Action  : The code is kept as close to the F77 code as possible
      !!              Check value: atg(40,40,10000) = 3.255974e-4
      !!
      !! References : N. P. Fotonoff and R.C. Millard jr., 
      !!              Algoritms for computation of fundamental
      !!              properties of seawater
      !!              Unesco technical papers in marine science 44
      !!              Unesco 1983
      !!
      !! History :
      !!        !  07-05 (K. Mogensen) Original code based on the F77 code.
      !!----------------------------------------------------------------------

      !! * Arguments

      REAL(KIND=wp), INTENT(IN) :: p_s    ! Salinity in PSU
      REAL(KIND=wp), INTENT(IN) :: p_t    ! Temperature in centigrades
      REAL(KIND=wp), INTENT(IN) :: p_p    ! Pressure in decibars.

      !! * Local declarations
      
      REAL(KIND=wp) :: z_ds

      z_ds = p_s - 35.0
      atg = ((( -2.1687e-16 * p_t + 1.8676e-14 ) * p_t - 4.6206e-13 ) * p_p      &
         &  + (( 2.7759e-12 * p_t - 1.1351e-10 ) * z_ds + (( - 5.4481e-14 * p_t  &  
         &  + 8.733e-12 ) * p_t - 6.7795e-10 ) * p_t  + 1.8741e-8)) * p_p        &
         &  + ( -4.2393e-8 * p_t + 1.8932e-6 ) * z_ds                          &
         &  + (( 6.6228e-10 * p_t - 6.836e-8 ) * p_t + 8.5258e-6 ) * p_t + 3.5803e-5

   END FUNCTION atg
      
   REAL(KIND=wp) FUNCTION theta( p_s, p_t0, p_p0, p_pr )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION theta  ***
      !!          
      !! ** Purpose : Compute potential temperature
      !!
      !! ** Method  : A regression formula is used. 
      !!
      !! ** Action  : The code is kept as close to the F77 code as possible
      !!              Check value: theta(40,40,10000,0) = 36.89073
      !!
      !! References : N. P. Fotonoff and R.C. Millard jr., 
      !!              Algoritms for computation of fundamental
      !!              properties of seawater
      !!              Unesco technical papers in marine science 44
      !!              Unesco 1983
      !!
      !! History :
      !!        !  07-05 (K. Mogensen) Original code based on the F77 code.
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=wp), INTENT(IN) :: p_s
      REAL(KIND=wp), INTENT(IN) :: p_t0
      REAL(KIND=wp), INTENT(IN) :: p_p0
      REAL(KIND=wp), INTENT(IN) :: p_pr

      !! * Local declarations
      REAL(KIND=wp) :: z_p
      REAL(KIND=wp) :: z_t
      REAL(KIND=wp) :: z_h
      REAL(KIND=wp) :: z_xk
      REAL(KIND=wp) :: z_q

      z_p = p_p0
      z_t = p_t0
      z_h = p_pr - z_p
      z_xk = z_h * atg( p_s, z_t, z_p ) 
      Z_t = z_t + 0.5 * z_xk
      z_q = z_xk
      z_p = z_p + 0.5 * z_h
      z_xk = z_h * atg( p_s, z_t, z_p )
      z_t = z_t + 0.29289322 * ( z_xk - z_q )
      z_q = 0.58578644 * z_xk + 0.121320344 * z_q
      z_xk = z_h * atg( p_s, z_t, z_p )
      z_t = z_t + 1.707106781 * ( z_xk - z_q )
      z_q = 3.414213562 * z_xk - 4.121320244 * z_q
      z_p = z_p + 0.5 * z_h
      z_xk = z_h * atg( p_s, z_t, z_p )
      theta = z_t + ( z_xk - 2.0 * z_q ) / 6.0

   END FUNCTION theta

   REAL(KIND=wp) FUNCTION depth( p_p, p_lat )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION depth  ***
      !!          
      !! ** Purpose : Compute depth from pressure and latitudes
      !!
      !! ** Method  : A regression formula is used. 
      !!
      !! ** Action  : The code is kept as close to the F77 code as possible
      !!              Check value: depth(10000,30) = 9712.653
      !!
      !! References : N. P. Fotonoff and R.C. Millard jr., 
      !!              Algoritms for computation of fundamental
      !!              properties of seawater
      !!              Unesco technical papers in marine science 44
      !!              Unesco 1983
      !!
      !! History :
      !!        !  07-05 (K. Mogensen) Original code based on the F77 code.
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=wp), INTENT(IN) :: p_p     ! Pressure in decibars
      REAL(KIND=wp), INTENT(IN) :: p_lat   ! Latitude in degrees

      !! * Local declarations
      REAL(KIND=wp) :: z_x
      REAL(KIND=wp) :: z_gr
      
      z_x = SIN( p_lat / 57.29578 )
      z_x = z_x * z_x
      z_gr = 9.780318 * ( 1.0 + ( 5.2788e-3 + 2.36e-5 * z_x ) * z_x ) + 1.092e-6 * p_p
      depth = ((( -1.82e-15 * p_p + 2.279e-10 ) * p_p - 2.2512e-5 ) * p_p + 9.72659 ) * p_p
      depth = depth / z_gr

   END FUNCTION depth

   REAL(KIND=wp) FUNCTION p_to_dep( p_p, p_lat )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION p_to_dep  ***
      !!          
      !! ** Purpose : Compute depth from pressure and latitudes
      !!
      !! ** Method  : A regression formula is used. This version is less
      !!              accurate the "depth" but invertible.
      !!
      !! ** Action  : 
      !!
      !! References : P.M Saunders
      !!              Pratical conversion of pressure to depth
      !!              Journal of physical oceanography Vol 11, 1981, pp 573-574
      !!
      !! History :
      !!        !  07-05  (K. Mogensen) Original code
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=wp), INTENT(IN) :: p_p    ! Pressure in decibars
      REAL(KIND=wp), INTENT(IN) :: p_lat  ! Latitude in degrees

      !! * Local declarations
      REAL(KIND=wp) :: z_x
      REAL(KIND=wp) :: z_c1
      REAL(KIND=wp) :: z_c2

      z_x = SIN( p_lat / 57.29578 )
      z_x = z_x * z_x
      z_c1 = ( 5.92  + 5.25 * z_x ) * 1e-3 
      z_c2 = 2.21e-6 
      p_to_dep = (1 - z_c1)  * p_p - z_c2 * p_p * p_p

   END FUNCTION p_to_dep

   REAL(KIND=wp) FUNCTION dep_to_p( p_dep, p_lat )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION dep_to_p  ***
      !!          
      !! ** Purpose : Compute depth from pressure and latitudes
      !!
      !! ** Method  : The expression used in p_to_dep is inverted.
      !!
      !! ** Action  : 
      !!
      !! References : P.M Saunders
      !!              Pratical conversion of pressure to depth
      !!              Journal of physical oceanography Vol 11, 1981, pp 573-574
      !!
      !! History :
      !!        !  07-05  (K. Mogensen) Original code
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=wp), INTENT(IN) :: p_dep    ! Depth in meters
      REAL(KIND=wp), INTENT(IN) :: p_lat    ! Latitude in degrees

      !! * Local declarations
      REAL(KIND=wp) :: z_x
      REAL(KIND=wp) :: z_c1
      REAL(KIND=wp) :: z_c2
      REAL(KIND=wp) :: z_d

      z_x = SIN( p_lat / 57.29578 )
      z_x = z_x * z_x
      z_c1 = ( 5.92  + 5.25 * z_x ) * 1e-3 
      z_c2 = 2.21e-6 
      z_d = ( z_c1 - 1 ) * ( z_c1 - 1  ) - 4 * z_c2 * p_dep
      dep_to_p = (( 1 - z_c1 ) - SQRT( z_d )) / ( 2 * z_c2 )

   END FUNCTION dep_to_p

END MODULE obs_conv
