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

MODULE tide_mod
   !!======================================================================
   !!                       ***  MODULE  tide_mod  ***
   !! Compute nodal modulations corrections and pulsations
   !!======================================================================
   !! History :  1.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE daymod         ! calendar

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tide_harmo       ! called by tideini and diaharm modules
   PUBLIC   tide_init_Wave   ! called by tideini and diaharm modules

   INTEGER, PUBLIC, PARAMETER ::   jpmax_harmo = 19   !: maximum number of harmonic

   TYPE, PUBLIC ::    tide
      CHARACTER(LEN=4) ::   cname_tide
      REAL(wp)         ::   equitide
      INTEGER          ::   nutide
      INTEGER          ::   nt, ns, nh, np, np1, shift
      INTEGER          ::   nksi, nnu0, nnu1, nnu2, R
      INTEGER          ::   nformula
   END TYPE tide

   TYPE(tide), PUBLIC, DIMENSION(jpmax_harmo) ::   Wave   !:

   REAL(wp) ::   sh_T, sh_s, sh_h, sh_p, sh_p1             ! astronomic angles
   REAL(wp) ::   sh_xi, sh_nu, sh_nuprim, sh_nusec, sh_R   !
   REAL(wp) ::   sh_I, sh_x1ra, sh_N                       !

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , LOCEAN-IPSL (2010) 
   !! $Id: tide_mod.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tide_init_Wave
   !!----------------------------------------------------------------------
   !! History :  3.2  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------

   !             !! name_tide , equitide , nutide , nt , ns , nh , np , np1 , shift , nksi , nnu0 , nnu1 , nnu2 , R , formula !!
   !             !!           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !!
   Wave( 1) = tide(  'M2'     , 0.242297 ,    2   ,  2 , -2 ,  2 ,  0 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,    78   )
   Wave( 2) = tide(  'N2'     , 0.046313 ,    2   ,  2 , -3 ,  2 ,  1 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,    78   )
   Wave( 3) = tide( '2N2'     , 0.006184 ,    2   ,  2 , -4 ,  2 ,  2 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,    78   )
   Wave( 4) = tide(  'S2'     , 0.113572 ,    2   ,  2 ,  0 ,  0 ,  0 ,  0  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,     0   )
   Wave( 5) = tide(  'K2'     , 0.030875 ,    2   ,  2 ,  0 ,  2 ,  0 ,  0  ,    0  ,  0   ,  0   ,  0   , -2   , 0 ,   235   )
   !              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
   Wave( 6) = tide(  'K1'     , 0.142408 ,    1   ,  1 ,  0 ,  1 ,  0 ,  0  ,  -90  ,  0   ,  0   , -1   ,  0   , 0 ,   227   )
   Wave( 7) = tide(  'O1'     , 0.101266 ,    1   ,  1 , -2 ,  1 ,  0 ,  0  ,  +90  ,  2   , -1   ,  0   ,  0   , 0 ,    75   )
   Wave( 8) = tide(  'Q1'     , 0.019387 ,    1   ,  1 , -3 ,  1 ,  1 ,  0  ,  +90  ,  2   , -1   ,  0   ,  0   , 0 ,    75   )
   Wave( 9) = tide(  'P1'     , 0.047129 ,    1   ,  1 ,  0 , -1 ,  0 ,  0  ,  +90  ,  0   ,  0   ,  0   ,  0   , 0 ,    0    )
   !              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
   Wave(10) = tide(  'M4'     , 0.000000 ,    4   ,  4 , -4 ,  4 ,  0 ,  0  ,    0  ,  4   , -4   ,  0   ,  0   , 0 ,    1    )
   !              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
   Wave(11) = tide(  'Mf'     , 0.042017 ,    0   ,  0 ,  2 ,  0 ,  0 ,  0  ,    0  , -2   ,  0   ,  0   ,  0   , 0 ,   74    )
   Wave(12) = tide(  'Mm'     , 0.022191 ,    0   ,  0 ,  1 ,  0 , -1 ,  0  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,   73    )
   Wave(13) = tide(  'Msqm'   , 0.000667 ,    0   ,  0 ,  4 , -2 ,  0 ,  0  ,    0  , -2   ,  0   ,  0   ,  0   , 0 ,   74    )
   Wave(14) = tide(  'Mtm'    , 0.008049 ,    0   ,  0 ,  3 ,  0 , -1 ,  0  ,    0  , -2   ,  0   ,  0   ,  0   , 0 ,   74    )
   !              !           !          !        !    !    !    !    !     !       !      !      !      !      !   !         !
   Wave(15) = tide(  'S1'     , 0.000000 ,    1   ,  1 ,  0 ,  0 ,  0 ,  0  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,    0    )   
   Wave(16) = tide(  'MU2'    , 0.005841 ,    2   ,  2 , -4 ,  4 ,  0 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,   78    )
   Wave(17) = tide(  'NU2'    , 0.009094 ,    2   ,  2 , -3 ,  4 , -1 ,  0  ,    0  ,  2   , -2   ,  0   ,  0   , 0 ,   78    ) 
   Wave(18) = tide(  'L2'     , 0.006694 ,    2   ,  2 , -1 ,  2 , -1 ,  0  , +180  ,  2   , -2   ,  0   ,  0   , 0 ,  215    )
   Wave(19) = tide(  'T2'     , 0.006614 ,    2   ,  2 ,  0 , -1 ,  0 ,  1  ,    0  ,  0   ,  0   ,  0   ,  0   , 0 ,    0    )
   END SUBROUTINE tide_init_Wave


   SUBROUTINE tide_harmo( pomega, pvt, put , pcor, ktide ,kc)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(kc), INTENT(in ) ::   ktide            ! Indice of tidal constituents
      INTEGER                , INTENT(in ) ::   kc               ! Total number of tidal constituents
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pomega           ! pulsation in radians/s
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pvt, put, pcor   !
      !!----------------------------------------------------------------------
      !
      CALL astronomic_angle
      CALL tide_pulse( pomega, ktide ,kc )
      CALL tide_vuf  ( pvt, put, pcor, ktide ,kc )
      !
   END SUBROUTINE tide_harmo


   SUBROUTINE astronomic_angle
      !!----------------------------------------------------------------------
      !!  tj is time elapsed since 1st January 1900, 0 hour, counted in julian
      !!  century (e.g. time in days divide by 36525)
      !!----------------------------------------------------------------------
      REAL(wp) ::   cosI, p, q, t2, t4, sin2I, s2, tgI2, P1, sh_tgn2, at1, at2
      REAL(wp) ::   zqy , zsy, zday, zdj, zhfrac
      !!----------------------------------------------------------------------
      !
      zqy = AINT( (nyear-1901.)/4. )
      zsy = nyear - 1900.
      !
      zdj  = dayjul( nyear, nmonth, nday )
      zday = zdj + zqy - 1.
      !
      zhfrac = nsec_day / 3600.
      !
      !----------------------------------------------------------------------
      !  Sh_n Longitude of ascending lunar node
      !----------------------------------------------------------------------
      sh_N=(259.1560564-19.328185764*zsy-.0529539336*zday-.0022064139*zhfrac)*rad
      !----------------------------------------------------------------------
      ! T mean solar angle (Greenwhich time)
      !----------------------------------------------------------------------
      sh_T=(180.+zhfrac*(360./24.))*rad
      !----------------------------------------------------------------------
      ! h mean solar Longitude
      !----------------------------------------------------------------------
      sh_h=(280.1895014-.238724988*zsy+.9856473288*zday+.0410686387*zhfrac)*rad
      !----------------------------------------------------------------------
      ! s mean lunar Longitude
      !----------------------------------------------------------------------
      sh_s=(277.0256206+129.38482032*zsy+13.176396768*zday+.549016532*zhfrac)*rad
      !----------------------------------------------------------------------
      ! p1 Longitude of solar perigee
      !----------------------------------------------------------------------
      sh_p1=(281.2208569+.01717836*zsy+.000047064*zday+.000001961*zhfrac)*rad
      !----------------------------------------------------------------------
      ! p Longitude of lunar perigee
      !----------------------------------------------------------------------
      sh_p=(334.3837214+40.66246584*zsy+.111404016*zday+.004641834*zhfrac)*rad

      sh_N = MOD( sh_N ,2*rpi )
      sh_s = MOD( sh_s ,2*rpi )
      sh_h = MOD( sh_h, 2*rpi )
      sh_p = MOD( sh_p, 2*rpi )
      sh_p1= MOD( sh_p1,2*rpi )

      cosI = 0.913694997 -0.035692561 *cos(sh_N)

      sh_I = ACOS( cosI )

      sin2I   = sin(sh_I)
      sh_tgn2 = tan(sh_N/2.0)

      at1=atan(1.01883*sh_tgn2)
      at2=atan(0.64412*sh_tgn2)

      sh_xi=-at1-at2+sh_N

      IF( sh_N > rpi )   sh_xi=sh_xi-2.0*rpi

      sh_nu = at1 - at2

      !----------------------------------------------------------------------
      ! For constituents l2 k1 k2
      !----------------------------------------------------------------------

      tgI2 = tan(sh_I/2.0)
      P1   = sh_p-sh_xi

      t2 = tgI2*tgI2
      t4 = t2*t2
      sh_x1ra = sqrt( 1.0-12.0*t2*cos(2.0*P1)+36.0*t4 )

      p = sin(2.0*P1)
      q = 1.0/(6.0*t2)-cos(2.0*P1)
      sh_R = atan(p/q)

      p = sin(2.0*sh_I)*sin(sh_nu)
      q = sin(2.0*sh_I)*cos(sh_nu)+0.3347
      sh_nuprim = atan(p/q)

      s2 = sin(sh_I)*sin(sh_I)
      p  = s2*sin(2.0*sh_nu)
      q  = s2*cos(2.0*sh_nu)+0.0727
      sh_nusec = 0.5*atan(p/q)
      !
   END SUBROUTINE astronomic_angle


   SUBROUTINE tide_pulse( pomega, ktide ,kc )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_pulse  ***
      !!                      
      !! ** Purpose : Compute tidal frequencies
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in ) ::   kc       ! Total number of tidal constituents
      INTEGER , DIMENSION(kc), INTENT(in ) ::   ktide    ! Indice of tidal constituents
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pomega   ! pulsation in radians/s
      !
      INTEGER  ::   jh
      REAL(wp) ::   zscale
      REAL(wp) ::   zomega_T =  13149000.0_wp
      REAL(wp) ::   zomega_s =    481267.892_wp
      REAL(wp) ::   zomega_h =     36000.76892_wp
      REAL(wp) ::   zomega_p =      4069.0322056_wp
      REAL(wp) ::   zomega_n =      1934.1423972_wp
      REAL(wp) ::   zomega_p1=         1.719175_wp
      !!----------------------------------------------------------------------
      !
      zscale =  rad / ( 36525._wp * 86400._wp ) 
      !
      DO jh = 1, kc
         pomega(jh) = (  zomega_T * Wave( ktide(jh) )%nT   &
            &          + zomega_s * Wave( ktide(jh) )%ns   &
            &          + zomega_h * Wave( ktide(jh) )%nh   &
            &          + zomega_p * Wave( ktide(jh) )%np   &
            &          + zomega_p1* Wave( ktide(jh) )%np1  ) * zscale
      END DO
      !
   END SUBROUTINE tide_pulse


   SUBROUTINE tide_vuf( pvt, put, pcor, ktide ,kc )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_vuf  ***
      !!                      
      !! ** Purpose : Compute nodal modulation corrections
      !!
      !! ** Outputs : vt: Phase of tidal potential relative to Greenwich (radians)
      !!              ut: Phase correction u due to nodal motion (radians)
      !!              ft: Nodal correction factor
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in ) ::   kc               ! Total number of tidal constituents
      INTEGER , DIMENSION(kc), INTENT(in ) ::   ktide            ! Indice of tidal constituents
      REAL(wp), DIMENSION(kc), INTENT(out) ::   pvt, put, pcor   !
      !
      INTEGER ::   jh   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      DO jh = 1, kc
         !  Phase of the tidal potential relative to the Greenwhich 
         !  meridian (e.g. the position of the fictuous celestial body). Units are radian:
         pvt(jh) = sh_T * Wave( ktide(jh) )%nT    &
            &    + sh_s * Wave( ktide(jh) )%ns    &
            &    + sh_h * Wave( ktide(jh) )%nh    &
            &    + sh_p * Wave( ktide(jh) )%np    &
            &    + sh_p1* Wave( ktide(jh) )%np1   &
            &    +        Wave( ktide(jh) )%shift * rad
         !
         !  Phase correction u due to nodal motion. Units are radian:
         put(jh) = sh_xi     * Wave( ktide(jh) )%nksi   &
            &    + sh_nu     * Wave( ktide(jh) )%nnu0   &
            &    + sh_nuprim * Wave( ktide(jh) )%nnu1   &
            &    + sh_nusec  * Wave( ktide(jh) )%nnu2   &
            &    + sh_R      * Wave( ktide(jh) )%R

         !  Nodal correction factor:
         pcor(jh) = nodal_factort( Wave( ktide(jh) )%nformula )
      END DO
      !
   END SUBROUTINE tide_vuf


   RECURSIVE FUNCTION nodal_factort( kformula ) RESULT( zf )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kformula
      !
      REAL(wp) :: zf
      REAL(wp) :: zs, zf1, zf2
      !!----------------------------------------------------------------------
      !
      SELECT CASE( kformula )
      !
      CASE( 0 )                  !==  formule 0, solar waves
         zf = 1.0
         !
      CASE( 1 )                  !==  formule 1, compound waves (78 x 78)
         zf=nodal_factort(78)
         zf = zf * zf
         !
      CASE ( 2 )                 !==  formule 2, compound waves (78 x 0)  ===  (78) 
       zf1= nodal_factort(78)
       zf = nodal_factort( 0)
       zf = zf1 * zf
       !
      CASE ( 4 )                 !==  formule 4,  compound waves (78 x 235) 
         zf1 = nodal_factort( 78)
         zf  = nodal_factort(235)
         zf  = zf1 * zf
         !
      CASE ( 5 )                 !==  formule 5,  compound waves (78 *78 x 235)
         zf1 = nodal_factort( 78)
         zf  = nodal_factort(235)
         zf  = zf * zf1 * zf1
         !
      CASE ( 6 )                 !==  formule 6,  compound waves (78 *78 x 0)
         zf1 = nodal_factort(78)
         zf  = nodal_factort( 0)
         zf  = zf * zf1 * zf1 
         !
      CASE( 7 )                  !==  formule 7, compound waves (75 x 75)
         zf = nodal_factort(75)
         zf = zf * zf
         !
      CASE( 8 )                  !==  formule 8,  compound waves (78 x 0 x 235)
         zf  = nodal_factort( 78)
         zf1 = nodal_factort(  0)
         zf2 = nodal_factort(235)
         zf  = zf * zf1 * zf2
         !
      CASE( 9 )                  !==  formule 9,  compound waves (78 x 0 x 227)
         zf  = nodal_factort( 78)
         zf1 = nodal_factort(  0)
         zf2 = nodal_factort(227)
         zf  = zf * zf1 * zf2
         !
      CASE( 10 )                 !==  formule 10,  compound waves (78 x 227)
         zf  = nodal_factort( 78)
         zf1 = nodal_factort(227)
         zf  = zf * zf1
         !
      CASE( 11 )                 !==  formule 11,  compound waves (75 x 0)
!!gm bug???? zf 2 fois !
         zf = nodal_factort(75)
         zf = nodal_factort( 0)
         zf = zf * zf1
         !
      CASE( 12 )                 !==  formule 12,  compound waves (78 x 78 x 78 x 0) 
         zf1 = nodal_factort(78)
         zf  = nodal_factort( 0)
         zf  = zf * zf1 * zf1 * zf1
         !
      CASE( 13 )                 !==  formule 13, compound waves (78 x 75)
         zf1 = nodal_factort(78)
         zf  = nodal_factort(75)
         zf  = zf * zf1
         !
      CASE( 14 )                 !==  formule 14, compound waves (235 x 0)  ===  (235)
         zf  = nodal_factort(235)
         zf1 = nodal_factort(  0)
         zf  = zf * zf1
         !
      CASE( 15 )                 !==  formule 15, compound waves (235 x 75) 
         zf  = nodal_factort(235)
         zf1 = nodal_factort( 75)
         zf  = zf * zf1
         !
      CASE( 16 )                 !==  formule 16, compound waves (78 x 0 x 0)  ===  (78)
         zf  = nodal_factort(78)
         zf1 = nodal_factort( 0)
         zf  = zf * zf1 * zf1
         !
      CASE( 17 )                 !==  formule 17,  compound waves (227 x 0) 
         zf1 = nodal_factort(227)
         zf  = nodal_factort(  0)
         zf  = zf * zf1
         !
      CASE( 18 )                 !==  formule 18,  compound waves (78 x 78 x 78 )
         zf1 = nodal_factort(78)
         zf  = zf1 * zf1 * zf1
         !
      CASE( 19 )                 !==  formule 19, compound waves (78 x 0 x 0 x 0)  ===  (78)
!!gm bug2 ==>>>   here identical to formule 16,  a third multiplication by zf1 is missing
         zf  = nodal_factort(78)
         zf1 = nodal_factort( 0)
         zf = zf * zf1 * zf1
         !
      CASE( 73 )                 !==  formule 73
         zs = sin(sh_I)
         zf = (2./3.-zs*zs)/0.5021
         !
      CASE( 74 )                 !==  formule 74
         zs = sin(sh_I)
         zf = zs * zs / 0.1578
         !
      CASE( 75 )                 !==  formule 75
         zs = cos(sh_I/2)
         zf = sin(sh_I) * zs * zs / 0.3800
         !
      CASE( 76 )                 !==  formule 76
         zf = sin(2*sh_I) / 0.7214
         !
      CASE( 77 )                 !==  formule 77
         zs = sin(sh_I/2)
         zf = sin(sh_I) * zs * zs / 0.0164
         !
      CASE( 78 )                 !==  formule 78
         zs = cos(sh_I/2)
         zf = zs * zs * zs * zs / 0.9154
         !
      CASE( 79 )                 !==  formule 79
         zs = sin(sh_I)
         zf = zs * zs / 0.1565
         !
      CASE( 144 )                !==  formule 144
         zs = sin(sh_I/2)
         zf = ( 1-10*zs*zs+15*zs*zs*zs*zs ) * cos(sh_I/2) / 0.5873
         !
      CASE( 149 )                !==  formule 149
         zs = cos(sh_I/2)
         zf = zs*zs*zs*zs*zs*zs / 0.8758
         !
      CASE( 215 )                !==  formule 215
         zs = cos(sh_I/2)
         zf = zs*zs*zs*zs / 0.9154 * sh_x1ra
         !
      CASE( 227 )                !==  formule 227 
         zs = sin(2*sh_I)
         zf = sqrt( 0.8965*zs*zs+0.6001*zs*cos (sh_nu)+0.1006 )
         !
      CASE ( 235 )               !==  formule 235 
         zs = sin(sh_I)
         zf = sqrt( 19.0444*zs*zs*zs*zs + 2.7702*zs*zs*cos(2*sh_nu) + .0981 )
         !
      END SELECT
      !
   END FUNCTION nodal_factort


   FUNCTION dayjul( kyr, kmonth, kday )
      !!----------------------------------------------------------------------
      !!  *** THIS ROUTINE COMPUTES THE JULIAN DAY (AS A REAL VARIABLE)
      !!----------------------------------------------------------------------
      INTEGER,INTENT(in) ::   kyr, kmonth, kday
      !
      INTEGER,DIMENSION(12) ::  idayt, idays
      INTEGER  ::   inc, ji
      REAL(wp) ::   dayjul, zyq
      !
      DATA idayt/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
      !!----------------------------------------------------------------------
      !
      idays(1) = 0.
      idays(2) = 31.
      inc = 0.
      zyq = MOD( kyr-1900. , 4. )
      IF( zyq == 0.)   inc = 1.
      DO ji = 3, 12
         idays(ji)=idayt(ji)+inc
      END DO
      dayjul = idays(kmonth) + kday
      !
   END FUNCTION dayjul

   !!======================================================================
END MODULE tide_mod
