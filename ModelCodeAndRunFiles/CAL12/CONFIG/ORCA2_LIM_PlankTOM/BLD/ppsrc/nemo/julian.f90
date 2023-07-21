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

MODULE julian
   !!======================================================================
   !!                       ***  MODULE julian   ***
   !! Ocean          : Julian data utilities
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   jul2greg        : Convert relative time to date
   !!   greg2jul        : Convert date to relative time
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : &       ! Precision variables
      & wp, &
      & dp  
   !USE in_out_manager           ! I/O manager
   USE lib_mpp,  ONLY : &
      & ctl_warn, ctl_stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC jul2greg,        &  ! Convert relative time to date
      &   greg2jul            ! Convert date to relative time 
  
   !! $Id: julian.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
 
   RECURSIVE SUBROUTINE jul2greg( ksec, kminut, khour, kday, kmonth, kyear, &
      &                           prelday, krefdate )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE jul2greg  ***
      !!
      !! ** Purpose : Take the relative time in days and re-express in terms of
      !!              seconds, minutes, hours, days, month, year.
      !!
      !! ** Method  : Reference date : 19500101
      !!
      !! ** Action  :
      !!
      !! History
      !!      ! 06-04  (A. Vidard) Original
      !!      ! 06-05  (A. Vidard) Reformatted and refdate      
      !!      ! 06-10  (A. Weaver) Cleanup
      !!      ! 2014-09 (D. Lea) Change to use FLOOR to deal with negative prelday
      !!-----------------------------------------------------------------------

      ! * Arguments
      INTEGER, INTENT(IN), OPTIONAL :: &
         & krefdate
      INTEGER, INTENT(OUT) :: &
         & ksec,   &
         & kminut, &
         & khour,  &
         & kday,   &
         & kmonth, &
         & kyear
      REAL(KIND=dp), INTENT(IN) :: &
         & prelday

      !! * Local declarations
      INTEGER, PARAMETER :: &
         & jpgreg = 2299161, &
         & jporef = 2433283, &
         & jparef = 2415021
      INTEGER :: &
         & ijulian, &
         & ij1,     &
         & ija,     &
         & ijb,     &
         & ijc,     &
         & ijd,     &
         & ije,     &
         & isec,    &
         & imin,    &
         & ihou,    &
         & iday,    &
         & imon,    &
         & iyea,    &
         & iref
      REAL(KIND=wp) :: &
         & zday, &
         & zref
      CHARACTER(len=200) :: &
         & cerr

      ! Main computation
      IF ( PRESENT( krefdate ) ) THEN

         SELECT CASE ( krefdate )

         CASE( 0 ) 
            iref = jpgreg

         CASE( 19500101 )
            iref = jporef

         CASE( 19000101 )
            iref = jparef

         CASE DEFAULT
            WRITE(cerr,'(A,I8.8)')'jul2greg: Unknown krefdate:', krefdate
            CALL ctl_stop( cerr )

         END SELECT

      ELSE
         iref = jporef 
      ENDIF

      zday = prelday
      ksec = FLOOR( 86400. * MOD( zday, 1. ) )

      IF ( ksec < 0. ) ksec = 86400. + ksec

      khour  = ksec / 3600
      kminut = ( ksec - 3600 * khour ) / 60
      ksec   = MOD( ksec , 60 )

      ijulian = iref + INT( zday )
      IF ( zday < 0. ) ijulian = ijulian - 1

      ! If input date after 10/15/1582 :
      IF ( ijulian >= jpgreg ) THEN
	 ij1 = INT( ( DBLE( ijulian - 1867216 ) - 0.25 ) / 36524.25 )
	 ija = ijulian + 1 + ij1 - INT( ( 0.25 * ij1 ) )
      ELSE
	 ija = ijulian
      ENDIF

      ijb = ija + 1524
      ijc = INT( 6680. + ( DBLE ( ijb - 2439870 ) - 122.1 ) / 365.25 )
      ijd = 365 * ijc + INT( 0.25 * ijc )
      ije = INT( ( ijb - ijd ) / 30.6001 )
      kday = ijb - ijd - INT( 30.6001 * ije )
      kmonth = ije - 1
      IF ( kmonth > 12 ) kmonth = kmonth - 12
      kyear = ijc - 4715
      IF ( kmonth > 2 ) kyear = kyear - 1
      IF ( kyear <= 0 ) kyear = kyear - 1

   END SUBROUTINE jul2greg

   SUBROUTINE greg2jul( ksec, kmin, khour, kday, kmonth, kyear, pjulian, &
      &                 krefdate )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE greg2jul  ***
      !!
      !! ** Purpose : Produce the time relative to the current date and time.
      !!
      !! ** Method  : The units are days, so hours and minutes transform to
      !!              fractions of a day. 
      !!
      !!              Reference date : 19500101
      !! ** Action  :
      !!
      !! History :
      !!      ! 06-04  (A. Vidard) Original
      !!      ! 06-04  (A. Vidard) Reformatted
      !!      ! 06-10  (A. Weaver) Cleanup
      !!-----------------------------------------------------------------------

      ! * Arguments
      INTEGER, INTENT(IN) :: &
         & ksec,   &
         & kmin,   &
         & khour,  & 
         & kday,   &
         & kmonth, & 
         & kyear
      REAL(KIND=dp), INTENT(OUT) :: &
         & pjulian
      INTEGER, INTENT(IN), OPTIONAL :: &
         & krefdate

      !! * Local declarations
      INTEGER, PARAMETER :: &
         & jpgreg = 15 + 31 * ( 10 + 12 * 1582 ), &     ! Gregorian calendar introduction date
         & jporef = 2433283,                      &     ! Julian reference date: 19500101
         & jparef = 2415021,                      &     ! Julian reference date: 19000101
         & jpgref = 2299161                             ! Julian reference date start of Gregorian calender
      INTEGER :: &
         & ija,     &
         & ijy,     &
         & ijm,     &
         & ijultmp, &
         & ijyear,  &
         & iref
      CHARACTER(len=200) :: &
         & cerr

      IF ( PRESENT( krefdate ) ) THEN
         SELECT CASE ( krefdate )

         CASE( 0 ) 
            iref = jpgref

         CASE( 19500101 )
            iref = jporef

         CASE( 19000101 )
            iref = jparef

         CASE DEFAULT
            WRITE(cerr,'(A,I8.8)')'greg2jul: Unknown krefdate:', krefdate
            CALL ctl_stop( cerr )

         END SELECT

      ELSE
         iref = jporef 
      ENDIF

      ! Main computation
      ijyear = kyear
      IF ( ijyear < 0 ) ijyear = ijyear + 1
      IF ( kmonth > 2 ) THEN
	 ijy = ijyear
	 ijm = kmonth + 1
      ELSE
	 ijy = ijyear  - 1
	 ijm = kmonth + 13
      ENDIF
      ijultmp = INT( 365.25 * ijy ) + INT( 30.6001 * ijm ) + kday + 1720995
      IF ( kday + 31 * ( kmonth + 12 * ijyear ) >= jpgreg ) THEN
	 ija = INT( 0.01 * ijy )
	 ijultmp = ijultmp + 2 - ija + INT( 0.25 * ija )
      ENDIF
      pjulian = ( ijultmp - iref ) + ( ( 60 * khour + kmin ) * 60 + ksec ) / 86400.

   END SUBROUTINE greg2jul
   
END MODULE julian
