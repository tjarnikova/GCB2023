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

MODULE obs_utils
   !!======================================================================
   !!                       ***  MODULE obs_utils   ***
   !! Observation diagnostics: Utility functions
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   grt_cir_dis     : Great circle distance 
   !!   grt_cir_dis_saa : Great circle distance (small angle)
   !!   chkerr          : Error-message managment for NetCDF files
   !!   chkdim          : Error-message managment for NetCDF files
   !!   fatal_error     : Fatal error handling
   !!   ddatetoymdhms   : Convert YYYYMMDD.hhmmss to components
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce, ONLY : &        ! Precision variables
      & wp, &
      & dp, &
      & i8  
   USE in_out_manager           ! I/O manager
   USE lib_mpp                  ! For ctl_warn/stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC grt_cir_dis,     &  ! Great circle distance 
      &   grt_cir_dis_saa, &  ! Great circle distance (small angle)
      &   str_c_to_for,    &  ! Remove non-printable chars from string
      &   chkerr,          &  ! Error-message managment for NetCDF files
      &   chkdim,          &  ! Check if dimensions are correct for a variable
      &   fatal_error,     &  ! Fatal error handling
      &   warning,         &  ! Warning handling
      &   ddatetoymdhms       ! Convert YYYYMMDD.hhmmss to components
         
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_utils.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
 
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: grt_cir_dis.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   REAL(KIND=wp) FUNCTION grt_cir_dis( pa1, pa2, pb1, pb2, pc1, pc2 )
      !!----------------------------------------------------------------------
      !!                     *** FUNCTION grt_cir_dis ***
      !!
      !! ** Purpose : Great circle distance between pts (lat1,lon1) 
      !!               & (lat2,lon2)
      !!                   
      !! ** Method   : Geometry.
      !!
      !! History :
      !!        !  1995-12 (G. Madec, E. Durand, A. Weaver, N. Daget) Original 
      !!        !  2006-03 (A. Vidard) Migration to NEMOVAR 
      !!        !  2006-10 (A. Weaver) Cleanup
      !!----------------------------------------------------------------------
      
      !! * Arguments
      REAL(KIND=wp) :: pa1   !  sin(lat1)
      REAL(KIND=wp) :: pa2   !  sin(lat2)
      REAL(KIND=wp) :: pb1   !  cos(lat1) * cos(lon1)
      REAL(KIND=wp) :: pb2   !  cos(lat2) * cos(lon2)
      REAL(KIND=wp) :: pc1   !  cos(lat1) * sin(lon1)
      REAL(KIND=wp) :: pc2   !  cos(lat2) * sin(lon2)

      grt_cir_dis = &
         &  ASIN( SQRT( 1.0 - ( pa1 * pa2 + pb1 * pb2 + pc1 * pc2 )**2 ) )
      
   END FUNCTION grt_cir_dis

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: grt_cir_dis_saa.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   REAL(KIND=wp) FUNCTION grt_cir_dis_saa( pa, pb, pc )
      !!----------------------------------------------------------------------
      !!                     *** FUNCTION grt_cir_dis_saa ***
      !!
      !! ** Purpose : Great circle distance between pts (lat1,lon1) 
      !!               & (lat2,lon2) with a small-angle approximation
      !!
      !! ** Method  : Geometry
      !!
      !! ** Action  :
      !!
      !! History
      !!      !  95-12 (G. Madec, E. Durand, A. Weaver, N. Daget) Original 
      !!      !  06-03 (A. Vidard) Migration to NEMOVAR 
      !!      !  06-10 (A. Weaver) Cleanup
      !!----------------------------------------------------------------------
      
      !! * Arguments
      REAL(KIND=wp) :: pa   !  lon1 - lon2
      REAL(KIND=wp) :: pb   !  lat1 - lat2
      REAL(KIND=wp) :: pc   !  cos(lat2)

      grt_cir_dis_saa = SQRT( pa * pa + ( pb * pc )**2 )

   END FUNCTION grt_cir_dis_saa
 

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: str_c_to_for.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE str_c_to_for( cd_str )
      !!---------------------------------------------------------------------
      !!   
      !!                     *** ROUTINE str_c_to_for ***
      !! 
      !! ** Purpose : Loop over a string and replace all non-printable
      !!              ASCII characters with spaces assuming English 
      !!              characters only
      !!
      !! ** Method  : Loop over a string and replace all non-printable
      !!              ASCII characters with spaces assuming English 
      !!              characters only
      !!
      !! ** Action  : 
      !! 
      !! History :  
      !!        ! : 06-05 (K. Mogensen) Original
      !!        ! : 06-05 (A. Vidard) Cleaning up
      !!        ! : 06-10 (A. Weaver) More cleaning
      !!---------------------------------------------------------------------
      !! * Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: cd_str

      !! * Local declarations
      INTEGER :: &
         & ji

      DO ji = 1, LEN( cd_str )
         IF (     ( IACHAR( cd_str(ji:ji) ) > 128 ) &
            & .OR.( IACHAR( cd_str(ji:ji) ) < 32  ) ) cd_str(ji:ji) = ' '
      END DO

   END SUBROUTINE str_c_to_for

   SUBROUTINE chkerr( kstatus, cd_name, klineno )
      !!----------------------------------------------------------------------
      !!   
      !!                    *** ROUTINE chkerr ***
      !! 
      !! ** Purpose : Error-message managment for NetCDF files.
      !! 
      !! ** Method  : 
      !!
      !! ** Action  :
      !!
      !! History
      !!      ! 02-12  (N. Daget)  hdlerr
      !!      ! 06-04  (A. Vidard) f90/nemovar migration, change name
      !!      ! 06-10  (A. Weaver) Cleanup
      !!----------------------------------------------------------------------
      !! * Modules used
      USE netcdf             ! NetCDF library
      USE dom_oce, ONLY : &  ! Ocean space and time domain variables
         & nproc

      !! * Arguments
      INTEGER :: kstatus
      INTEGER :: klineno
      CHARACTER(LEN=*) :: cd_name
      
      !! * Local declarations
      CHARACTER(len=200) :: clineno

      ! Main computation
      IF ( kstatus /= nf90_noerr ) THEN
         WRITE(clineno,'(A,I8)')' at line number ', klineno
         CALL ctl_stop( ' chkerr', ' Netcdf Error in ' // TRIM( cd_name ), &
            &           clineno, nf90_strerror( kstatus ) )
      ENDIF

   END SUBROUTINE chkerr

   SUBROUTINE chkdim( kfileid, kvarid, kndim, kdim, cd_name, klineno )
      !!----------------------------------------------------------------------
      !!   
      !!                    *** ROUTINE chkerr ***
      !! 
      !! ** Purpose : Error-message managment for NetCDF files.
      !! 
      !! ** Method  : 
      !!
      !! ** Action  :
      !!
      !! History
      !!      ! 07-03  (K. Mogenen + E. Remy) Initial version
      !!----------------------------------------------------------------------
      !! * Modules used
      USE netcdf             ! NetCDF library
      USE dom_oce, ONLY : &  ! Ocean space and time domain variables
         & nproc

      !! * Arguments
      INTEGER :: kfileid       ! NetCDF file id   
      INTEGER :: kvarid        ! NetCDF variable id   
      INTEGER :: kndim         ! Expected number of dimensions
      INTEGER, DIMENSION(kndim) :: kdim      ! Expected dimensions
      CHARACTER(LEN=*) :: cd_name            ! Calling routine name
      INTEGER ::  klineno      ! Calling line number

      !! * Local declarations
      INTEGER :: indim
      INTEGER, ALLOCATABLE, DIMENSION(:) :: &
         & idim,ilendim
      INTEGER :: ji
      LOGICAL :: llerr
      CHARACTER(len=200) :: clineno

      CALL chkerr( nf90_inquire_variable( kfileid, kvarid, ndims=indim ), &
         &         cd_name, klineno )

      ALLOCATE(idim(indim),ilendim(indim))

      CALL chkerr( nf90_inquire_variable( kfileid, kvarid, dimids=idim ), &
         &         cd_name, klineno )

      DO ji = 1, indim
         CALL chkerr( nf90_inquire_dimension( kfileid, idim(ji), &
            &                                 len=ilendim(ji) ), &
            &         cd_name, klineno )
      END DO
      
      IF ( indim /= kndim ) THEN
         WRITE(clineno,'(A,I8)')' at line number ', klineno
         CALL ctl_stop( ' chkdim',  &
            &           ' Netcdf no dim error in ' // TRIM( cd_name ), &
            &           clineno )
      ENDIF

      DO ji = 1, indim
         IF ( ilendim(ji) /= kdim(ji) ) THEN
            WRITE(clineno,'(A,I8)')' at line number ', klineno
            CALL ctl_stop( ' chkdim',  &
               &           ' Netcdf dim len error in ' // TRIM( cd_name ), &
               &           clineno )
         ENDIF
      END DO
         
      DEALLOCATE(idim,ilendim)

   END SUBROUTINE chkdim
   
  SUBROUTINE fatal_error( cd_name, klineno )
      !!----------------------------------------------------------------------
      !!
      !!                    *** ROUTINE fatal_error ***
      !!
      !! ** Purpose : Fatal error handling
      !!
      !! ** Method  :
      !!
      !! ** Action  :
      !!
      !! History
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER :: klineno
      CHARACTER(LEN=*) :: cd_name
      !! * Local declarations
      CHARACTER(len=200) :: clineno

      WRITE(clineno,'(A,I8)')' at line number ', klineno
      CALL ctl_stop( ' fatal_error', ' Error in ' // TRIM( cd_name ), &
         &           clineno)
      
   END SUBROUTINE fatal_error

   SUBROUTINE warning( cd_name, klineno )
      !!----------------------------------------------------------------------
      !!
      !!                    *** ROUTINE warning ***
      !!
      !! ** Purpose : Warning handling
      !!
      !! ** Method  :
      !!
      !! ** Action  :
      !!
      !! History
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER :: klineno
      CHARACTER(LEN=*) :: cd_name
      !! * Local declarations
      CHARACTER(len=200) :: clineno

      WRITE(clineno,'(A,I8)')' at line number ', klineno
      CALL ctl_warn( ' warning', ' Potential problem in ' // TRIM( cd_name ), &
         &           clineno)
      
   END SUBROUTINE warning

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ddatetoymdhms.h90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE ddatetoymdhms( ddate, kyea, kmon, kday, khou, kmin, ksec )
      !!----------------------------------------------------------------------
      !!
      !!                    *** ROUTINE ddatetoymdhms ***
      !!
      !! ** Purpose : Convert YYYYMMDD.hhmmss to components
      !!
      !! ** Method  :
      !!
      !! ** Action  :
      !!
      !! History
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      real(wp), INTENT(IN) :: ddate
      INTEGER, INTENT(OUT) :: kyea
      INTEGER, INTENT(OUT) :: kmon
      INTEGER, INTENT(OUT) :: kday
      INTEGER, INTENT(OUT) :: khou
      INTEGER, INTENT(OUT) :: kmin
      INTEGER, INTENT(OUT) :: ksec
      !! * Local declarations
      INTEGER :: iyymmdd
      INTEGER :: ihhmmss
      
      iyymmdd = INT( ddate )
      ihhmmss = INT( ( ddate - iyymmdd ) * 1000000 )
      kyea = iyymmdd/10000
      kmon = iyymmdd / 100 - 100 * kyea
      kday = MOD( iyymmdd, 100 )
      khou = ihhmmss/10000
      kmin = ihhmmss / 100 - 100 * khou
      ksec = MOD( ihhmmss, 100 )
      
   END SUBROUTINE ddatetoymdhms

END MODULE obs_utils
