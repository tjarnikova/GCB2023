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

MODULE obs_prof_io
   !!======================================================================
   !!                       ***  MODULE obs_prof_io  ***
   !! Observation operators : I/O for ENACT and Coriolis files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_enactfile    :  Read a obfbdata structure from an ENACT file.
   !!   read_coriofile    :  Read a obfbdata structure from an Coriolis file.
   !!----------------------------------------------------------------------
   USE par_kind
   USE obs_utils
   USE obs_fbm
   USE obs_conv
   USE netcdf
   IMPLICIT NONE

   INTEGER, PARAMETER :: imaxlev = 10000

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_prof_io.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obsprof_io.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE read_enactfile( cdfilename, inpfile, kunit, ldwp, ldgrid )
      !!---------------------------------------------------------------------
      !!
      !!                     ** ROUTINE read_enactfile **
      !!
      !! ** Purpose : Read from file the profile ENACT observations.
      !!
      !! ** Method  : The data file is a NetCDF file. 
      !!
      !! ** Action  :
      !!
      !! History : 
      !!          ! 09-01 (K. Mogensen) Original based on old versions
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(LEN=*) :: cdfilename ! Input filename
      TYPE(obfbdata)   :: inpfile    ! Output obfbdata structure
      INTEGER          :: kunit      ! Unit for output
      LOGICAL          :: ldwp       ! Print info
      LOGICAL          :: ldgrid     ! Save grid info in data structure
      !! * Local declarations
      INTEGER :: iobs                ! Number of observations
      INTEGER :: ilev                      ! Number of levels
      INTEGER :: i_file_id
      INTEGER :: i_obs_id
      INTEGER :: i_lev_id
      INTEGER :: i_phi_id
      INTEGER :: i_lam_id
      INTEGER :: i_depth_id
      INTEGER :: i_var_id
      INTEGER :: i_pl_num_id
      INTEGER :: i_reference_date_time_id
      INTEGER :: i_format_version_id
      INTEGER :: i_juld_id
      INTEGER :: i_data_type_id
      INTEGER :: i_wmo_inst_type_id
      INTEGER :: i_qc_var_id
      INTEGER :: i_dc_ref_id
      INTEGER :: i_qc_flag_id
      CHARACTER(LEN=40) :: cl_fld_lam
      CHARACTER(LEN=40) :: cl_fld_phi
      CHARACTER(LEN=40) :: cl_fld_depth 
      CHARACTER(LEN=40) :: cl_fld_var_tp
      CHARACTER(LEN=40) :: cl_fld_var_s
      CHARACTER(LEN=40) :: cl_fld_var_ti
      CHARACTER(LEN=40) :: cl_fld_var_juld_qc
      CHARACTER(LEN=40) :: cl_fld_var_pos_qc
      CHARACTER(LEN=40) :: cl_fld_var_depth_qc 
      CHARACTER(LEN=40) :: cl_fld_var_qc_t
      CHARACTER(LEN=40) :: cl_fld_var_qc_s
      CHARACTER(LEN=40) :: cl_fld_var_prof_qc_t
      CHARACTER(LEN=40) :: cl_fld_var_prof_qc_s
      CHARACTER(LEN=40) :: cl_fld_reference_date_time
      CHARACTER(LEN=40) :: cl_fld_juld
      CHARACTER(LEN=40) :: cl_fld_data_type
      CHARACTER(LEN=40) :: cl_fld_pl_num
      CHARACTER(LEN=40) :: cl_fld_format_version
      CHARACTER(LEN=40) :: cl_fld_wmo_inst_type
      CHARACTER(LEN=40) :: cl_fld_qc_flags_profiles
      CHARACTER(LEN=40) :: cl_fld_qc_flags_levels

      CHARACTER(LEN=14), PARAMETER :: cl_name = 'read_enactfile'
      CHARACTER(LEN=16)            :: cl_data_type = ''
      CHARACTER(LEN=4 )            :: cl_format_version = ''
      INTEGER, DIMENSION(1) :: istart1, icount1
      INTEGER, DIMENSION(2) :: istart2, icount2
      CHARACTER(len=imaxlev) :: clqc
      CHARACTER(len=1) :: cqc
      INTEGER :: ji, jk
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iqc1
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iqc2

      !-----------------------------------------------------------------------
      ! Initialization
      !-----------------------------------------------------------------------

      cl_fld_lam                 = 'LONGITUDE'
      cl_fld_phi                 = 'LATITUDE'
      cl_fld_depth               = 'DEPH_CORRECTED'
      cl_fld_reference_date_time = 'REFERENCE_DATE_TIME'
      cl_fld_juld                = 'JULD'
      cl_fld_data_type           = 'DATA_TYPE'
      cl_fld_format_version      = 'FORMAT_VERSION'
      cl_fld_wmo_inst_type       = 'WMO_INST_TYPE'
      cl_fld_pl_num              = 'PLATFORM_NUMBER'

      cl_fld_var_qc_t            = 'POTM_CORRECTED_QC'
      cl_fld_var_prof_qc_t       = 'PROFILE_POTM_QC'
      cl_fld_var_tp              = 'POTM_CORRECTED'
      cl_fld_var_qc_s            = 'PSAL_CORRECTED_QC'
      cl_fld_var_prof_qc_s       = 'PROFILE_PSAL_QC'
      cl_fld_var_s               = 'PSAL_CORRECTED'
      cl_fld_var_depth_qc        = 'DEPH_CORRECTED_QC'
      cl_fld_var_juld_qc         = 'JULD_QC'
      cl_fld_var_pos_qc          = 'POSITION_QC'
      cl_fld_var_ti              = 'TEMP'
      cl_fld_qc_flags_profiles   = 'QC_FLAGS_PROFILES'
      cl_fld_qc_flags_levels     = 'QC_FLAGS_LEVELS'

      icount1(1) = 1 

      !-----------------------------------------------------------------------
      ! Open file
      !-----------------------------------------------------------------------

      CALL chkerr( nf90_open( TRIM( cdfilename ), nf90_nowrite, &
            &      i_file_id ),           cl_name, 113 )

      !-----------------------------------------------------------------------
      ! Read the heading of the file
      !-----------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_data_type,      &
         &         i_data_type_id ),      cl_name, 120 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_data_type_id,        &
         &         cl_data_type ),        cl_name, 122 )
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_format_version, &
         &         i_format_version_id ), cl_name, 124 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_format_version_id,   &
         &         cl_format_version ),   cl_name, 126 )
      
      CALL str_c_to_for( cl_data_type )
      CALL str_c_to_for( cl_format_version )
      
      IF(ldwp)WRITE(kunit,*)
      IF(ldwp)WRITE(kunit,*) ' read_enactfile :' 
      IF(ldwp)WRITE(kunit,*) ' ~~~~~~~~~~~~~~~~'
      IF(ldwp)WRITE(kunit,*) '               Data type           = ', &
         &                TRIM( ADJUSTL( cl_data_type ) )
      IF(ldwp)WRITE(kunit,*) '               Format version      = ',  &
         &                TRIM( ADJUSTL( cl_format_version ) )
      
      IF ( ( ( INDEX( cl_data_type,"ENACT v1.0" ) == 1 ) .OR.   &
         &   ( INDEX( cl_data_type,"ENACT v1.4" ) == 1 ) .OR.   &
         &   ( INDEX( cl_data_type,"ENACT v1.5" ) == 1 ) .OR.   &
         &   ( INDEX( cl_data_type,"ENSEMBLES EN3 v1" ) == 1 )      ) &
         &   .AND.                                              &
         &   ( INDEX( cl_format_version,"2.0"   ) == 1 ) ) THEN
         IF(ldwp)WRITE(kunit,*)'               Valid input file'
      ELSE
         CALL fatal_error( 'Invalid input file', 147 )
      ENDIF

      !---------------------------------------------------------------------
      ! Read the number of observations and levels to allocate arrays
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_dimid        ( i_file_id, 'N_PROF', i_obs_id ),         &
         &         cl_name, 155 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_obs_id, len = iobs ),     &
         &         cl_name, 157 )
      CALL chkerr( nf90_inq_dimid        ( i_file_id, 'N_LEVELS', i_lev_id ),     &
         &         cl_name, 159 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_lev_id, len = ilev ), &
         &         cl_name, 161 )
      IF(ldwp)WRITE(kunit,*) '               No. of data records = ', iobs
      IF(ldwp)WRITE(kunit,*) '               No. of levels       = ', ilev
      IF(ldwp)WRITE(kunit,*) 
      IF (ilev > imaxlev) THEN
         CALL fatal_error( 'Increase imaxlev in obs_prof_io.F90', 166 )
      ENDIF

      !---------------------------------------------------------------------
      ! Allocate arrays
      !---------------------------------------------------------------------

      CALL init_obfbdata( inpfile )
      CALL alloc_obfbdata( inpfile, 2, iobs, ilev, 0, 1, ldgrid )
      inpfile%cname(1) = 'POTM'
      inpfile%cname(2) = 'PSAL'
      inpfile%coblong(1) = 'Potential temperature'
      inpfile%coblong(2) = 'Practical salinity'
      inpfile%cobunit(1) = 'Degrees Celsius'
      inpfile%cobunit(2) = 'PSU'
      inpfile%cextname(1) = 'TEMP'
      inpfile%cextlong(1) = 'Insitu temperature'
      inpfile%cextunit(1) = 'Degrees Celsius'

      !---------------------------------------------------------------------
      ! Read the QC atributes
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_qc_t, i_qc_var_id ),                         &
         &         cl_name, 190 )        
      istart2(1) = 1
      icount2(2) = 1
      icount2(1) = ilev
      DO ji = 1, iobs 
         istart2(2) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, clqc,                                   &
            &                         start = istart2, count = icount2),                              &
            &         cl_name, 198 )
         DO jk = 1, ilev
            inpfile%ivlqc(jk,ji,1) = IACHAR( clqc(jk:jk) ) - IACHAR( '0' )
         END DO
      END DO
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_qc_s, i_qc_var_id ),                         &
         &         cl_name, 204 )
      DO ji = 1, iobs
         istart2(2) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, clqc,                                   &
            &                         start = istart2, count = icount2),                              &
            &         cl_name, 209 )
         DO jk = 1, ilev
            inpfile%ivlqc(jk,ji,2) = IACHAR( clqc(jk:jk) ) - IACHAR( '0' )
         END DO
      END DO
      ! No depth QC in files
      DO ji = 1, iobs
         DO jk = 1, ilev
            inpfile%idqc(jk,ji)  = 1
            inpfile%idqcf(:,jk,ji) = 0
         END DO
      END DO

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_prof_qc_t, i_qc_var_id ),                    &
         &         cl_name,  223 )
      DO ji = 1,iobs
         istart1(1) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                    &
            &                         start = istart1, count = icount1),                              &
            &         cl_name, 228 )
         inpfile%ivqc(ji,1) = IACHAR( cqc ) - IACHAR( '0' )
      END DO
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_prof_qc_s, i_qc_var_id ),                    & 
         &         cl_name, 232 ) 
      DO ji = 1,iobs
         istart1(1) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                    &
            &                         start = istart1, count = icount1),                              &
            &         cl_name, 237 )
         inpfile%ivqc(ji,2) = IACHAR( cqc ) - IACHAR( '0' )
      END DO
!!      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_juld_qc, i_qc_var_id ),                       &
!!         &         cl_name, 241 ) 
!!      !DO ji = 1,iobs
!!         istart1(1) = ji
!!         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                    &
!!            &                         start = istart1, count = icount1),                              &
!!            &         cl_name, 246 )
!!         inpfile%itqc(ji)    = IACHAR( cqc ) - IACHAR( '0' )
!!         inpfile%itqcf(:,ji) = 0
!!      END DO
      ! Since the flags are not set in the ENACT files we reset them to 0
      inpfile%itqc(:)    = 1      
      inpfile%itqcf(:,:) = 0
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_pos_qc, i_qc_var_id ),                       &
         &         cl_name, 254 ) 
      DO ji = 1,iobs
         istart1(1) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                    &
            &                         start = istart1, count = icount1),                              &
            &         cl_name, 259 )
         inpfile%ipqc(ji)    = IACHAR( cqc ) - IACHAR( '0' )
         inpfile%ipqcf(:,ji) = 0
      END DO
      DO ji = 1,iobs
         inpfile%ioqc(ji) = MIN( inpfile%ivqc(ji,1), inpfile%ivqc(ji,2) )
      END DO
      IF ( nf90_inq_varid( i_file_id, cl_fld_qc_flags_profiles, i_qc_flag_id ) == nf90_noerr ) THEN
         ALLOCATE( &
            & iqc1(iobs) &
            & )
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_flag_id, iqc1 ),                                &
            &         cl_name, 271 )
         DO ji = 1,iobs
            inpfile%ioqcf(1,ji)   = iqc1(ji)
            inpfile%ivqcf(1,ji,:) = iqc1(ji)
            inpfile%ioqcf(2,ji)   = 0
            inpfile%ivqcf(2,ji,:) = 0
         END DO
         DEALLOCATE( &
            & iqc1 &
            & )
      ELSE
         IF(ldwp) WRITE(kunit,*)'No QC profile flags in file'
         inpfile%ioqcf(:,:)   = 0
         inpfile%ivqcf(:,:,:) = 0
      ENDIF
      IF ( nf90_inq_varid( i_file_id, cl_fld_qc_flags_levels, i_qc_flag_id ) == nf90_noerr ) THEN
         ALLOCATE( &
            & iqc2(ilev,iobs) &
            & )
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_flag_id, iqc2 ),                                &
            &         cl_name, 291 )
         DO ji = 1,iobs
            DO jk = 1,ilev
               inpfile%ivlqcf(1,jk,ji,:) = iqc2(jk,ji)
               inpfile%ivlqcf(2,jk,ji,:) = 0
            END DO
         END DO
         DEALLOCATE( &
            & iqc2 &
            & )
      ELSE
         IF(ldwp) WRITE(kunit,*)'No QC level flags in file'
         inpfile%ivlqcf(:,:,:,:) = 0
      ENDIF

      !---------------------------------------------------------------------
      ! Read the time/position variables 
      !---------------------------------------------------------------------
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_juld, i_juld_id ),                               &
         &         cl_name, 311 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_juld_id, inpfile%ptim ),                              &
         &         cl_name, 313 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_depth, i_depth_id ),                             &
            &         cl_name, 316 )         
      CALL chkerr( nf90_get_var  ( i_file_id, i_depth_id, inpfile%pdep ),                             &
         &         cl_name, 318 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_phi, i_phi_id ),                                 &
         &         cl_name, 321 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_phi_id, inpfile%pphi ),                               &
         &         cl_name, 323 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_lam, i_lam_id ),                                 &
         &         cl_name, 326 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_lam_id, inpfile%plam ),                               &
         &         cl_name, 328 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_reference_date_time, i_reference_date_time_id ), &
         &         cl_name, 331 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_reference_date_time_id, inpfile%cdjuldref ),          &
         &         cl_name, 333 )
      
      !---------------------------------------------------------------------
      ! Read the platform information
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_wmo_inst_type, i_wmo_inst_type_id ),             &
         &         cl_name, 340 )          
      CALL chkerr( nf90_get_var  ( i_file_id, i_wmo_inst_type_id, inpfile%cdtyp ),                    &
         &         cl_name, 342 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_pl_num, i_pl_num_id ),                           &
         &         cl_name, 345 )         
      CALL chkerr( nf90_get_var  ( i_file_id, i_pl_num_id, inpfile%cdwmo ),                           &
         &         cl_name, 347 )         

      !---------------------------------------------------------------------
      ! Read the variables
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_tp, i_var_id ),                              &
         &         cl_name, 354 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, inpfile%pob(:,:,1) ),                         &
         &         cl_name, 356 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_s, i_var_id ),                               &
         &         cl_name, 359 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, inpfile%pob(:,:,2) ),                         &
         &         cl_name, 361 )
 
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_ti, i_var_id ),                              &
         &         cl_name, 364 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, inpfile%pext(:,:,1) ),                        &
         &         cl_name, 366 )

      !---------------------------------------------------------------------
      ! Close file
      !---------------------------------------------------------------------

      CALL chkerr( nf90_close( i_file_id ),           cl_name, 372 )

      !---------------------------------------------------------------------
      ! Set file indexes
      !---------------------------------------------------------------------
      DO ji = 1, inpfile%nobs
         inpfile%kindex(ji) = ji
      END DO

   END SUBROUTINE read_enactfile

   SUBROUTINE read_coriofile( cdfilename, inpfile, kunit, ldwp, ldgrid )
      !!---------------------------------------------------------------------
      !!
      !!                     ** ROUTINE read_coriofile **
      !!
      !! ** Purpose : Read from file the profile CORIO observations.
      !!
      !! ** Method  : The data file is a NetCDF file. 
      !!
      !! ** Action  :
      !!
      !! History : 
      !!          ! 09-01 (K. Mogensen) Original based on old versions
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(LEN=*) :: cdfilename ! Input filename
      TYPE(obfbdata)   :: inpfile    ! Output enactfile structure
      INTEGER          :: kunit      ! Unit for output
      LOGICAL          :: ldwp       ! Print info
      LOGICAL          :: ldgrid     ! Save grid info in data structure
      INTEGER  :: &
         & iobs, &
         & ilev
      !! * Local declarations
      INTEGER :: &
         & i_file_id,                &
         & i_obs_id,                 &
         & i_lev_id,                 &
         & i_phi_id,                 & 
         & i_lam_id,                 &
         & i_depth_id,               &
         & i_pres_id,                &
         & i_var_id,                 &
         & i_pl_num_id,              &
         & i_format_version_id,      &
         & i_juld_id,                &
         & i_data_type_id,           &
         & i_wmo_inst_type_id,       &
         & i_qc_var_id,              &
         & i_dc_ref_id
      CHARACTER(LEN=40) :: & 
         & cl_fld_lam,                 &
         & cl_fld_phi,                 &
         & cl_fld_depth,               & 
         & cl_fld_depth_qc,            & 
         & cl_fld_pres,                & 
         & cl_fld_pres_qc,             & 
         & cl_fld_var_t,               &
         & cl_fld_var_s,               &
         & cl_fld_var_ti,              &
         & cl_fld_var_pos_qc,          &
         & cl_fld_var_qc_t,            &
         & cl_fld_var_qc_s,            &
         & cl_fld_var_prof_qc_t,       &
         & cl_fld_var_prof_qc_s,       &
         & cl_fld_dc_ref,              &
         & cl_fld_juld,                &
         & cl_fld_pl_num,              &
         & cl_fld_wmo_inst_type
      CHARACTER(LEN=14), PARAMETER :: &
         & cl_name = 'read_coriofile'
      CHARACTER(LEN=4 )            :: &
         & cl_format_version = ''
      INTEGER, DIMENSION(1) :: &
         & istart1, icount1
      INTEGER, DIMENSION(2) :: &
         & istart2, icount2
      CHARACTER(len=imaxlev) :: &
         & clqc
      CHARACTER(len=1) :: &
         & cqc
      CHARACTER(len=256) :: &
         & cdjulref
      INTEGER :: &
         & ji, jk
      INTEGER :: &
         & iformat
      LOGICAL :: & 
         & lsal
      REAL(fbdp), DIMENSION(:,:), ALLOCATABLE :: &
         & zpres
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
         & ipresqc
      CHARACTER(len=256) :: &
         & cerr
      !-----------------------------------------------------------------------
      ! Initialization
      !-----------------------------------------------------------------------

      icount1(1) = 1 
      lsal = .TRUE.

      !-----------------------------------------------------------------------
      ! Open file
      !-----------------------------------------------------------------------

      CALL chkerr( nf90_open( TRIM( cdfilename ), nf90_nowrite, &
            &      i_file_id ),           cl_name, 480 )

      !-----------------------------------------------------------------------
      ! Check format and set variables accordingly
      !-----------------------------------------------------------------------

      IF ( ( nf90_inq_dimid( i_file_id, 'N_PROF', i_obs_id ) == nf90_noerr ) .AND. &
         & ( nf90_inq_dimid( i_file_id, 'N_LEVELS', i_lev_id ) == nf90_noerr ) ) THEN
         iformat = 1
      ELSEIF ( ( nf90_inq_dimid( i_file_id, 'mN_PROF', i_obs_id ) == nf90_noerr ) .AND. &
         & ( nf90_inq_dimid( i_file_id, 'mN_ZLEV', i_lev_id ) == nf90_noerr ) ) THEN
         iformat = 2
      ELSE
         WRITE(cerr,'(2A)')'Invalid data format in ',cl_name
         CALL fatal_error( cerr, 494 )
      ENDIF
      IF ( iformat == 1 ) THEN
         cl_fld_lam                 = 'LONGITUDE'
         cl_fld_phi                 = 'LATITUDE'  
         cl_fld_depth               = 'DEPH'
         cl_fld_depth_qc            = 'DEPH_QC'
         cl_fld_pres                = 'PRES'
         cl_fld_pres_qc             = 'PRES_QC'
         cl_fld_juld                = 'JULD'
         cl_fld_wmo_inst_type       = 'WMO_INST_TYPE'
         cl_fld_dc_ref              = 'DC_REFERENCE'
         cl_fld_pl_num              = 'PLATFORM_NUMBER'
         cl_fld_var_qc_t            = 'TEMP_QC'
         cl_fld_var_prof_qc_t       = 'PROFILE_TEMP_QC'
         cl_fld_var_t               = 'TEMP'
         cl_fld_var_qc_s            = 'PSAL_QC'
         cl_fld_var_prof_qc_s       = 'PROFILE_PSAL_QC'
         cl_fld_var_s               = 'PSAL'
         cl_fld_var_pos_qc          = 'POSITION_QC'
      ELSEIF ( iformat==2 ) THEN
         cl_fld_lam                 = 'LONGITUDE'
         cl_fld_phi                 = 'LATITUDE'  
         cl_fld_depth               = 'DEPH'
         cl_fld_depth_qc            = 'QC_DEPH'
         cl_fld_pres                = 'PRES'
         cl_fld_pres_qc             = 'QC_PRES'
         cl_fld_juld                = 'JULD'
         cl_fld_wmo_inst_type       = 'INST_TYPE'
         cl_fld_dc_ref              = 'REFERENCE'
         cl_fld_pl_num              = 'PLATFORM_NUMBER'
         cl_fld_var_qc_t            = 'QC_TEMP'
         cl_fld_var_prof_qc_t       = 'Q_PROFILE_TEMP'
         cl_fld_var_t               = 'TEMP'
         cl_fld_var_qc_s            = 'QC_PSAL'
         cl_fld_var_prof_qc_s       = 'Q_PROFILE_PSAL'
         cl_fld_var_s               = 'PSAL'
         cl_fld_var_pos_qc          = 'Q_POSITION'
      ENDIF

      !-----------------------------------------------------------------------
      ! Read the heading of the file
      !-----------------------------------------------------------------------

      IF(ldwp)WRITE(kunit,*)
      IF(ldwp)WRITE(kunit,*) ' read_coriofile :' 
      IF(ldwp)WRITE(kunit,*) ' ~~~~~~~~~~~~~~~~'
      IF(ldwp)WRITE(kunit,*) '               Format version      = ', iformat

      !---------------------------------------------------------------------
      ! Read the number of observations and levels to allocate arrays
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inquire_dimension( i_file_id, i_obs_id, len = iobs ), &
         &         cl_name, 548 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_lev_id, len = ilev ), &
         &         cl_name, 550 )
      IF(ldwp)WRITE(kunit,*) '               No. of data records = ', iobs
      IF(ldwp)WRITE(kunit,*) '               No. of levels       = ', ilev
      IF(ldwp)WRITE(kunit,*) 
      IF (ilev > imaxlev) THEN
         CALL fatal_error( 'Increase imaxlev in obs_prof_io.F90', 555 )
      ENDIF

      !---------------------------------------------------------------------
      ! Allocate arrays
      !---------------------------------------------------------------------

      CALL init_obfbdata( inpfile )
      CALL alloc_obfbdata( inpfile, 2, iobs, ilev, 0, 1, ldgrid )
      inpfile%cname(1) = 'POTM'
      inpfile%cname(2) = 'PSAL'
      inpfile%coblong(1) = 'Potential temperature'
      inpfile%coblong(2) = 'Practical salinity'
      inpfile%cobunit(1) = 'Degrees Celsius'
      inpfile%cobunit(2) = 'PSU'
      inpfile%cextname(1) = 'TEMP'
      inpfile%cextlong(1) = 'Insitu temperature'
      inpfile%cextunit(1) = 'Degrees Celsius'
      ALLOCATE( &
         & zpres(ilev,iobs),  &
         & ipresqc(ilev,iobs) &
         & )
      !---------------------------------------------------------------------
      ! Get julian data reference (iformat==2)
      !---------------------------------------------------------------------

      IF (iformat==2) THEN
         CALL chkerr ( nf90_get_att( i_file_id, nf90_global, &
            &                        "Reference_date_time", cdjulref ), &
            &         cl_name, 584 )
         inpfile%cdjuldref = cdjulref(7:10)//cdjulref(4:5)// &
            & cdjulref(1:2)//cdjulref(12:13)//cdjulref(15:16)//cdjulref(18:19)
      ENDIF

      !---------------------------------------------------------------------
      ! Read the QC attributes
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_qc_t, i_qc_var_id ),                                &
         &         cl_name, 594 )
      istart2(1) = 1
      icount2(2) = 1
      icount2(1) = ilev
      DO ji = 1, iobs 
         istart2(2) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, clqc,                                   &
            &                         start = istart2, count = icount2),                              &
            &         cl_name, 602 )
         DO jk = 1, ilev
            inpfile%ivlqc(jk,ji,1) = IACHAR( clqc(jk:jk) ) - IACHAR( '0' )
         END DO
      END DO
      IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_s, i_qc_var_id ) == nf90_noerr ) THEN
         DO ji = 1, iobs
            istart2(2) = ji
            CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, clqc,                                &
               &                         start = istart2, count = icount2),                           &
               &         cl_name, 612 )
            DO jk = 1, ilev
               inpfile%ivlqc(jk,ji,2) = IACHAR( clqc(jk:jk) ) - IACHAR( '0' )
            END DO
         END DO
      ELSE
         inpfile%ivlqc(:,:,2) = 4
         inpfile%pob(:,:,2) = fbrmdi
         lsal = .FALSE.
      ENDIF

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_prof_qc_t, i_qc_var_id ),                    &
         &         cl_name,  624 )
      DO ji = 1,iobs
         istart1(1) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                    &
            &                         start = istart1, count = icount1),                              &
            &         cl_name, 629 )
         inpfile%ivqc(ji,1) = IACHAR( cqc ) - IACHAR( '0' )
      END DO
      IF (lsal) THEN
         CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_prof_qc_s, i_qc_var_id ),                 & 
            &         cl_name, 634 ) 
         DO ji = 1,iobs
            istart1(1) = ji
            CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                 &
               &                         start = istart1, count = icount1),                           &
               &         cl_name, 639 )
            inpfile%ivqc(ji,2) = IACHAR( cqc ) - IACHAR( '0' )
         END DO
      ELSE
         inpfile%ivqc(:,2) = 4
      ENDIF
      DO ji = 1,iobs
         inpfile%ioqc(ji) = MIN( inpfile%ivqc(ji,1), inpfile%ivqc(ji,2) )
      END DO
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_pos_qc, i_qc_var_id ),                       &
         &         cl_name, 649 ) 
      DO ji = 1, iobs
         istart1(1) = ji
         CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, cqc,                                    &
            &                         start = istart1, count = icount1),                              &
            &         cl_name, 654 )
         inpfile%ipqc(ji)  = IACHAR( cqc ) - IACHAR( '0' )
      END DO
      
      !---------------------------------------------------------------------
      ! Read the time/position variables 
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_juld, i_juld_id ),                               &
         &         cl_name, 663 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_juld_id, inpfile%ptim ),                              &
         &         cl_name, 665 )
      IF (iformat==1) THEN
         CALL chkerr ( nf90_get_att( i_file_id, i_juld_id, &
            &                        "units", cdjulref ), &
            &         cl_name, 669 )
         inpfile%cdjuldref = cdjulref(12:15)//cdjulref(17:18)// &
            & cdjulref(20:21)//cdjulref(23:24)//cdjulref(26:27)//&
            & cdjulref(29:30)
      ENDIF
      
      IF ( nf90_inq_varid( i_file_id, cl_fld_depth, i_depth_id ) == nf90_noerr ) THEN
         CALL chkerr( nf90_get_var  ( i_file_id, i_depth_id, inpfile%pdep ),                          &
            &         cl_name, 677 )
         CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_depth_qc, i_qc_var_id ),                      & 
            &         cl_name, 679 ) 
         DO ji = 1, iobs
            istart2(2) = ji
            CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, clqc,                                &
               &                         start = istart2, count = icount2),                           &
               &         cl_name, 684 )
            DO jk = 1, ilev
               inpfile%idqc(jk,ji) = IACHAR( clqc(jk:jk) ) - IACHAR( '0' )
            END DO
         END DO
      ELSE
         inpfile%pdep(:,:) = fbrmdi
         inpfile%idqc(:,:) = 4
      ENDIF

      IF ( nf90_inq_varid( i_file_id, cl_fld_pres, i_pres_id ) == nf90_noerr ) THEN
         CALL chkerr( nf90_get_var  ( i_file_id, i_pres_id, zpres ),                                  &
            &         cl_name, 696 )
         CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_pres_qc, i_qc_var_id ),                       & 
            &         cl_name, 698 ) 
         DO ji = 1, iobs
            istart2(2) = ji
            CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, clqc,                                &
               &                         start = istart2, count = icount2),                           &
               &         cl_name, 703 )
            DO jk = 1, ilev
               ipresqc(jk,ji) = IACHAR( clqc(jk:jk) ) - IACHAR( '0' )
            END DO
         END DO
      ELSE
         zpres(:,:) = fbrmdi
         ipresqc(:,:) = 4
      ENDIF
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_phi, i_phi_id ),                                 &
         &         cl_name, 714 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_phi_id, inpfile%pphi ),                               &
         &         cl_name, 716 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_lam, i_lam_id ),                                 &
         &         cl_name, 719 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_lam_id, inpfile%plam ),                               &
         &         cl_name, 721 )
      
      !---------------------------------------------------------------------
      ! Read the platform information
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_wmo_inst_type, i_wmo_inst_type_id ),             &
         &         cl_name, 728 )          
      CALL chkerr( nf90_get_var  ( i_file_id, i_wmo_inst_type_id, inpfile%cdtyp ),                    &
         &         cl_name, 730 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_pl_num, i_pl_num_id ),                           &
         &         cl_name, 733 )         
      CALL chkerr( nf90_get_var  ( i_file_id, i_pl_num_id, inpfile%cdwmo ),                           &
         &         cl_name, 735 )         

      
      !---------------------------------------------------------------------
      ! Read the variables
      !---------------------------------------------------------------------

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_t, i_var_id ),                               &
         &         cl_name, 743 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, inpfile%pext(:,:,1) ),                        &
         &         cl_name, 745 )

      IF (lsal) THEN      
         CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_var_s, i_var_id ),                            &
            &         cl_name, 749 )
         CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, inpfile%pob(:,:,2) ),                      &
            &         cl_name, 751 )
      ENDIF

      !---------------------------------------------------------------------
      ! Close file
      !---------------------------------------------------------------------
         
      CALL chkerr( nf90_close( i_file_id ),           cl_name, 758 )

      !---------------------------------------------------------------------
      ! Set file indexes
      !---------------------------------------------------------------------
      DO ji = 1, inpfile%nobs
         inpfile%kindex(ji) = ji
      END DO
      
      !---------------------------------------------------------------------
      !  Coriolis data conversion from insitu to potential temperature
      !---------------------------------------------------------------------
      !---------------------------------------------------------------------
      !  Convert pressure to depth if depth not present
      !---------------------------------------------------------------------
      DO ji = 1, inpfile%nobs
         IF ( inpfile%pphi(ji) < 9999.0 ) THEN
            DO jk = 1, inpfile%nlev
               IF ( inpfile%pdep(jk,ji) >= 9999.0 ) THEN
                  IF ( zpres(jk,ji) < 9999.0 ) THEN
                     inpfile%pdep(jk,ji) = &
                        & p_to_dep( REAL(zpres(jk,ji),wp), REAL(inpfile%pphi(ji),wp) )
                     inpfile%idqc(jk,ji) = ipresqc(jk,ji)
                  ENDIF
               ENDIF
            END DO
         ENDIF
      END DO
      
      !---------------------------------------------------------------------
      !  Convert depth to pressure if pressure not present
      !---------------------------------------------------------------------
      DO ji = 1, inpfile%nobs
         IF ( inpfile%pphi(ji) < 9999.0 ) THEN
            DO jk = 1, inpfile%nlev
               IF ( zpres(jk,ji) >= 9999.0 ) THEN
                  IF ( inpfile%pdep(jk,ji) < 9999.0 ) THEN
                     zpres(jk,ji) = dep_to_p( REAL(inpfile%pdep(jk,ji),wp), &
                        &                            REAL(inpfile%pphi(ji),wp) )
                     ipresqc(jk,ji) = inpfile%idqc(jk,ji)
                  ENDIF
               ENDIF
            END DO
         ENDIF
      END DO
      
      !---------------------------------------------------------------------
      !  Convert insitu temperature to potential temperature if
      !  salinity, insitu temperature and pressure are present
      !---------------------------------------------------------------------
      DO ji = 1, inpfile%nobs
         DO jk = 1, inpfile%nlev
            IF (( inpfile%pob(jk,ji,2) < 9999.0 ) .AND. &
               &( inpfile%pext(jk,ji,1) < 9999.0 ) .AND. &
               &( zpres(jk,ji) < 9999.0 ) ) THEN
               inpfile%pob(jk,ji,1) = potemp( REAL(inpfile%pob(jk,ji,2), wp),  &
                  &                           REAL(inpfile%pext(jk,ji,1), wp), &
                  &                           REAL(zpres(jk,ji),wp),  &
                  &                           0.0_wp )
            ELSE
               inpfile%pob(jk,ji,1) = fbrmdi
            ENDIF
         END DO
      END DO

      !---------------------------------------------------------------------
      !  Initialize flags since they are not in the CORIOLIS input files
      !---------------------------------------------------------------------

      inpfile%ioqcf(:,:)      = 0
      inpfile%ipqcf(:,:)      = 0
      inpfile%itqcf(:,:)      = 0
      inpfile%idqcf(:,:,:)    = 0
      inpfile%ivqcf(:,:,:)    = 0
      inpfile%ivlqcf(:,:,:,:) = 0

   END SUBROUTINE read_coriofile

END MODULE obs_prof_io
