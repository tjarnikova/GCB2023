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

MODULE obs_sst_io
   !!======================================================================
   !!                       ***  MODULE obs_sst_io  ***
   !! Observation operators : I/O for GHRSST files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_sstfile    :  Read a obfbdata structure from an GHRSST file
   !!----------------------------------------------------------------------
   USE par_kind
   USE obs_utils
   USE obs_fbm
   USE julian
   USE netcdf

   IMPLICIT NONE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_sst_io.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obssst_io.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE read_ghrsst( cdfilename, inpfile, kunit, ldwp, ldgrid )
      !!---------------------------------------------------------------------
      !!
      !!                     ** ROUTINE read_ghrsst **
      !!
      !! ** Purpose : Read from file the GHRSST observations.
      !!
      !! ** Method  : The data file is a NetCDF file. 
      !!
      !! ** Action  :
      !!
      !! References : 
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
      CHARACTER(LEN=12),PARAMETER :: cpname = 'read_ghrsst'
      INTEGER :: i_file_id        ! netcdf IDS
      INTEGER :: i_time_id
      INTEGER :: i_ni_id
      INTEGER :: i_data_id
      INTEGER :: i_var_id
      INTEGER :: i_data           ! Number of data per parameter in current file
      INTEGER :: i_time           ! Number of reference times in file
      INTEGER, DIMENSION(:), POINTER :: &
         & i_reftime        ! Reference time in file in seconds since 1/1/1981.
      INTEGER, DIMENSION(:,:), POINTER :: &
         & i_dtime, &       ! Offset in seconds since reference time
         & i_qc,    &       ! Quality control flag.
         & i_type           ! Type of SST measurement.            
      REAL(wp), DIMENSION(:), POINTER :: &
         & z_phi,   &       ! Latitudes
         & z_lam            ! Longitudes
      REAL(wp), DIMENSION(:,:), POINTER :: &
         & z_sst            ! SST data     
      INTEGER, PARAMETER :: imaxdim = 2      ! Assumed maximum for no. dims. in file
      INTEGER, DIMENSION(2) :: idims         ! Dimensions in file
      INTEGER :: iilen            ! Length of netCDF attributes
      INTEGER :: itype            ! Typeof netCDF attributes
      REAL(KIND=wp) :: zsca           ! Scale factor
      REAL(KIND=wp) :: zoff           ! Offset for data in netcdf file
      REAL(KIND=wp) :: z_offset       ! Offset for time conversion
      REAL(KIND=wp) :: zfill          ! Fill value in netcdf file
      CHARACTER (len=33) :: creftime  ! Reference time of file
      INTEGER :: i_refyear            ! Integer version of reference time
      INTEGER :: i_refmonth
      INTEGER :: i_refday
      INTEGER :: i_refhour
      INTEGER :: i_refmin
      INTEGER :: i_refsec
      INTEGER :: ichunk
      INTEGER :: jtim
      INTEGER :: jobs
      INTEGER :: iobs

      CALL chkerr( nf90_open( TRIM( cdfilename ), nf90_nowrite, &
         &                i_file_id, chunksize=ichunk), cpname, 70 )
      
      ! Get the netCDF dimensions
      
      CALL chkerr( nf90_inq_dimid( i_file_id, 'time', i_time_id ),  &
         &       cpname, 75 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_time_id, &
         &                              len = i_time ),  &
         &      cpname, 78 )
      
      CALL chkerr( nf90_inq_dimid( i_file_id, 'ni', i_ni_id ),  &
         &       cpname, 81 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_ni_id, &
         &                              len = i_data ),  &
         &      cpname, 84 )
      
      
      ! Allocate NetCDF variables
      
      ALLOCATE( &
         & i_reftime    ( i_time                 ), &
         & i_dtime      ( i_data,i_time  ), &      
         & i_qc         ( i_data,i_time  ), & 
         & i_type       ( i_data,i_time  ), &     
         & z_phi        ( i_data                 ), &   
         & z_lam        ( i_data                 ), &  
         & z_sst        ( i_data,i_time  )  &
         & )
      
      ! Get reference time of file which is in seconds since 1981/1/1 00:00. 
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'time', i_var_id ), & 
         &       cpname, 102 )
      idims(1) = i_time
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 104 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, i_reftime),&
         &       cpname, 106 )
      IF (nf90_inquire_attribute( i_file_id, i_var_id, "units") &
         &                      == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     "units",creftime), cpname,  110 )
      ELSE
         creftime = "seconds since 1981-01-01 00:00:00"
      ENDIF
      READ(creftime(15:18),*)i_refyear
      READ(creftime(20:21),*)i_refmonth
      READ(creftime(23:24),*)i_refday
      READ(creftime(26:27),*)i_refhour
      READ(creftime(29:30),*)i_refmin
      READ(creftime(32:33),*)i_refsec
      !Work out offset in days between reference time and 1/1/1950.
      CALL greg2jul( i_refsec, i_refmin, i_refhour, i_refday, &
         &           i_refmonth, i_refyear, z_offset)
      
      ! Get list of times for each ob in seconds relative to reference time
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'sst_dtime', i_var_id ), & 
         &       cpname, 127 )
      idims(1) = i_data
      idims(2) = i_time
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 130 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, i_dtime),&
         &       cpname, 132 )
      zsca = 1.0      
      IF (nf90_inquire_attribute( i_file_id, i_var_id, "scale_factor") &
         &                      == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     "scale_factor",zsca), cpname,  137 )
      ENDIF
      zoff = 0.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, "add_offset") &
         &                      == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     "add_offset",zoff), cpname,  143 )
      ENDIF
      i_dtime(:,:) = NINT((zsca * FLOAT(i_dtime(:,:))) &
         &                 + zoff)
      
      ! Get longitudes
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'lon', i_var_id ), & 
         &         cpname, 151 )
      idims(1) = i_data
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 153 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, z_lam),   &
         &         cpname, 155 )
      
      ! Get latitudes
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'lat', i_var_id ), & 
         &         cpname, 160 )
      idims(1) = i_data
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 162 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, z_phi),   &
         &         cpname, 164 )
      
      ! Get SST data
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'sea_surface_temperature', &
         &                         i_var_id ), & 
         &         cpname, 170 )
      idims(1) = i_data
      idims(2) = i_time
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 173 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, z_sst), &
         &       cpname, 175 )
      zoff = 0.
      IF (nf90_inquire_attribute( i_file_id, i_var_id, "scale_factor") &
         &                      == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     "scale_factor",zsca), cpname, 180 )
      ENDIF
      zsca = 1.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, "scale_factor") &
         &                       == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                       "scale_factor",zsca), cpname, 186 )
      ENDIF
      zfill = 0.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, "_FillValue") &
         &                       == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                       "_FillValue",zfill), cpname, 192 )
      ENDIF
      WHERE(z_sst(:,:) /=  zfill)
         z_sst(:,:) = (zsca * z_sst(:,:)) + zoff
      ELSEWHERE
         z_sst(:,:) = fbrmdi
      END WHERE
      
      ! Get QC flag
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'confidence_flag', i_var_id ), & 
         &       cpname, 203 )
      idims(1) = i_data
      idims(2) = i_time
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 206 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, i_qc),   &
            &       cpname, 208 )
      
      ! Get SST obs type
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'data_source', i_var_id ), & 
         &         cpname, 213 )
      idims(1) = i_data
      idims(2) = i_time
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 216 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, i_type),   &
         &       cpname, 218 )
      
      ! Close the file
      
      CALL chkerr( nf90_close( i_file_id ), cpname, 222 )

      ! Fill the obfbdata structure

      ! Allocate obfbdata
      
      iobs = i_data * i_time
      CALL init_obfbdata( inpfile )
      CALL alloc_obfbdata( inpfile, 1, iobs, 1, 0, 0, ldgrid )
      inpfile%cname(1) = 'SST'

      ! Fill the obfbdata structure from input data

      inpfile%cdjuldref = "19500101000000"
      iobs = 0
      DO jtim = 1, i_time
         DO jobs = 1, i_data
            iobs = iobs + 1
            ! Characters
            WRITE(inpfile%cdwmo(iobs),'(A6,A2)') 'ghrsst','  '
            WRITE(inpfile%cdtyp(iobs),'(I4)') i_type(jobs,jtim)
            ! Real values
            inpfile%plam(iobs)         = z_lam(jobs)
            inpfile%pphi(iobs)         = z_phi(jobs)
            inpfile%pob(1,iobs,1)      = z_sst(jobs,jtim)
            inpfile%ptim(iobs)         = &
               & REAL(i_reftime(jtim))/(60.*60.*24.) + &
               & z_offset + REAL(i_dtime(jobs,jtim))/(60.*60.*24.)
            inpfile%pdep(1,iobs)       = 0.0
            ! Integers
            inpfile%kindex(iobs)       = iobs
            IF ( z_sst(jobs,jtim) == fbrmdi ) THEN
               inpfile%ioqc(iobs)      = 4
               inpfile%ivqc(iobs,1)    = 4 
               inpfile%ivlqc(1,iobs,1) = 4
            ELSE 
               inpfile%ioqc(iobs)      = i_qc(jobs,jtim)
               inpfile%ivqc(iobs,1)    = i_qc(jobs,jtim)
               inpfile%ivlqc(1,iobs,1) = 1
            ENDIF
            inpfile%ipqc(iobs)         = 0 
            inpfile%ipqcf(:,iobs)      = 0
            inpfile%itqc(iobs)         = 0
            inpfile%itqcf(:,iobs)      = 0
            inpfile%ivqcf(:,iobs,1)    = 0
            inpfile%ioqcf(:,iobs)      = 0
            inpfile%idqc(1,iobs)       = 0
            inpfile%idqcf(1,1,iobs)    = 0
            inpfile%ivlqcf(:,1,iobs,1) = 0
         END DO
      END DO

      DEALLOCATE( &
         & i_reftime, &
         & i_dtime,   &
         & i_qc,      &
         & i_type,    &
         & z_phi,     &
         & z_lam,     &
         & z_sst      &
         & )

   END SUBROUTINE read_ghrsst



END MODULE obs_sst_io
