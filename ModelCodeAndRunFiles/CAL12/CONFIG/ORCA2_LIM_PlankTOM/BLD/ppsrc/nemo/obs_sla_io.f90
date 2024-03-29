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

MODULE obs_sla_io
   !!======================================================================
   !!                       ***  MODULE obs_sla_io  ***
   !! Observation operators : I/O for AVISO SLA files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_enactfile    :  Read a obfbdata structure from an AVISO SLA file
   !!----------------------------------------------------------------------
   USE par_kind
   USE obs_utils
   USE obs_fbm
   use obs_sla_types
   USE netcdf

   IMPLICIT NONE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_sla_io.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obssla_io.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE read_avisofile( cdfilename, inpfile, kunit, ldwp, ldgrid )
      !!---------------------------------------------------------------------
      !!
      !!                     ** ROUTINE read_avisofile **
      !!
      !! ** Purpose : Read from file the AVISO SLA observations.
      !!
      !! ** Method  : The data file is a NetCDF file. 
      !!
      !! ** Action  :
      !!
      !! References : http://www.aviso.oceanobs.com
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
      CHARACTER(LEN=14),PARAMETER :: cpname = 'read_avisofile'
      INTEGER :: i_file_id     ! netcdf IDS
      INTEGER :: i_tracks_id
      INTEGER :: i_cycles_id
      INTEGER :: i_data_id
      INTEGER :: i_var_id
      INTEGER, PARAMETER :: imaxdim = 2    ! Assumed maximum for no. dims. in file
      INTEGER, DIMENSION(2) :: idims       ! Dimensions in file
      INTEGER :: iilen         ! Length of netCDF attributes
      INTEGER :: itype         ! Typeof netCDF attributes
      REAL(fbdp) :: zsca       ! Scale factor
      REAL(fbdp) :: zfill      ! Fill value
      CHARACTER(len=3) ::  cmission      ! Mission global attribute
      INTEGER :: itracks       ! Maximum number of passes in file
      INTEGER :: icycles       ! Maximum number of cycles for each pass
      INTEGER :: idata         ! Number of data per parameter in current file
      REAL(fbdp) :: zdeltat    ! Time gap getween two measurements in seconds
      INTEGER, DIMENSION(:), POINTER :: & 
         & iptracks,     &  ! List of passes contained in current file
         & ipnbpoints,   &  ! Number of points per pass
         & ipdataindexes    ! Index of data in theoretical profile
      INTEGER, DIMENSION(:,:), POINTER :: & 
         & ipcycles         ! List of cycles per pass
      REAL(fbdp), DIMENSION(:), POINTER :: &
         & zphi,   &        ! Latitudes
         & zlam             ! Longitudes
      REAL(fbdp), DIMENSION(:,:), POINTER :: &
         & zbegindates      ! Date of point with index 0
      REAL(fbdp) :: zbeginmiss    ! Missing data for dates
      REAL(fbsp), DIMENSION(:,:), POINTER :: &
         & zsla             ! SLA data
      REAL(fbdp), DIMENSION(:), POINTER :: &
         & zjuld            ! Julian date
      LOGICAL, DIMENSION(:), POINTER :: &
         & llskip           ! Skip observation
      CHARACTER(len=14) :: cdjuldref    ! Julian data reference
      INTEGER :: imission   ! Mission number converted from Mission global 
                            ! netCDF atttribute.
      CHARACTER(len=255) :: ctmp
      INTEGER :: iobs
      INTEGER :: jl
      INTEGER :: jm
      INTEGER :: jj
      INTEGER :: ji
      INTEGER :: jk
      INTEGER :: jobs
      INTEGER :: jcycle

      ! Open the file

      CALL chkerr( nf90_open( TRIM( cdfilename ), nf90_nowrite, i_file_id ), &
         &         cpname, 81 )

      ! Get the netCDF dimensions
      
      CALL chkerr( nf90_inq_dimid( i_file_id, 'Tracks', i_tracks_id ),  &
         &         cpname, 86 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_tracks_id, &
         &                                 len = itracks ),  &
         &         cpname, 89 )
      
      CALL chkerr( nf90_inq_dimid( i_file_id, 'Cycles', i_cycles_id ),  &
         &         cpname, 92 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_cycles_id, &
         &                                 len = icycles ),  &
         &      cpname, 95 )
      
      CALL chkerr( nf90_inq_dimid( i_file_id, 'Data', i_data_id ),  &
         &         cpname, 98 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_data_id, &
         &                                 len = idata ),  &
         &         cpname, 101 )

      ! Allocate memory for input data

      ALLOCATE( &
         & iptracks     ( itracks          ), &  
         & ipnbpoints   ( itracks          ), &  
         & ipcycles     ( icycles, itracks ), &
         & zphi         ( idata            ), &   
         & zlam         ( idata            ), &  
         & zbegindates  ( icycles, itracks ), &
         & ipdataindexes( idata            ), & 
         & zsla         ( icycles, idata   ), &
         & zjuld        ( idata*icycles    ), &
         & llskip       ( idata*icycles    )  &
         & )

      ! Get time gap getween two measurements in seconds
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'DeltaT', i_var_id ), & 
         &         cpname, 121 )
      CALL chkdim( i_file_id, i_var_id, 0, idims, cpname, 122 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zdeltat), &
         &         cpname, 124 )
      zsca = 1.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, 'scale_factor') &
         &                      == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     'scale_factor',zsca), cpname,  129)
      ENDIF
      zdeltat = zsca * zdeltat
         
      ! Get List of passes contained in current file      
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'Tracks', i_var_id ), & 
         &         cpname, 136 )
      idims(1) = itracks
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 138 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, iptracks), &
         &         cpname, 140 )
      
      ! Get Number of points per pass
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'NbPoints', i_var_id ), & 
         &         cpname, 145 )
      idims(1) = itracks
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 147 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, ipnbpoints),&
         &         cpname, 149 )
      
      ! Get list of cycles per pass
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'Cycles', i_var_id ), & 
         &         cpname, 154 )
      idims(1) = icycles
      idims(2) = itracks
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 157 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, ipcycles),   &
         &         cpname, 159 )

      ! Get longitudes

      CALL chkerr( nf90_inq_varid( i_file_id, 'Longitudes', i_var_id ), & 
         &         cpname, 164 )
      idims(1) = idata
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 166 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zlam),   &
         &         cpname, 168 )
      zsca = 1.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, 'scale_factor') &
         &                         == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     'scale_factor',zsca), cpname,  173)
      ENDIF
      zlam(:) = zsca * zlam(:)
      
      ! Get latitudes
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'Latitudes', i_var_id ), & 
         &         cpname, 180 )
      idims(1) = idata
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 182 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zphi),   &
         &         cpname, 184 )
      zsca = 1.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, 'scale_factor') &
         &                       == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     'scale_factor',zsca), cpname,  189)
      ENDIF
      zphi(:) = zsca * zphi(:)
      
      ! Get date of point with index 0
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'BeginDates', i_var_id ), & 
         &        cpname, 196 )
      idims(1) = icycles
      idims(2) = itracks
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 199 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zbegindates), &
         &         cpname, 201 )
      CALL chkerr( nf90_inquire_attribute( i_file_id, i_var_id, &
         &                                 'units', len = iilen,  &
         &                                  xtype = itype), cpname, 204 )
      IF (( itype /= NF90_CHAR ) .OR. ( iilen > 255 )) THEN
         CALL fatal_error('Error decoding BeginDates unit', 206 )
      ENDIF
      CALL chkerr( nf90_get_att( i_file_id, i_var_id, 'units', &
         &                       ctmp ), cpname, 209 )
      jl=1
      DO jm = 1, LEN(TRIM(ctmp))
         IF ((ctmp(jm:jm)>='0').AND.(ctmp(jm:jm)<='9')) THEN
            cdjuldref(jl:jl) = ctmp(jm:jm)
            jl=jl+1
         ENDIF
         IF (jl>14) EXIT
      END DO
      CALL chkerr( nf90_inquire_attribute( i_file_id, i_var_id, '_FillValue', &
         &                                  xtype = itype), cpname, 219 )
      IF ( itype /= NF90_DOUBLE ) THEN
         CALL fatal_error('Error decoding BeginDates missing data', 221 )
      ENDIF
      CALL chkerr( nf90_get_att( i_file_id, i_var_id, '_FillValue', &
         &                       zbeginmiss ), cpname, 224 )

      ! Get indices of data in theoretical profile
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'DataIndexes', i_var_id ), & 
         &         cpname, 229 )
      idims(1) = idata
      CALL chkdim( i_file_id, i_var_id, 1, idims, cpname, 231 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, ipdataindexes),   &
         &       cpname, 233 )
      
      ! Get SLA data
      
      CALL chkerr( nf90_inq_varid( i_file_id, 'SLA', i_var_id ), & 
         &         cpname, 238 )
      idims(1) = icycles
      idims(2) = idata
      CALL chkdim( i_file_id, i_var_id, 2, idims, cpname, 241 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zsla),   &
         &         cpname, 243 )
      zsca = 1.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, 'scale_factor') &
         &                        == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                       'scale_factor',zsca), cpname, 248 )
      ENDIF
      zfill = 0.0
      IF (nf90_inquire_attribute( i_file_id, i_var_id, '_FillValue') &
         &                      == nf90_noerr) THEN
         CALL chkerr( nf90_get_att( i_file_id, i_var_id, &
            &                     '_FillValue',zfill), cpname,  254 )
      ENDIF
      WHERE(zsla(:,:) /=  zfill)
         zsla(:,:) = zsca * zsla(:,:)
      ELSEWHERE
         zsla(:,:) = fbrmdi
      END WHERE
      
      ! Get the global Mission netCDF attribute
      
      cmission='   '
      CALL chkerr( nf90_inquire_attribute( i_file_id, nf90_global, &
         &                                 'Mission', len = iilen,  &
         &                                 xtype = itype), cpname, 267 )
      IF (( itype /= NF90_CHAR ) .OR. ( iilen > 3 )) THEN
         CALL fatal_error('Error decoding Mission global attribute', 269 )
      ENDIF
      CALL chkerr( nf90_get_att( i_file_id, nf90_global, 'Mission', &
         &                       cmission ), cpname, 272 )
      
      ! Convert it to an integer
      imission = 0
      DO jm = 1, imaxmissions
         IF (cmission==cmissions(jm)) THEN
            imission = jm
            EXIT
         ENDIF
      END DO
      
      ! Close the file
      
      CALL chkerr( nf90_close( i_file_id ), cpname, 285 )

      ! Compute Julian dates for all observations
      
      jm = 0
      jl = 0
      DO jj = 1, itracks
         DO ji = 1, ipnbpoints(jj)
            jl = jl + 1
            DO jk = 1, icycles
               jm = jm + 1
               IF (zbegindates(jk,jj)==zbeginmiss) THEN
                  llskip(jm) = .TRUE.
                  zjuld(jm)  = fbrmdi
               ELSE
                  llskip(jm) = .FALSE.
                  zjuld(jm)  = zbegindates(jk,jj) + &
                     &         (ipdataindexes(jl) * zdeltat / 86400._wp )
               ENDIF
            END DO
         END DO
      END DO

      ! Get rid of missing data

      jm = 0
      DO jobs = 1, idata
         DO jcycle = 1, icycles
            jm = jm + 1
            IF (zsla(jcycle,jobs) == fbrmdi) llskip(jm) = .TRUE.
         END DO
      END DO
      
      ! Allocate obfbdata

      iobs = COUNT( .NOT.llskip(:) )
      CALL init_obfbdata( inpfile )
      CALL alloc_obfbdata( inpfile, 1, iobs, 1, 0, 0, ldgrid )
      inpfile%cname(1) = 'SLA'

      ! Fill the obfbdata structure from input data

      inpfile%cdjuldref = cdjuldref
      iobs = 0
      jm = 0
      DO jobs = 1, idata
         DO jcycle = 1, icycles
            jm = jm + 1
            IF (llskip(jm)) CYCLE
            iobs = iobs + 1
            ! Characters
            WRITE(inpfile%cdwmo(iobs),'(A3,A5)') cmissions(imission), '     '
            WRITE(inpfile%cdtyp(iobs),'(I4)') imission
            ! Real values
            inpfile%plam(iobs)         = zlam(jobs)
            inpfile%pphi(iobs)         = zphi(jobs)
            inpfile%pob(1,iobs,1)      = zsla(jcycle,jobs)
            inpfile%ptim(iobs)         = zjuld(jm)
            inpfile%pdep(1,iobs)       = 0.0
            ! Integers
            inpfile%kindex(iobs)       = iobs
            inpfile%ioqc(iobs)      = 1
            inpfile%ivqc(iobs,1)    = 1
            inpfile%ivlqc(1,iobs,1) = 1
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


      ! Deallocate memory for input data

      DEALLOCATE( &
         & iptracks,      &
         & ipnbpoints,    &
         & ipcycles,      &
         & zphi,          &
         & zlam,          &  
         & zbegindates,   &
         & ipdataindexes, &
         & zsla,          &
         & zjuld,         &
         & llskip         &
         & )

   END SUBROUTINE read_avisofile


END MODULE obs_sla_io
