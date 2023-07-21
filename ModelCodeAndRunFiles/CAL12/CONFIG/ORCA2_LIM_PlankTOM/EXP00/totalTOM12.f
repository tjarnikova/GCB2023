      PROGRAM totals
C
C  Calcuates totals from model output files.
C 5 versions of model output are recognised:
C        0) programs called with first argument PIT; 
C no MAC,COC,PIC,PHA,FIX,BAC
C        1) program called with version TOM5
C  no MAC,PIC,PHA,FIX,BAC
C        2) program called with version TOM5, dmspp found in dia3d file
C  no MAC,PIC,PHA,FIX,BAC
C        3) program called with version TOM10
C all 10 pft present
C        4) program c1alled with version TOM10, dmspp found in dia3d file
C all 10 pft present
C Output to two files; totals2.output contians biomass
C totals1.output contains everythig else;
C    columns of zeroes will be present if ExpCO3 or nitrfix not found 
C  in dia3d file.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT NONE
      INCLUDE 'netcdf.inc'
      INTEGER status,ncmsk,maskid
      INTEGER jpi,jpj,jpk,jpl,jpm,jpn,jpo,jpo3,stat3d
C      PARAMETER(jpi=362,jpj=292,jpk=75,jpl=15,jpm=13)
      PARAMETER(jpi=182,jpj=149,jpk=31,jpl=12,jpm=15,jpn=15,jpo=6)
      PARAMETER(jpo3=4)!number of 2d variables that are in dia3d
      INTEGER ji,jj,jk,jl,jm,itime,klev(jpo)
      REAL area(jpi,jpj),vol(jpi,jpj,jpk)
      REAL raass
      REAL*8 rat2d(jpo),conv2d(jpo),rat3d(jpn),conv3d(jpn),biotot(jpm)
      CHARACTER*50 filein,fileou
      character*4 unts(jpm),year
      character*6 sim
      character*9 version,var2d(jpo),var3d(jpn),varname(jpm)
      LOGICAL lwp,lexist
      lwp = .FALSE.
c      lwp = .TRUE.
      raass = 3600.*24.*365.
      data varname/"PRO","PTE","MES","GEL","MAC",
     &"DIA","MIX","COC","PIC","PHA","FIX","BAC","DOC","POC","GOC"/
      DATA var2d/'EXP',"ExpCO3","ExpARA","sinksil",'Cflx','Oflx'/
      DATA klev/11,11,11,11,1,1/
      DATA conv2d/12e-15,12e-15,12e-15,1e-12,12e-15,1e-12/
      DATA var3d/"PPT","GRAMIC","GRAPTE",'GRAMES',"GRAGEL",'GRAMAC',
     &'GRAMICPHY',"GRAMESPHY","denitr","nitrfix","DELO2","prococ",
     &"proara","loscal","losara"/
      DATA conv3d/12e-15,12e-15,12e-15,12e-15,12e-15,12e-15,
     &12e-15,12e-15,14e-12,14e-12,1e-12,12e-15,12e-15,
     &12e-15,12e-15/
      rat2d = 0.
      rat3d = 0.
      OPEN(UNIT=3,FILE='total.arg',STATUS='OLD')
      read( 3,*) version
      read( 3,*) sim
      read( 3,*) year
      write(*,*) version,sim, year
      do jm = 1,jpm
         unts(jm)="PgC"
      enddo
C             s/h   h/d d/y  g/mol Pg/g monthsinfile
        status = nf_open("basin_mask.nc",nf_noclobber,ncmsk)
        if (status.ne.0) STOP 'can not open basin_mask.nc'
        status = nf_inq_varid(ncmsk,'AREA',maskid)
        if (status.ne.0) write(*,*) "Problem in AREA",status
        status = nf_get_var_real(ncmsk,maskid,area)
        if (status.ne.0) write(*,*) "Problem reading AREA",status
        status = nf_inq_varid(ncmsk,'VOLUME',maskid)
        status = nf_get_var_real(ncmsk,maskid,vol)
        status = nf_close(ncmsk)
         filein = "ORCA2_1m_"//year//"0101_"//year//"1231_diad_T.nc"
        stat3d =nf_open(filein,nf_noclobber,ncmsk)
        if (stat3d.ne.0) then
          filein = "ORCA2_1m_"//year//"0101_"//year//"1231_dia3d_T.nc"
          WRITE(*,*) ' no diad, try dia3d',filein
          status =nf_open(filein,nf_noclobber,ncmsk)
          if (status.ne.0) write(*,*)  'can not open ',filein, status
c             stop
        else
          WRITE(*,*) ' opened ',filein
        endif
        status = nf_inq_varid(ncmsk,'time_counter',maskid)
        if (status.ne.0) write(*,*) "Problem in time_counter",status
        status = nf_inq_dimlen(ncmsk,maskid,itime)
          WRITE(*,*) ' set itime to 12',ncmsk,maskid,itime
        itime=12
C seconds in a year * atomic weight of Carbon * 1e-15 for Petagrams/
C   no of time intervals summed (gives PgC for one year)
      DO jm = 1, jpo3
        call sum_over_area(rat2d(jm),area,ncmsk,var2d(jm),itime,jpi,jpj,
     &    klev(jm))
      END DO
      DO jm = 1, jpn
        call sum_over_vol(rat3d(jm),vol,ncmsk,var3d(jm),itime,jpi,jpj,
     &    jpk)
      END DO
! Cflx
      IF (stat3d.ne.0) then
        status = nf_close(ncmsk)
        filein = "ORCA2_1m_"//year//"0101_"//year//"1231_dia2d_T.nc"
        stat3d =nf_open(filein,nf_noclobber,ncmsk)
      ENDIF
      if (stat3d.eq.0) THEN
       DO jm = jpo3+1, jpo
        call sum_over_area(rat2d(jm),area,ncmsk,var2d(jm),itime,jpi,jpj,
     &    klev(jm))
       END DO
      ENDIF
      status = nf_close(ncmsk)
      filein="ORCA2_1m_"//year//"0101_"//year//"1231_ptrc_T.nc"
      status =nf_open(filein,nf_noclobber,ncmsk)
       if (status.ne.0) then
         write(*,*)  'can not open ', filein
       else
         write(*,*)  'opened ', filein
       endif
        status = nf_inq_varid(ncmsk,'time_counter',maskid)
        if (status.ne.0) write(*,*) "Problem in time",status
!        status = nf_inq_dimlen(ncmsk,maskid,itime)
      do jm=1,jpm
        call sum_over_vol(biotot(jm),vol,ncmsk,
     &                 varname(jm), itime,jpi,jpj,jpk) 
        biotot(jm) = biotot(jm)/float(itime)
      end do
      status = nf_close(ncmsk)
      WRITE(*,'(A)')
     &  'year pp   gmic   gpte  gmes  ggel  gmac  gmip gmep  dnit  fnit &
     &delo2'
      WRITE(*,'(A,23(F5.1))') year,
     &  (rat3d(jm)*raass*conv3d(jm)/float(itime),jm=1,jpn)
      INQUIRE (FILE='totals1.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=2,FILE='totals1.output',FORM='FORMATTED',STATUS='OLD',
     &    POSITION='APPEND')
      ELSE
        OPEN(UNIT=2,FILE='totals1.output',FORM='FORMATTED')
        WRITE(2,'(A)')  'year  pp   gmic   gpte  gmes  ggel  gmac  gmip
     & gmep  dnit  fnit  delo2'      
        WRITE(2,'(A)')
     &  '      PgC/y PgC/y PgC/y PgC/y PgC/y PgC/y PgC/y PgC/y PgC/y
     &TgN/y TgN/y TmolO/y'
      ENDIF
      WRITE(*,'(5x,15(A6))') (varname(jm),jm=1,jpm)
      WRITE(*,'(A,12(f6.3),3(f7.3))') year,
     &  (biotot(jm)*1000.*12*1e-15,jm=1,jpm)
      INQUIRE (FILE='totals2.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=3,FILE='totals2.output',FORM='FORMATTED',STATUS='OLD',
     &    POSITION='APPEND')
      ELSE
        OPEN(UNIT=3,FILE='totals2.output',FORM='FORMATTED')
        WRITE(3,'(5x,(15(A6)))') (varname(jm),jm=1,jpm)
        WRITE(3,'(5x,(15(A6)))') (unts(jm),jm=1,jpm)
      ENDIF
      INQUIRE (FILE='totals3.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=10,FILE='totals3.output',FORM='FORMATTED',
     &    STATUS='OLD',POSITION='APPEND')
       ELSE
        OPEN(UNIT=10,FILE='totals3.output',FORM='FORMATTED')
        WRITE(10,'(A)')
     &  'year  exp   co3@100m ara@100m si@100m CO2f  O2f  '
         WRITE(10,'(A)')
     &  '      PgC/y PgC/y    PgC/y    Tmol/y  PgC/y Tmol/y'
      ENDIF
      INQUIRE (FILE='totals4.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=5,FILE='totals4.output',FORM='FORMATTED',
     &    STATUS='OLD',POSITION='APPEND')
       ELSE
        OPEN(UNIT=5,FILE='totals4.output',FORM='FORMATTED')
        WRITE(5,'(7A6)')
     &  'year',(TRIM(var3d(jm)),jm=12,jpn)  
         WRITE(5,'(A)')
     &  '      PgC/y PgC/y PgC/y PgC/y PgC/y PgC/y'
      ENDIF

      WRITE(2,'(A, 9(F6.2),2f7.1)') year,
     &  (rat3d(jm)*raass*conv3d(jm)/float(itime),jm=1,11)
      WRITE(3,'(A,12(f6.3),3(f7.3))') year,
     &  (biotot(jm)*1000.*12*1e-15,jm=1,jpm)
      write(10,'(a,(F7.2),5(F8.2))') year,
     &  (rat2d(jm)*raass*conv2d(jm)/float(itime),jm=1,jpo)
      WRITE(5,'(A, 9(F6.2))') year,
     &  (rat3d(jm)*raass*conv3d(jm)/float(itime),jm=12,jpn)
      close(5)
      close(10)
      close(4)
      close(3)
      CLOSE(2)
      END
      subroutine sum_over_area(total,area,ncmsk,var,itime,jpi,jpj,jk)
      real*8 total
      integer jpi,jpj,jk,itime
      real area(jpi,jpj),dta(jpi,jpj,itime)
      character*9 var
      integer ji,jj,jl,ncmsk,start(4),count(4)
      data start/1,1,1,1/
      start(3)=jk
      count(1)=jpi
      count(2)=jpj
      count(3)=1
      count(4)=itime
!      write(*,*) jpi,jpj,itime
       total=0.
       status = nf_inq_varid(ncmsk,var,maskid)
       IF (status .EQ. 0) THEN
        IF (jk .EQ. 1) THEN
          status = nf_get_var_real(ncmsk,maskid,dta)
        ELSE
          status = nf_get_vara_real(ncmsk,maskid,start,count,dta)
        ENDIF
        IF (status .NE. 0) STOP "sumoverarea"
        DO jl = 1, itime
            DO jj = 2,jpj-1
              DO ji = 2,jpi-1
                IF (dta(ji,jj,jl) .LT. 1. .and.dta(ji,jj,jl) .gt. -1.)   &
     &            total = total + dta(ji,jj,jl)*area(ji,jj)
              END DO
            END DO
        END DO
       ELSE
        write(*,*) "no ",var
       ENDIF 
       write(*,*) var,total*3600.*24.*365.*12.e-15/float(itime)
       return
       end
      subroutine sum_over_vol( total,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk)
      real*8 total
      integer jpi,jpj,jpk
      real vol(jpi,jpj,jpk),units,raass
      real,dimension(:,:,:,:), allocatable :: dta
      character*9 var
      integer ji,jj,jk,jl,itime,maskid
      raass = 3600.*24.*365.
      units = raass*12.*1e-15/float(itime)
      total=0.
      allocate (dta(jpi,jpj,jpk,itime))
      dta=0.
      status = nf_inq_varid(ncmsk,adjustl(trim(var)),maskid)
      if (status.ne.0) write(*,*) "Problem in",var,status
      status = nf_get_var_real(ncmsk,maskid,dta)
      DO jl = 1, itime
         DO jk = 1, jpk-1 
          DO jj = 2,jpj-1
            DO ji = 2,jpi-1
              IF (dta(ji,jj,jk,jl) .LT. 1.
     &       .and.dta(ji,jj,jk,jl) .gt. 0.  ) then
               total=total+ dta(ji,jj,jk,jl)*vol(ji,jj,jk)
              endif
            END DO
          END DO
        END DO
      END DO
      write(*,*) var,total,total*units
      deallocate (dta)
      return
      end
