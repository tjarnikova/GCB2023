
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
      INTEGER ncou,vardims(2),varid(4),status,ncmsk,maskid,vardim,vartim
      INTEGER jpi,jpj,jpk,jpl,lmax,jpm
      INTEGER ivn, pit,tom5,tom5dms,tom10,tom10dms
      PARAMETER(jpi=182,jpj=149,jpk=31,jpl=15,lmax=365,jpm=13)
      PARAMETER(pit=0,tom5=1,tom5dms=2,tom10=3,tom10dms=4)
      INTEGER ji,jj,jk,jl,jm,jmv,itime
    
      REAL area(jpi,jpj),vol(jpi,jpj,jpk)
      REAL ppn(jpi,jpj,jpk,jpl),ppd(jpi,jpj,jpk,jpl)
      REAL dmp(jpi,jpj,jpk,jpl),dmsflx(jpi,jpj,jpk)
      REAL oflx(jpi,jpj,lmax),cflx(jpi,jpj,lmax),rtl(lmax)
      REAL unass
      REAL*8 dmspp,dms,dmd,sumvol
      REAL*8 mortal(3),respz(3),grarem(3),grapoc(3)
      REAL*8 gramt(3),micgra(2),mesgra(2),macgra(3)
      REAL*8 expco3,sumnit,total,export,sumdelo,co3100,sil100
    
   
      REAL units,spval,raass
      REAL*8 ppna,grames,gramic,gramac,gramph(3),gramphsurf(3)
      REAL*8  sink2,resph1,mes2goc,mac2goc,goc2mes,goc2mac
      REAL*8  goc2mic,gocprod
    
      REAL*8 sumco2,sumo2,sumo,sumc
      REAL*8 biotot(jpm)
      REAL*8 sumdms,sumd
      CHARACTER*50 filein,fileou
      character*10 var
      character*3 varname(jpm)
      character*10, allocatable :: varnm(:)
      character*4, allocatable :: unts(:)
      character*4 sim,year
      character*10 version
      LOGICAL lwp,lexist
      lwp = .FALSE.
c      lwp = .TRUE.
      spval = 1.e35
      raass = 3600.*24.*365.
      data varname/"PRO","MES","MAC",
     &"DIA","MIX","COC","PIC",
     &"PHA","FIX","BAC","DOC",
     &"POC","GOC"/
      
      dms=0.
      dmd=0.
      dmspp=0.
      sumvol=0.
      sil100=0.
      co3100=0.
     
      sumnit=0.
      sumo2=0.
      OPEN(UNIT=3,FILE='total.arg',STATUS='OLD')
      read( 3,*) version
      read( 3,*) sim
      read( 3,*) year
      write(*,*) version,sim, year
      IF (trim(adjustl(version)) .EQ. "PIT") THEN
        ivn=0
       write(*,*) "PIT version ", version
        jmv=7
        allocate(varnm(7))
        allocate(unts(7))
        do jm = 1,2
         write(varnm(jm),'(A10)')varname(jm)
      
         unts(jm)="PgC"
        enddo
C No MAC
        do jm = 3,4
         write(varnm(jm),'(A10)')varname(jm+1)
       
          unts(jm)="PgC"
        enddo
        do jm = 5,7
C No COC, PIC,PHA,FIX,BAC
        
         write(varnm(jm),'(A10)')varname(jm+6)
          unts(jm)="PgC"
        enddo
      ELSEIF (trim(adjustl(version)) .EQ. "TOM5") THEN
       ivn=1
       write(*,*) "TOM5 version ", version
        jmv=8
        allocate(varnm(jmv))
        allocate(unts(jmv))
        do jm = 1,2
         write(varnm(jm),'(A10)')varname(jm)
       
         unts(jm)="PgC"
        enddo
        do jm = 3,5
C No MAC 
         write(varnm(jm),'(A10)')varname(jm+1)
       
         unts(jm)="PgC"
        enddo       
        do jm = 6,8
C No PIC,PHA,FIX,BAC
         write(varnm(jm),'(A10)')varname(jm+5)
         unts(jm)="PgC"
        enddo        
      elseif (trim(adjustl(version)) .EQ. "TOM10" .or.
     &        trim(adjustl(version)) .EQ. "PlankTOM33"  ) THEN 
       ivn=3
       write(*,*) "TOM10 version ", version
       jmv=jpm
       allocate(varnm(jmv))
       allocate(unts(jmv))
      
       do jm = 1,jmv
         write(varnm(jm),'(A10)')varname(jm)
         unts(jm)="PgC"
       enddo
       write(*,*) (varnm(jm),jm=1,jmv)
      ENDIF
C             s/h   h/d d/y  g/mol Pg/g monthsinfile
        status = nf_open("basin_mask.nc"
     &    ,nf_noclobber,ncmsk)
        if (status.ne.0) STOP 'can not open basin_mask.nc'
        status = nf_inq_varid(ncmsk,'AREA',maskid)
        if (status.ne.0) write(*,*) "Problem in AREA",status
        status = nf_get_var_real(ncmsk,maskid,area)
        if (status.ne.0) write(*,*) "Problem reading AREA",status
        status = nf_inq_varid(ncmsk,'VOLUME',maskid)
        status = nf_get_var_real(ncmsk,maskid,vol)
        status = nf_close(ncmsk)
        filein = trim(adjustl(sim))//'_'//year//"_dia3d.nc"
        status =nf_open(filein,nf_noclobber,
     &  ncmsk)
        if (status.ne.0) then
             write(*,*)  'can not open ',filein
             stop
        endif
        status = nf_inq_dimid(ncmsk,'time_counter',maskid)
        if (status.ne.0) write(*,*) "Problem in time",status
        status = nf_inq_dimlen(ncmsk,maskid,itime)
C seconds in a year * atomic weight of Carbon * 1e-15 for Petagrams/
C   no of time intervals summed (gives PgC for one year)
      units = raass*12.*1e-15/float(itime)
      ppna = 0.
      export = 0.
      gramac=0.
      grames=0.
      gramic=0.
C total grazing
      gramt = 0.
      grapoc = 0.
      grarem = 0.
      mortal=0.
      respz=0.
C grazing on phytoplankton
      gramph = 0.
      gramphsurf = 0.
C microzoo grazing on poc and bacteria
      micgra=0.
C mesozoo grazing on oc, and micr
      mesgra=0.
C macrozoo grazing on POC, GOC and zooplankton
      macgra=0.
      expco3 = 0.
C extra diagnostics for Roisin May 09
      sink2=0.
      resph1=0.
      mes2goc=0.
      mac2goc=0.
      goc2mes=0.
      goc2mac=0.
      goc2mic=0.
      gocprod=0.
      write(var,'("PPT")')
      call sum_over_vol( ppna,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
   
      status = nf_inq_varid(ncmsk,'PMO2',maskid)
      IF (status.ne.0) THEN
        write(*,*) "Problem in PMO2, try EXP",status
        status = nf_inq_varid(ncmsk,'EXP',maskid)
        if (status.ne.0) write(*,*) "Problem in EXP",status
        status = nf_get_var_real(ncmsk,maskid,ppd)
        DO jl = 1, itime
          DO jj = 2,jpj-1
            DO ji = 2,jpi-1
              IF (ppd(ji,jj,10,jl) .LT. 1) export = export+
     &          (ppd(ji,jj,10,jl))*area(ji,jj)
            END DO
          END DO
        END DO
      ELSE
        status = nf_get_var_real(ncmsk,maskid,ppd)
        status = nf_inq_varid(ncmsk,'PMO',maskid)
        if (status.ne.0)then
          write(*,*) "Problem in PMO",status
        else
         status = nf_get_var_real(ncmsk,maskid,ppn)
         DO jl = 1, itime
          DO jj = 2,jpj-1
            DO ji = 2,jpi-1
              IF (ppd(ji,jj,10,jl) .LT. 1) export = export+
     &          (ppd(ji,jj,10,jl)+ppn(ji,jj,10,jl))*area(ji,jj)
            END DO
          END DO
         END DO
        endif
      ENDIF
      status = nf_inq_varid(ncmsk,'sinksil',maskid)
        if (status.ne.0) then 
          write(*,*) "Problem in sinksil",status
        else
         status = nf_get_var_real(ncmsk,maskid,ppd)
         if (status.ne.0) then 
          write(*,*) "Problem in sinksil",status
         else
          DO jl = 1, itime
          DO jj = 2,jpj-1
            DO ji = 2,jpi-1
              IF (ppd(ji,jj,10,jl) .LT. 1) sil100 = sil100+
     &          (ppd(ji,jj,10,jl))*area(ji,jj)
            END DO
          END DO
          END DO
          sil100=sil100*raass*1.e-12/float(itime)
         endif
        endif
        status = nf_inq_varid(ncmsk,'ExpCO3',maskid)
        if (status.ne.0) then
          write(*,*) "Problem in ExpCO3",status
        else
         status = nf_get_var_real(ncmsk,maskid,ppd)
         DO jl = 1, itime
          DO jj = 2,jpj-1
            DO ji = 2,jpi-1
              IF (ppd(ji,jj,10,jl) .LT. 1) co3100 = co3100+
     &          (ppd(ji,jj,10,jl))*area(ji,jj)
            END DO
          END DO
         END DO
        endif
      write(var,'("GRAMES")')
      call sum_over_vol( grames,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
      write(var,'("GRAMIC")')
      call sum_over_vol( gramic,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)

      status = nf_inq_varid(ncmsk,'GRAMAC',maskid)
      IF (status.ne.0) THEN
        if (version .EQ. "TOM10") write(*,*) "Problem in GRAMAC",status
      ELSE
       write(var,'("GRAMAC")')
       call sum_over_vol( gramac,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
      ENDIF
      do jm =1,3
        if (jm .eq. 1 ) then
           write(var,50) 
        elseif (jm .eq. 2 ) then
          write(var,51) 
        elseif (jm .eq. 3 ) then
          write(var,52)
        endif
 50     format('GRAMICPH')
 51     format('GRAMESPH')
 52     format('GRAMACPH')
      status = nf_inq_varid(ncmsk,var,maskid)
      IF (status.ne.0) THEN     
        if (jm .eq. 1 ) then
           write(var,60) 
        elseif (jm .eq. 2 ) then
          write(var,61) 
        elseif (jm .eq. 3 ) then
          write(var,62)
        endif
 60     format('GRAMICPHY')
 61     format('GRAMESPHY')
 62     format('GRAMACPHY')

        status = nf_inq_varid(ncmsk,var,maskid)
        IF (status.ne.0) THEN
            write(*,*) "Problem in ",var, status
        ENDIF
      ENDIF
      
      call sum_over_vol( gramph(jm),vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
 
      enddo


      status = nf_inq_varid(ncmsk,'ExpCO3',maskid)
      IF (status .EQ. 0) THEN
        status = nf_get_var_real(ncmsk,maskid,ppd)
        DO jl = 1, itime
            DO jj = 2,jpj-1
              DO ji = 2,jpi-1
                IF (ppd(ji,jj,10,jl) .LT. 1) expco3 = expco3+
     &            (ppd(ji,jj,10,jl))*area(ji,jj)
              END DO
            END DO
        END DO
      ELSE
        write(*,*) "no ExpCO3"
      ENDIF   
 
      write(var,'("denitr")')
      call sum_over_vol( sumnit,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
      write(var,'("nitrfix")')
      call sum_over_vol( total,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
      sumnit=sumnit+total
      sumnit=sumnit*raass*14.*1e-12/float(itime)
      write(var,'("DEL02")')
      call sum_over_vol( sumdelo,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)

      sumdelo=sumdelo*raass*1e-12/float(itime)
 
      status = nf_inq_varid(ncmsk,'DMSPP',maskid)
      if ( status .ne. 0 ) then 
        status = nf_inq_varid(ncmsk,'dmspp',maskid) 
      endif    
      IF (status .EQ. 0) THEN
         ivn=ivn+1
         status = nf_get_var_real(ncmsk,maskid,dmp)
         DO jl = 1, itime
            DO jj = 2,jpj-1
              DO ji = 2,jpi-1
                IF (dmp(ji,jj,1,jl) .LT. 1) dmspp = dmspp+
     &            (dmp(ji,jj,1,jl))*area(ji,jj)
                sumvol=sumvol+area(ji,jj)
              END DO
            END DO
         END DO
         dmspp=dmspp*1e9/(sumvol*itime)
         sumvol=0
      ELSE
          write(*,*) "no DMSPP"
      endif 

      status = nf_inq_varid(ncmsk,'sinking2',maskid)
        if (status.ne.0) then  
         write(*,*) "Problem in sinking2",status
        else 
         status = nf_get_var_real(ncmsk,maskid,ppd)
         DO jl = 1, itime
          DO jj = 2,jpj-1
            DO ji = 2,jpi-1
              IF (ppd(ji,jj,10,jl) .LT. 1) sink2 = sink2+
     &          (ppd(ji,jj,10,jl))*area(ji,jj)
            END DO
          END DO
         END DO
       endif
c      write(var,'("sinking2")')
c      call sum_over_vol( sink2,vol,ncmsk,var,
c     &                    itime,jpi,jpj,jpl)
      write(var,'("resphy")')
      call sum_over_vol( resph1,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)    
      write(var,'("mesotoGOC")')
      call sum_over_vol( mes2goc,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
   
      write(var,'("macrtoGOC")')
      call sum_over_vol( mac2goc,vol,ncmsk,var,
     &                   itime,jpi,jpj,jpk,jpl)
   
      write(var,'("GOCtomeso")')
      call sum_over_vol( goc2mes,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
      write(var,'("GOCtomacr")')
      call sum_over_vol( goc2mac,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
  
      write(var,'("GOCtomicr")')
      call sum_over_vol( goc2mic,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
    
      write(var,'("GOCprodn")')
      call sum_over_vol( gocprod,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C open netcdf output file
C -----------------------
C
      fileou = "PISCES_flx.nc"
      status = nf_create(fileou,NF_CLOBBER,ncou)
      if(status .NE. 0) STOP ' can not open PISCES_flx.nc '
      status = nf_def_dim(ncou,'TIME',NF_UNLIMITED,vardim)
      if (status .NE. 0) then
            write(*,*) 'can not define unlimited dimesnion ',status
            stop
        endif
      status = nf_def_var(ncou,'TIME',NF_FLOAT,1,vardim,vartim)
      if (status .NE. 0) then
            write(*,*) 'can not define time ',status
            stop
        endif
      status = nf_def_var(ncou,'Cflx',nf_float,1,vardim,varid(1))
      if (status .NE. 0) then
            write(*,*) 'can not define Cflx ',status
            stop
        endif
      status = nf_put_att_text(ncou, vartim, 'units', 33, 
     &  'seconds since 1948-01-01 00:00:00')
      status = nf_put_att_text(ncou, varid(1), 'units', 4,'Pg/y')
      status = nf_put_att_real(ncou, varid(1), '_FillValue', nf_float
     &, 1,spval)
      status = nf_def_var(ncou,'Oflx',nf_float,1,vardim,varid(2))
          if (status .NE. 0) then
            write(*,*) 'can not define Oflx ',status
            stop
        endif
      status = nf_put_att_text(ncou, varid(2), 'units', 8,'TmolO2/y')
      status = nf_put_att_real(ncou, varid(2), '_FillValue', nf_float
     &, 1,spval)
      status = nf_enddef(ncou)
      status = nf_close(ncmsk)
        filein = trim(adjustl(sim))//'_'//year//"_dia2d.nc"
        status = nf_open(filein,nf_noclobber,ncmsk)
        if (status .NE. 0) then
            write(*,*) 'can not open ',filein
            stop
        endif
        status = nf_inq_varid(ncmsk,'time_counter',maskid)
        if (status.ne.0) write(*,*) "Problem in time",status
   
        status = nf_get_var_real(ncmsk,maskid,rtl)
! Cflx
        status = nf_inq_varid(ncmsk,'Cflx',maskid)
        if (status.ne.0)then
          write(*,*) "Problem in Cflx",status
        else   
         status = nf_get_var_real(ncmsk,maskid,cflx)
         sumco2=0.
         IF(lwp) write(*,*) "time,rtl ",itime,rtl
         DO jl = 1, itime
          sumc=0.
          DO jj = 2, jpj-1
            DO ji = 2, jpi-1
              IF ((area(ji,jj) .LT. 1e11) .AND. 
     &          (cflx(ji,jj,jl) .LT. 1e19)) THEN
                sumc=sumc+cflx(ji,jj,jl)*area(ji,jj)*raass*12e-15
              ENDIF
            END DO
           END DO
           sumco2=sumco2+sumc/itime
c          write(*,*) 'about to write data ',jl,vartim,rtl(jl),sumc
 
         
           status = nf_put_vara_double(ncou,varid(1),jl,1,sumc)
           IF (status .NE. 0) STOP 39
          status = nf_put_vara_real(ncou,vartim,jl,1,rtl(jl))
          if (status .NE. 0) then
            write(*,*) 'can not write time ',status
            stop 38
          endif
         END DO
        endif
! Oflx
        status = nf_inq_varid(ncmsk,'Oflx',maskid)
        if (status.ne.0) then
            write(*,*) "Problem in Oflx",status
        else
         status = nf_get_var_real(ncmsk,maskid,oflx)
         sumo2=0.
         IF(lwp) write(*,*) "time,rtl ",itime,rtl
         DO jl = 1, itime
          sumo=0.
          DO jj = 2, jpj-1
            DO ji = 2, jpi-1
              IF ((area(ji,jj) .LT. 1e11) .AND. 
     &          (oflx(ji,jj,jl) .LT. 1e19)) THEN
                sumo=sumo+oflx(ji,jj,jl)*area(ji,jj)*raass*1e-12
              ENDIF
            END DO
          END DO
          sumo2=sumo2+sumo/itime
c          IF(lwp) write(*,*) 'about to write data ',jl,rtl(jl),sumo
          status = nf_put_vara_double(ncou,varid(1),jl,1,sumo)
          IF (status .NE. 0) STOP 39
         END DO
        endif
        status = nf_inq_varid(ncmsk,'DMSFLX',maskid)
        if (status.ne.0) then
           status = nf_inq_varid(ncmsk,'dmsflx',maskid)
        endif
        if (status.eq.0) then
         status = nf_get_var_real(ncmsk,maskid,dmsflx)       
         if ( status .eq. 0 ) then 
          sumdms=0.
      
          DO jl = 1, itime
          sumd=0.
          DO jj = 2, jpj-1
            DO ji = 2, jpi-1
              IF ((area(ji,jj) .LT. 1e11) .AND. 
     &          (cflx(ji,jj,jl) .LT. 1e19)) THEN
                sumd=sumd+dmsflx(ji,jj,jl)*area(ji,jj)*raass*1e-12
              ENDIF
            END DO
          END DO
          sumdms=sumdms+sumd/itime
          END DO
         endif
        endif
      do jm =1,3
        if (jm .eq. 1 ) then
           write(var,30)     
        elseif (jm .eq. 2 ) then
          write(var,31) 
        elseif (jm .eq. 3 ) then
          write(var,32) 
        endif
 30     format('gramicph')
 31     format('gramesph')
 32     format('gramacph')
   
      status = nf_inq_varid(ncmsk,var,maskid)
      IF (status.ne.0) THEN
       if (jm .eq. 1 ) then
           write(var,40) 
        elseif (jm .eq. 2 ) then
          write(var,41) 
        elseif (jm .eq. 3 ) then
          write(var,42) 
        endif
 40     format('gramicphy')
 41     format('gramesphy')
 42     format('gramacphy')
      endif
        call sum_over_area( gramphsurf(jm),area,ncmsk,
     &      var,itime,jpi,jpj,jpl) 

      enddo
C mortality,respiration,remineralisation, Grazing to POC
      DO jm = 1,3 
      if ( jm .gt. 1 ) then     
        write(var,10) jm
 10   format( 'respz',i1)  
      else 
       write(var,20 )
 20    format('respz')
      endif
        call sum_over_area( respz(jm),area,ncmsk,
     &      var,itime,jpi,jpj,jpl) 


      END DO
      DO jm = 1,3 
      if ( jm .gt. 1 ) then     
        write(var,11) jm
 11   format( 'tortz',i1)  
      else 
       write(var,21 )
 21    format('tortz')
      endif
        call sum_over_area( mortal(jm),area,ncmsk,
     &      var,itime,jpi,jpj,jpl) 

      END DO
      DO jm = 1,3 
      if ( jm .gt. 1 ) then     
        write(var,12) jm
 12     format( 'grarem',i1)  
      else 
       write(var,22 )
 22    format('grarem')
      endif
        call sum_over_area( grarem(jm),area,ncmsk,
     &      var,itime,jpi,jpj,jpl) 

      END DO
      DO jm = 1,3 
      if ( jm .gt. 1 ) then     
        write(var,13) jm
 13     format( 'grapoc',i1)  
      else 
       write(var,23 )
 23    format('grapoc')
      endif
       call sum_over_area( grapoc(jm),area,ncmsk,
     &      var,itime,jpi,jpj,jpl) 

 
      END DO
      DO jm = 1,3 
      if ( jm .eq. 1 ) then     
        write(var,70) 
      elseif ( jm .eq. 2 ) then     
        write(var,71) 
      elseif ( jm .eq. 3 ) then     
        write(var,72) 
      endif 
 70     format( 'gramit')  
 71     format( 'gramet')  
 72     format( 'gramat')  
       call sum_over_area( gramt(jm),area,ncmsk,var,itime,jpi,jpj,jpl) 

      END DO
      DO jm = 1,2 
      if ( jm .eq. 1 ) then     
        write(var,80) 
      elseif ( jm .eq. 2 ) then     
        write(var,81) 

      endif 
 80     format( 'grazm')  
 81     format( 'grabac')  

       call sum_over_area( micgra(jm),area,ncmsk,var,itime,jpi,jpj,jpl) 

      END DO
      DO jm = 1,2 
      if ( jm .eq. 1 ) then     
        write(var,90) 
      elseif ( jm .eq. 2 ) then     
        write(var,91) 

      endif 
 90     format( 'grazoc')  
 91     format( 'grazz')  

       call sum_over_area( mesgra(jm),area,ncmsk,var,itime,jpi,jpj,jpl) 

      END DO
      DO jm = 1,3 
      if ( jm .eq. 1 ) then     
        write(var,100) 
      elseif ( jm .eq. 2 ) then     
        write(var,101) 
      elseif ( jm .eq. 3 ) then     
        write(var,102)
      endif 
 100     format( 'grampoc')  
 101     format( 'gramgoc')
 102     format( 'gramaczoo')

       call sum_over_area( macgra(jm),area,ncmsk,var,itime,jpi,jpj,jpl) 

      END DO
      status = nf_close(ncmsk)
      filein=trim(adjustl(sim))//'_'//year//"_ptrc.nc"
      status =nf_open(filein,nf_noclobber,ncmsk)
       if (status.ne.0) then
         write(*,*)  'can not open ', filein
       endif
        status = nf_inq_dimid(ncmsk,'time_counter',maskid)
        if (status.ne.0) write(*,*) "Problem in time",status
        status = nf_inq_dimlen(ncmsk,maskid,itime)
      do jm=1,jmv
        biotot(jm) = 0.
!        write(*,*) varnm(jm)
        call sum_over_vol(biotot(jm),vol,ncmsk,
     &                 varnm(jm), itime,jpi,jpj,jpk,jpl) 
        biotot(jm) = biotot(jm)/float(itime)
 
      end do
      status = nf_inq_varid(ncmsk,'DMS',maskid)
      if ( status .ne. 0 ) then 
        status = nf_inq_varid(ncmsk,'dms',maskid)
      endif
      IF (status .EQ. 0) THEN
         status = nf_get_var_real(ncmsk,maskid,dmp)
         DO jl = 1, itime
            DO jj = 2,jpj-1
              DO ji = 2,jpi-1
                IF (dmp(ji,jj,1,jl) .LT. 1) dms = dms+
     &            (dmp(ji,jj,1,jl))*area(ji,jj)
                 sumvol=sumvol+area(ji,jj)
              END DO
            END DO
         END DO
         dms=dms*1e9/(sumvol*itime)
         sumvol=0
      ELSE
        write(*,*) "no DMS"
      endif
     
      status = nf_inq_varid(ncmsk,'DMD',maskid)
      IF (status .NE. 0) THEN
       status = nf_inq_varid(ncmsk,'dmd',maskid)
      ENDIF
      IF (status .EQ. 0) THEN
        status = nf_get_var_real(ncmsk,maskid,dmp)
        DO jl = 1, itime
            DO jj = 2,jpj-1
              DO ji = 2,jpi-1
                IF (dmp(ji,jj,1,jl) .LT. 1) dmd = dmd+
     &            (dmp(ji,jj,1,jl))*vol(ji,jj,1)
                  sumvol=sumvol+vol(ji,jj,1)
              END DO
            END DO
        END DO
        dmd=dmd*1e9/(sumvol*itime)
      ELSE
        write(*,*) "no DMD"
      endif
      status = nf_close(ncmsk)
      status = nf_close(ncou)
      write(*,*) "version number is",ivn
      WRITE(*,'(A54,13(A5))')
     &  'year pp   exp  CO2f gmac gmes gmic gmap gmep gmip co3x nitr ',
     &  (varname(jm),jm=1,jmv)
      WRITE(*,'(A,9(F5.1),e10.3,13(f5.2))') year,
     &  (ppna)*units,export*units,sumco2,gramac*units,grames*units,
     &  gramic*units,gramph(3)*units,gramph(2)*units,gramph(1)*units,
     &  expco3*units,sumnit,(biotot(jm)*12*1e-12,jm=1,jmv)
      INQUIRE (FILE='totals1.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=2,FILE='totals1.output',FORM='FORMATTED',STATUS='OLD',
     &    POSITION='APPEND')
      ELSE
        OPEN(UNIT=2,FILE='totals1.output',FORM='FORMATTED')
        if ( ivn .eq. tom10dms ) then 
          WRITE(2,'(A)')
     &  'year  pp    exp  CO2f  gmic   gmes  gmac  gmip  gmep  gmap  co3
     &x  nitr dmsflx delo2 O2flx'
          WRITE(2,'(5x, A)')
     &    ' PgC   PgC  PgC/y PgC    PgC   PgC   PgC   PgC   PgC   PgC        
     &TgN   Tmol   Tmol TmolO/y'

        elseif ( ivn .eq. tom5dms ) then 
          WRITE(2,'(A)')
     &  'year  pp    exp  CO2f  gmic  gmes  gmip  gmep  co3x  dmsf  O2f
     &'
          WRITE(2,'(5x, A)') ' PgC   PgC  PgC/y PgC   PgC   PgC   PgC   
     &  PgC   Tmol  TmolO/y    '
        elseif   ( ivn .eq. tom5 ) then
          WRITE(2,'(A)')
     &  'year  pp    exp  CO2f  gmic  gmes  gmip  gmep  co3x O2f'
          WRITE(2,'(5x, A)') ' PgC   PgC  PgC/y  PgC   PgC   PgC   PgC 
     &   PgC  TmolO/y'
        elseif   ( ivn .eq. tom10 ) then
          WRITE(2,'(A)')
     &  'year  pp    exp  CO2f  gmic   gmes  gmac  gmip  gmep  gmap co3x
     & nitr  delo2    O2f'      
          WRITE(2,'(5x, A)')
     &    ' PgC   PgC  PgC/y PgC    PgC   PgC   PgC   PgC   PgC   PgC         
     &TgN   TmolO  TmolO/y'
        endif
      ENDIF
      INQUIRE (FILE='totals2.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=3,FILE='totals2.output',FORM='FORMATTED',STATUS='OLD',
     &    POSITION='APPEND')
      ELSE
        OPEN(UNIT=3,FILE='totals2.output',FORM='FORMATTED')
        if ( ivn .eq. pit .or. ivn .eq. tom5 ) then
        WRITE(3,'(6x,(13(A6)))')
     &    (varnm(jm),jm=1,jmv)
        WRITE(3,'(5x,(13(A6)))')
     &    (unts(jm),jm=1,jmv)
        elseif ( ivn .eq. tom5dms ) then
        WRITE(3,'(5x,5(1x,A5),3(2x,A5),A)')
     &    (trim(adjustl(varnm(jm))),jm=1,jmv),"    DMS    DMD   DMSPP"
        WRITE(3,'(5x,5(1x,A5),3(2x,A5),A)'),
     &    (trim(adjustl(unts(jm))),jm=1,jmv), "  Gmol/l Gmol/l Gmol/l "
 
        else  if ( ivn .eq. tom10dms ) then 
         write(*,*) 'ivn is tom10dms'
c        WRITE(3,'10(1x,A5),3(2x,A5)')
c     &   " year ",(varnm(jm),jm=1,jmv),"    DMS   DMD   DMSPP"
c        WRITE(3,'5x,10(1x,A5),3(2x,A5)'),
c     &   (unts(jm),jm=1,jmv),"   Gmol/l Gmol/l Gmol/l "

        write(*,*) jmv, (varnm(jm),jm=1,jmv)
        WRITE(3,'(A6,10(1x,A5),3(2x,A5),A21)'),
     &   " year ", 
     &   (trim(adjustl(varnm(jm))),jm=1,jmv),"    DMS   DMD   DMSPP"
        WRITE(3,'(A6,10(1x,A5),3(2x,A5),A21)'),
     &   "     ",(unts(jm),jm=1,jmv),"   Gmol/l Gmol/l Gmol/l "
        else 

           WRITE(3,'(7x,(10(1x,A5)),3(2x,A5))'), (varnm(jm),jm=1,jmv)
           WRITE(3,'(6x,(10(1x,A5)),3(2x,A5))'), (unts(jm),jm=1,jmv)
        endif
      ENDIF
      INQUIRE (FILE='totals3.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=4,FILE='totals3.output',FORM='FORMATTED',STATUS='OLD',
     &    POSITION='APPEND')
      ELSE
        OPEN(UNIT=4,FILE='totals3.output',FORM='FORMATTED')
        if ( ivn .eq. tom10 .or. ivn .eq. tom10dms ) then 
          WRITE(4,'(A)')
     &  'year     resp to DOC     |     mortality      |  grazing to POC
     &    |  remineralisation'
          WRITE(4,'(A)')
     &  '       pro   mes   mac   |  pro   mes   mac   |  pro   mes   ma 
     &c   |  pro   mes   mac'
         WRITE(4,'( A)')
     &  '       PgC   PgC   PgC      PgC   PgC   PgC   |  PgC   PgC   Pg
     &C   |  PgC   PgC   PgC'
         endif
       ENDIF
       INQUIRE (FILE='totals4.output',EXIST=lexist)
      IF(lexist) THEN
        OPEN(UNIT=9,FILE='totals4.output',FORM='FORMATTED',STATUS='OLD',
     &    POSITION='APPEND')
      ELSE
        OPEN(UNIT=9,FILE='totals4.output',FORM='FORMATTED')
        if ( ivn .eq. tom10 .or. ivn .eq. tom10dms ) then 
          WRITE(9,'(A)')
     &  'year     total grazing   |  grazing on phyto  |  pro        |  
     &mes         |  mac '
          WRITE(9,'(A)')
     &  '       pro   mes   mac   |  pro   mes   mac   |  poc  bac   |    
     &oc   pro    |  poc   goc   zoo'
         WRITE(9,'(A)')
     &  '       PgC   PgC   PgC   |  PgC   PgC   PgC   |  PgC  PgC   |  
     &PgC  PgC    |  PgC   PgC   PgC '
         endif
       ENDIF
  
       if ( ivn .eq. tom10 .or. ivn .eq. tom10dms ) then 
         INQUIRE (FILE='totals5.output',EXIST=lexist)
         IF(lexist) THEN
           OPEN(UNIT=10,FILE='totals5.output',FORM='FORMATTED',
     &    STATUS='OLD',POSITION='APPEND')
          ELSE
            OPEN(UNIT=10,FILE='totals5.output',FORM='FORMATTED')
       
            WRITE(10,'(A)')
     &  'year  si@100m  co3@100m  GOC sink resp(dia+N2fix) mes2GOC mac2GOC 
     &GOC2mes GOC2mac GOC2mic GOCprod'
            WRITE(10,'(A)')
     &  '        TmolSi     PgC?     PgC         PgC          PgC     PgC  
     &   PgC     PgC   PgC     PgC'
          ENDIF
        endif

      if ( ivn .eq. tom5 .or. ivn .eq. pit  ) then 
        WRITE(2,'(A,7(F6.2),f6.2,f7.1)') year,
     &  (ppna)*units,export*units,sumco2,gramic*units,
     &  grames*units,gramph(1)*units,gramph(2)*units,
     &  expco3*units,sumo2
      elseif ( ivn .eq. tom5dms   ) then 
        WRITE(2,'(A,8(F6.2),f6.2,f7.1)') year,
     &  (ppna)*units,export*units,sumco2,gramic*units,
     &  grames*units,gramph(1)*units,gramph(2)*units,
     &  expco3*units,sumdms,sumo2
      elseif (  ivn .eq. tom10dms  ) then 
        WRITE(2,'(A,10(F6.2),f6.0,f6.2,1x,f6.0, f7.1)') year,
     &  (ppna)*units,export*units,sumco2,gramic*units,grames*units,
     &  gramac*units,gramph(1)*units,gramph(2)*units,gramph(3)*units,
     &  expco3*units,sumnit, sumdms,sumdelo,sumo2
      else 
        WRITE(2,'(A, 10(F6.2),f6.0,f6.0,f7.1)') year,
     &  (ppna)*units,export*units,sumco2,gramic*units,grames*units,
     &  gramac*units,gramph(1)*units,gramph(2)*units,gramph(3)*units,
     &  expco3*units,sumnit,sumdelo,sumo2
      endif
      if ( ivn .eq. pit .or. ivn .eq. tom5 ) then
          WRITE(3,'(A,2x,13(f6.3))') year,
     &           (biotot(jm)*12*1e-12,jm=1,jmv)
      elseif ( ivn .eq. tom5dms ) then
          WRITE(3,'(A,2x,5(f6.3),6(f7.3))') year,
     &       (biotot(jm)*12*1e-12,jm=1,8),
     &       dms,dmd,dmspp
      elseif ( ivn .eq. tom10dms ) then
          WRITE(3,'(A,2x,10(f6.3),6(f7.3))') year,
     &       (biotot(jm)*12*1e-12,jm=1,jmv),
     &       dms,dmd,dmspp
      else
         WRITE(3,'(A,2x,10(f6.3),3(f7.3))') year,
     &       (biotot(jm)*12*1e-12,jm=1,jmv)
      endif
      if ( ivn .eq. tom10 .or. ivn .eq. tom10dms ) then 
       write(4,'(a,1x,4(3f6.3,3x))') year, (units*respz(jm),jm=1,3),
     &        (units*mortal(jm),jm=1,3),
     &     (units*grapoc(jm),jm=1,3), (units*grarem(jm),jm=1,3)
      endif 
      if ( ivn .eq. tom10 .or. ivn .eq. tom10dms ) then 
       write(9,'(a,1x,2(3f6.3,3x),2(2f6.3,3x),3f6.3,3x)') 
     &          year, (units*gramt(jm),jm=1,3),
     &   (units*gramphsurf(jm),jm=1,3), (units*micgra(jm),jm=1,2),
     &  (units*mesgra(jm),jm=1,2), (units*macgra(jm),jm=1,3)
       write(10,'(a,3x,f6.3,3x,f6.3,3x,f6.3,5x,f6.3,8x,6(f6.3,1x))') 
     &          year, sil100,co3100*units,
     &         units*sink2, units*resph1, units*mes2goc,
     &          units*mac2goc, units*goc2mes, units*goc2mac, 
     &          units*goc2mic, units*gocprod
      endif 
      close(9)
      close(4)
      close(3)
      CLOSE(2)
 500  continue
      END
      subroutine sum_over_area( total,area,ncmsk,var,itime,jpi,jpj,jpl)
      real*8 total
      
      real area(jpi,jpj),dta(jpi,jpj,itime)
      character*10 var
      integer ji,jj,jl,itime,ncmsk
!      write(*,*) jpi,jpj,jpl,itime
       total=0.
       write(*,*) 'find variable ', var
       status = nf_inq_varid(ncmsk,var,maskid)
       IF (status .EQ. 0) THEN
        status = nf_get_var_real(ncmsk,maskid,dta)
        DO jl = 1, itime
            DO jj = 2,jpj-1
              DO ji = 2,jpi-1
                IF (dta(ji,jj,jl) .LT. 1.) total = total +
     &            dta(ji,jj,jl)*area(ji,jj)
              END DO
            END DO
        END DO
       ELSE
        write(*,*) "no ",var
     
       ENDIF 
!      write(*,*) total
       return
       end
      subroutine sum_over_vol( total,vol,ncmsk,var,
     &                    itime,jpi,jpj,jpk,jpl)
      real*8 total
      integer jpi,jpj,jpk,jpl
      real vol(jpi,jpj,jpk)
      real,dimension(:,:,:,:), allocatable :: dta
      character*10 var
      integer ji,jj,jk,jl,itime,maskid
      total=0.
      write(*,*)jpi,jpj,jpk,itime
      write(*,*) 'allocating ', jpi*jpk*jpj*itime*8
      allocate (dta(jpi,jpj,jpk,itime))
      dta=0.
      write(*,*) 'Looking for ', var
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
!      write(*,*) total
      deallocate (dta)
      return
      end
