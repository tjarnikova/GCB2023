       SUBROUTINE bgcatm
#if defined key_planktom && key_trc_atmco2 

!
!-----------------------------------------------------------------------
!
!      Calculates the atmospheric CO2 concentration
!
!-----------------------------------------------------------------------
!
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      IMPLICIT NONE
!
      INTEGER iyy,is
      REAL zyr,slope,inter
!
      iyy = ndastp/10000
      zyr = float(iyy) + (float(nday_year)-0.5)/365.
!
      IF (zyr .lt. yrco2(1) .OR. zyr .gt. yrco2(nmaxrec)) THEN
        write(numout,*)'Caution: The date is outside tabulated values.'
        write(numout,*)'zyr, zyr(min) zyr(max) = ', zyr, yrco2(1),yrco2(nmaxrec)
        write(numout,*)'approximate pCO2'
        pco2at = 278+0.00009*max(zyr-1870.,0.)**2.847
#  if defined key_trc_n2o
        pn2oat = rn_atmn2o
#  endif
#  if defined key_trc_ch4
        CALL FLUSH(numout)
        STOP 'bgcatm pCH4'
#  endif
      ENDIF
      DO is = 1, nmaxrec-1
        IF(zyr.ge.yrco2(is) .AND. zyr.le.yrco2(is+1)) THEN
          slope = (sipco2(is+1) - sipco2(is)) /(yrco2(is+1) - yrco2(is))
          inter = sipco2(is)-slope*yrco2(is)
          pco2at = slope*zyr + inter
#  if defined key_trc_n2o
          slope = (sipn2o(is+1) - sipn2o(is)) /(yrco2(is+1) - yrco2(is))
          inter = sipn2o(is)-slope*yrco2(is)
          pn2oat = slope*zyr + inter
#  endif
#  if defined key_trc_ch4
          slope = (sipch4(is+1) - sipch4(is)) /(yrco2(is+1) - yrco2(is))
          inter = sipch4(is)-slope*yrco2(is)
          pch4at = slope*zyr + inter
#  endif
          EXIT
        ENDIF
      END DO
      IF(lwp) WRITE(numout,*) 'Atmospheric xCO2 xN2O xCH4 concentration: ',&
     &  ndastp,zyr,pco2at,pn2oat,pch4at
#endif
      RETURN
      END
