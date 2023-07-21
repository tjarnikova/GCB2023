      SUBROUTINE bgcsed
#if defined key_planktom && defined key_top
!CC---------------------------------------------------------------------
!CC
!CC                       ROUTINE bgcsed
!CC                     *****************
!CC
!CC  PURPOSE :
!CC  ---------
!CC         Compute loss of organic matter in the sediments. This
!CC         is by no way a sediment model. The loss is simply 
!CC         computed to balance the input from rivers and dust
!CC
!C   METHOD :
!C   -------
!C      
!C
!C   INPUT :
!C   -----
!C
!C   OUTPUT :                   : no
!C   ------
!C
!C   WORKSPACE :
!C   ---------
!C
!C   EXTERNAL :
!C   --------
!C
!C   MODIFICATIONS:
!C   --------------
!C      original  : O. Aumont 
!C----------------------------------------------------------------------
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE lib_mpp
      USE dom_oce , ONLY :   mbathy     =>   mbathy     !: number of ocean level (=0,  & 1, ... , jpk-1) 
      IMPLICIT NONE
!C
!C local declarations
!C ==================
      INTEGER ji, jj, ikt
      REAL dummyv(5),sumsed(12),botcor(10),litres
#  include "domzgr_substitute.h90"
!
!    Loss of biogenic silicon and organic carbon in the sediments. 
!    First, the total inventory is computed.
!    -------------------------------------------------------------------
!
      sumsed=0.
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          ikt=max(mbathy(ji,jj)-1,1)
          litres=1000.*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,ikt)
          sumsed(1)=sumsed(1)+trnsed(ji,jj,jpdsi)*litres
          sumsed(2)=sumsed(2)+trnsed(ji,jj,jpgoc)*litres
          sumsed(7)=sumsed(7)+trnsed(ji,jj,jppoc)*litres
          sumsed(3)=sumsed(3)+trnsed(ji,jj,jpbfe)*litres
          sumsed(8)=sumsed(8)+trnsed(ji,jj,jpsfe)*litres
          sumsed(4)=sumsed(4)+trn(ji,jj,ikt,jpdoc)*litres
          sumsed(5)=sumsed(5)+trnsed(ji,jj,jpcal)*litres
          sumsed(10)=sumsed(10)+trnsed(ji,jj,jpara)*litres
          sumsed(11)=sumsed(11)+trnsed(ji,jj,jpgon)*litres
          sumsed(12)=sumsed(12)+trn(ji,jj,ikt,jpdic)*litres
        END DO
      END DO
      IF( lk_mpp ) CALL mpp_sum( sumsed,12)
!
!    Then the surface inputs of each element are scaled against the 
!    sediment inventories, and the bottom concentrations are multiplied
!    by this scaling factor. Thus, the amount lost in the sediments equals
!    the supply at the surface (dust+rivers)
!    -------------------------------------------------------------------
!
      DO ji = 1, 5
        dummyv(ji)=sumsed(ji)+sumsed(ji+5)
        IF (extinp(ji) .LE. dummyv(ji)) THEN
          botcor(ji)  =1.-extinp(ji)/dummyv(ji)
        ELSE
          botcor(ji)  =0.
          sedcor(ji)  =sedcor(ji)+extinp(ji)-dummyv(ji)
        ENDIF
      END DO
      IF (extinp(6)-sumsed(7)*ratn2c .LE. sumsed(11)) THEN
        botcor(6)  =1.-max(extinp(6)-sumsed(7)*ratn2c,0.)/sumsed(11)
        sedcor(6)  =sedcor(6)+min(extinp(6)-sumsed(7)*ratn2c,0.)
      ELSE
        botcor(6)  =0.
        sedcor(6)  =sedcor(6)+extinp(6)-sumsed(7)*ratn2c-sumsed(11)
      ENDIF
      IF (extinp(7)-sumsed(7)*ratn2c .LE. sumsed(12)) THEN
        botcor(7)  =1.-extinp(7)/sumsed(12)
      ELSE
        botcor(7)  =0.
        sedcor(7)  =sedcor(7)+extinp(7)-sumsed(12)
      ENDIF
      IF (lwp) write(numout,31) "sediment inventory (Tmol):", (dummyv(ji) &
     &  *1.e-12,ji=1,5),sumsed(11)*1.e-12,sumsed(12)*1.e-12
31    FORMAT(a26,7f10.3)
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          ikt=max(mbathy(ji,jj)-1,1)
          trnsed(ji,jj,jpdsi)=trnsed(ji,jj,jpdsi)*botcor(1)
          trnsed(ji,jj,jpgoc)=trnsed(ji,jj,jpgoc)*botcor(2)
          trnsed(ji,jj,jppoc)=trnsed(ji,jj,jppoc)*botcor(2)
          trnsed(ji,jj,jpbfe)=trnsed(ji,jj,jpbfe)*botcor(3)
          trnsed(ji,jj,jpsfe)=trnsed(ji,jj,jpsfe)*botcor(3)
          trn(ji,jj,ikt,jpdoc)=trn(ji,jj,ikt,jpdoc)*botcor(4)
          trnsed(ji,jj,jpcal)=trnsed(ji,jj,jpcal)*botcor(5)
          trnsed(ji,jj,jpara)=trnsed(ji,jj,jpara)*botcor(5)
          trnsed(ji,jj,jpgon)=trnsed(ji,jj,jpgon)*botcor(6)
          trn(ji,jj,ikt,jpdic)=trn(ji,jj,ikt,jpdic)*botcor(7)
        END DO
      END DO
#endif
      RETURN
      END
