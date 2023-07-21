       SUBROUTINE bgcprg(kt)
#if defined key_planktom &&  defined key_top
!!!---------------------------------------------------------------------

!!!
!!!                       ROUTINE bgcprg
!!!                     *****************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     Call Biological sources and sinks subroutines
!!!
!!   METHOD :
!!   -------
!!      
!!
!!   INPUT :
!!   -----
!!      argument
!!              ktask           : task identificator
!!              kt              : time step
!!      common
!!              all the common defined in opa
!!
!!
!!   OUTPUT :                   : no
!!   ------
!!
!!   WORKSPACE :
!!   ---------
!!
!!   EXTERNAL :
!!   --------
!!
!!   MODIFICATIONS:
!!   --------------
!!      original  : E. Maier-Reimer (GBC 1993): h3cprg
!!      additions : O. Aumont (1998)
!!      additions : C. Le Quere (1999)
!!      additions : O. Aumont (2001)
!!      additions : O. Aumont , EK (05/2001): add h3cadj 
!!----------------------------------------------------------------------
!! parameters and commons
!! ======================
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE lbclnk
      USE lib_mpp
      IMPLICIT NONE
!!----------------------------------------------------------------------
!! local declarations
!! ==================
      INTEGER ktask, kt, jn
      INTEGER ji, jj, jk
      REAL(wp) total(jptra), totaf(jptra),totic,totif,totec,totef, units
!
      IF(lwp) WRITE(numout,*) 'hjdjdh bgcprg_s ',trn(nlci-1,nlcj-1,max(2,mbathy(nlci-1,nlcj-1)),jpgoc)
      IF (kt .EQ. nit000) THEN
        qcumul = 0.
      END IF
        total = 0.
        DO jn = 1, jptra
          DO jk = 1, jpk
            DO jj = 2, nlcj-1
              DO ji = 2, nlci-1
                total(jn) = total(jn)+trn(ji,jj,jk,jn)*volumt(ji,jj,jk)
              END DO
            END DO
          END DO
        END DO
#  if defined key_mpp_mpi
        CALL mpp_sum(total,jptra)
#  endif
        totic=total(jpdoc)
        DO jn = jppoc, jpdia+jppft-1
          totic=totic+total(jn)
        END DO
        totif=total(jpsfe)+total(jpbfe)
        DO jn = jpbac, jpbac+jpzft
          totif=totif+total(jn)*ferat3
        END DO
        DO jn = jpdfe, jpdfe+jppft-1
          totif=totif+total(jn)
        END DO
      IF (kt .EQ. nit000) THEN
        IF (lwp) WRITE(numout,*) 'bgcprg start C,P,Fe,Si,O2,jptra' &
     &    ,totic+total(jpdic),totic+total(jppo4) &
     &    ,totif+total(jpfer),total(jpsil)+total(jpbsi)+total(jpdsi) &
     &    ,total(jpoxy)-rato2c*totic,total
      END IF
!
! Compute chemical variables
! --------------------------
!
          CALL bgcche
!......................................................................
!
! Interpolate chemical variables
! ------------------------------
!
          CALL bgcint(kt)
!
!......................................................................
!
! Compute CaCO3 saturation
! ------------------------
!
          CALL bgclys
!
!......................................................................
!
! Compute biology (POC)
! ------------------------------------
!
          CALL bgcbio(kt)
!
!......................................................................
!
! Close budgets
! ------------------------------------
!
          CALL bgcsed
          IF ( kt == nitend ) THEN
#  if defined key_mpp_mpi
            CALL mpp_sum(sedcor,6)
#  endif
            units = rfact/raass*100.
            IF(lwp) WRITE(numout,*)' bottom water correction lacked ',  &
     &        sedcor,' mol Si, C, Fe, DOC, Alk, which is ',             &
     &      (sedcor(ji)/extinp(ji)*units,ji=1,6),'% of external inputs.'
#    if defined key_trc_dms
            IF(lwp) WRITE(numout,*)' dms sinks converted to sources ', &
     &          dms_snkcount, ' times'
             IF(lwp) WRITE(numout,*)' dmd sinks converted to sources ', &
     &          dmd_snkcount, ' times'
#    endif
          ENDIF
!
!......................................................................
!
! Compute surface fluxes
! ----------------------
!
          CALL bgcflx
!      IF (kt .EQ. nit000) THEN
        totaf = 0.
        DO jn = 1, jptra
          DO jk = 1, jpk-1
            DO jj = 2, nlcj-1
              DO ji = 2, nlci-1
                totaf(jn) = totaf(jn)+trn(ji,jj,jk,jn)*volumt(ji,jj,jk)
              END DO
            END DO
          END DO
        END DO
#  if defined key_mpp_mpi
        CALL mpp_sum(totaf,jptra)
#  endif
        DO jn = 1, jptra
          qcumul(jn)=qcumul(jn)+totaf(jn)-total(jn)
        END DO
!      END IF
      IF (kt .EQ. nitend) THEN
        totec=totaf(jpdoc)
        DO jn = jppoc, jpdia+jppft-1
          totec=totec+totaf(jn)
        END DO
        totef=total(jpsfe)+total(jpbfe)
        DO jn = jpbac, jpbac+jpzft
          totef=totef+totaf(jn)*ferat3
        END DO
        DO jn = jpdfe, jpdfe+jppft-1
          totef=totef+totaf(jn)
        END DO
       IF (kt .EQ. nitend) THEN
        IF (lwp) WRITE(numout,*) 'bgcprg end C,P,Fe,Si,O2,jptra' &
     &    ,totec+totaf(jpdic),totec+totaf(jppo4) &
     &    ,totef+totaf(jpfer),totaf(jpsil)+totaf(jpbsi)+totaf(jpdsi) &
     &    ,totaf(jpoxy)-rato2c*totec,totaf
        totec=qcumul(jpdoc)
        DO jn = jppoc, jpdia+jppft-1
          totec=totec+qcumul(jn)
        END DO
        totef=total(jpsfe)+total(jpbfe)
        DO jn = jpbac, jpbac+jpzft
          totef=totef+qcumul(jn)*ferat3
        END DO
        DO jn = jpdfe, jpdfe+jppft-1
          totef=totef+qcumul(jn)
        END DO
        IF (lwp) WRITE(numout,*) 'bgcprg budget C,P,Fe,Si,O2,jptra' &
     &    ,totec+qcumul(jpdic),totec+qcumul(jppo4) &
     &    ,totef+qcumul(jpfer),qcumul(jpsil)+qcumul(jpbsi)+qcumul(jpdsi) &
     &    ,qcumul(jpoxy)-rato2c*totec,qcumul
       END IF
!        IF (lwp) WRITE(numout,*) 'bgcprg end C,P,Fe,O2,jptra' &
!     &    ,totof+totaf(jpdic)-totoc-total(jpdic)-qcumul(1) &
!     &    ,totof+totaf(jppo4)-totoc-total(jppo4) &
!     &    ,totaf(jpsfe)+totaf(jpbfe)+(totaf(jpmic) &
!     &    +totaf(jpmes))*ferat3+totaf(jpfer) &
!     &    +totaf(jpdfe)+totaf(jpnfe)+totaf(jpcfe) &
!     &    -(total(jpsfe)+total(jpbfe)+(total(jpmic) &
!     &    +total(jpmes))*ferat3+total(jpfer) &
!     &    +total(jpdfe)+total(jpnfe)+total(jpcfe)) &
!     &    ,totaf(jpoxy)-rato2c*(totof)-total(jpoxy)+rato2c*totoc-qcumul(2),totaf
      END IF
      DO jn=1 , jptra
        CALL lbc_lnk(trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(tra(:,:,:,jn), 'T', 1. )
      END DO
      trb(:,:,:,:)=trn(:,:,:,:)
#endif
!
      RETURN
      END

