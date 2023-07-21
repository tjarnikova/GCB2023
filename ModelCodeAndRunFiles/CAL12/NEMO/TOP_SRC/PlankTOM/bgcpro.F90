#if defined key_planktom && defined key_top
      SUBROUTINE bgcpro
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE bgcpro
!!!                     *****************
!!!
!!!  PURPOSE :
!!!  ---------
!!!         Compute the phytoplankton production depending on
!!!         light, temperature and nutrient availability
!!!
!!   METHOD :
!!   -------
!!      
!!
!!   INPUT :
!!   -----
!!      argument
!!              None
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
!!      original  p4zprod : O. Aumont (2002) 
!!      modifications     : E. Buitenhuis(2002)
!!----------------------------------------------------------------------
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
!!      USE sbc_oce , ONLY :   qsr_oce        =>    qsr_oce        !: penetrative solar radiation (w m-2) 
      USE traqsr  , ONLY : ln_qsr_sms
      USE iom
      IMPLICIT NONE
! local declarations
      INTEGER ji, jj, jk, jl
      REAL silpot1,silpot2,silfac,pislopen,ysopt
      REAL parlux,xchl,ekg,ekr,xlim1,xlim2,xlim3,xlim4(jpdia:jpdia+jppft-1)
      REAL xlim6(jpdia:jpdia+jppft-1),dinlim
      REAL xlim7,xlim8
! local dgom variables
      REAL pcphot(jpdia:jpdia+jppft-1),quopfe(jpdia:jpdia+jppft-1)
      REAL rhochl(jpdia:jpdia+jppft-1),vcfer(jpdia:jpdia+jppft-1)
      REAL perfrm(jpdia:jpdia+jppft-1)
      REAL pctnut,docpro
      REAL zmask,zkrdphy
!
#  include "domzgr_substitute.h90"
! etot always calculated in traqsr.F90
      DO jl = jpdia, jpdia+jppft-1
        pcmax(jl) = rn_mumpft(jl)*(1.+rn_resphy(jl))/rjjss
      END DO
!
      DO jl = jpdia+1, jpdia+jppft-1
        xlim5(jl)=1.
      END DO
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
!    Input of river carbon is balanced by carbon removal from the
!    bottom water layer in bgcsed. Because the iron cycle is not
!    closed, this does not work for iron, so particulate organic
!    carbon tends to accumulate in the bottom layer. To set a bound
!    on that, the sfe/poc and bfe/goc is kept below a maximum ratio.
!    ---------------------------------------------------------------
!
!            trn(ji,jj,jk,jpsfe)=trn(ji,jj,jk,jppoc)*min(1.E-3, &
!     &        trn(ji,jj,jk,jpsfe)/(rtrn+trn(ji,jj,jk,jppoc)))
!            trn(ji,jj,jk,jpbfe)=trn(ji,jj,jk,jpgoc)*min(1.E-3, &
!     &        trn(ji,jj,jk,jpbfe)/(rtrn+trn(ji,jj,jk,jpgoc)))

!            trn(ji,jj,jk,jpsfe)=trn(ji,jj,jk,jppoc)*max(rtrn, &
!     &        trn(ji,jj,jk,jpsfe)/(rtrn+trn(ji,jj,jk,jppoc)))
!            trn(ji,jj,jk,jpbfe)=trn(ji,jj,jk,jpgoc)*max(rtrn, &
!    &        trn(ji,jj,jk,jpbfe)/(rtrn+trn(ji,jj,jk,jpgoc)))


!
!
!      Michaelis-Menten Limitation term for nutrients with threshold for growth at 0.1*k-half
!      -----------------
! Diatoms
             jl = jpdia
!
             xlim4(jl)=(trn(ji,jj,jk,jppo4)-rn_kmpphy(jl)*rn_nutthe)&
     &                /(trn(ji,jj,jk,jppo4)+rn_kmpphy(jl)*(1.-rn_nutthe))
!
             xlim5(jl)=(trn(ji,jj,jk,jpsil)-rn_sildia*rn_nutthe)    &
     &                /(trn(ji,jj,jk,jpsil)+rn_sildia*(1.-rn_nutthe))
!
             xlim6(jl)=(trn(ji,jj,jk,jpdin)-rn_kmnphy(jl)*rn_nutthe)&
     &                /(trn(ji,jj,jk,jpdin)+rn_kmnphy(jl)*(1.-rn_nutthe))
!
           DO jl = jpdia+1, jpdia+jppft-2
!
             xlim4(jl)=(trn(ji,jj,jk,jppo4)-rn_kmpphy(jl)*rn_nutthe)&
     &                /(trn(ji,jj,jk,jppo4)+rn_kmpphy(jl)*(1.-rn_nutthe))
!
             xlim6(jl)=(trn(ji,jj,jk,jpdin)-rn_kmnphy(jl)*rn_nutthe)&
     &                /(trn(ji,jj,jk,jpdin)+rn_kmnphy(jl)*(1.-rn_nutthe))
!
           END DO
!    Michaelis-Menten Limitation term for nutrients for bacteria -
! for use in degradation of DMS
            xlim1=trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+rn_kmpbac)
            xlim2=trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)+rn_kmfbac)
            xlim3=trn(ji,jj,jk,jpdoc)/(trn(ji,jj,jk,jpdoc)+rn_kmobac)
            xlimbac(ji,jj,jk)=min(xlim1,xlim2,xlim3)
! N2 fixers
           jl = jpfix
! michaelis menten P
!
             xlim4(jl)=(trn(ji,jj,jk,jppo4)-rn_kmpphy(jl)*rn_nutthe)&
     &                /(trn(ji,jj,jk,jppo4)+rn_kmpphy(jl)*(1.-rn_nutthe))
             dinlim =  (trn(ji,jj,jk,jpdin)-rn_kmnphy(jl)*rn_nutthe)&
     &                /(trn(ji,jj,jk,jpdin)+rn_kmnphy(jl)*(1.-rn_nutthe))
!
! michaelis menten N
!
             xlim6(jl)= dinlim +rn_munfix*(1.-dinlim)
             dinpft(ji,jj,jk,jl)=dinlim/(xlim6(jl)+rtrn)*tmask(ji,jj,jk)
!
             DO jl = jpdia, jpdia+jppft-1
!
! quota model for Fe-light
!
              stofoo(ji,jj,jk,jl,2) = trn(ji,jj,jk,jl+jppft)     &
     &          /(trn(ji,jj,jk,jl)+rtrn)
              stofoo(ji,jj,jk,jl,3) = trn(ji,jj,jk,jl+2*jppft)   &
     &          /(trn(ji,jj,jk,jl)+rtrn)
            quopfe(jl) =max(min(stofoo(ji,jj,jk,jl,2),rn_qmaphy(jl)),rn_qmiphy(jl))
            xlim1 = (rn_rhfphy(jl)*rn_qmaphy(jl)-rn_qmaphy(jl))*    &
     &        (rn_qmaphy(jl)-quopfe(jl))/                           &
     &        (rn_qmaphy(jl)-rn_qmiphy(jl))+rn_qmaphy(jl)
            xlim2 = trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)+       &
     &        rn_kmfphy(jl))
            xlim3 =min((quopfe(jl)-rn_qmiphy(jl))                       &
     &        /(rn_qopphy(jl)-rn_qmiphy(jl)),1.)*(1.+rn_nutthe)-rn_nutthe
!            xlim7 = (4.6*etot(ji,jj,jk)-rn_tliphy(jl))*rn_nutthe/       &
!     &        (4.6*etot(ji,jj,jk)+rn_tliphy(jl)*(1.-rn_nutthe))
            xlimpft(ji,jj,jk,jl)=min(xlim4(jl),xlim5(jl),xlim6(jl),xlim3)
! 
! Fe uptake rate 
!
            vcfer(jl) = rn_mumpft(jl)*(1.+rn_resphy(jl))                   &
     &        *xlim1*min(xlim4(jl),xlim5(jl),xlim6(jl),xlim2)
!
! 4.6 micromol photons/J at 550 nm
!
            perfrm(jl)=max(rn_alpphy(jl)*stofoo(ji,jj,jk,jl,3)          &
     &        *(4.6*etot(ji,jj,jk)-rn_tliphy(jl)),0.)
            docpro=rn_docphy(jl)+(1.-xlimpft(ji,jj,jk,jl))*rn_domphy(jl)
           pctnut=pcmax(jl)*xlimpft(ji,jj,jk,jl)*tgfunc(ji,jj,jk,jl)

! light limitation
            xlim8      = (1.-exp(-perfrm(jl)/(pctnut+rtrn)))
            pcphot(jl) = pctnut*xlim8
            rhochl(jl)=rn_thmphy(jl)*pcphot(jl)/(perfrm(jl)+rtrn)
!
! synthesis rates
!
            prophy(ji,jj,jk,jl,3) = rhochl(jl)*pcphot(jl)               &
     &        *trn(ji,jj,jk,jl)*rfact
            prophy(ji,jj,jk,jl,2) = vcfer(jl)*xlim8*tgfunc(ji,jj,jk,jl) &
     &        *trn(ji,jj,jk,jl)*rfact/rjjss
            prophy(ji,jj,jk,jl,1) = pcphot(jl)*trn(ji,jj,jk,jl)*rfact
            docphy(ji,jj,jk,jl)   = prophy(ji,jj,jk,jl,1)*docpro
           END DO
!
!    FE/C and Si/C of diatoms
!    ------------------------
!    Si/C increases with iron stress and silicate availability
!    Si/C is increased for very high Si concentrations
!    to mimic the very high ratios observed in the Southern Ocean
!    (silpot2)
!    ET20220901: Although at high SiO3 Si:C increases (silpot2)
!     it stops decreasing at limiting SiO3 (doi:10.5194/bg-7-657-2010)
!     The growth-limiting SiO3 term was therefore moved from ysopt to silpot1
!     doi:10.4319/lo.2004.49.4.1134 supports, but doi:10.3354/meps195071 contradicts silpot1
            silpot1=1.+rn_ferbsi*min(1.,trn(ji,jj,jk,jpsil)/rn_sildia)* &
     &        (1.-min(trn(ji,jj,jk,jpfer)/rn_kmfphy(jpdia),1.))
            silpot2=rn_silbsi*trn(ji,jj,jk,jpsil)/(trn(ji,jj,jk,jpsil)+rn_kmsbsi)
            silfac=max(silpot1,silpot2)
            ysopt=rn_bsidia*silfac
            prorca3(ji,jj,jk) = prophy(ji,jj,jk,jpdia,1)*ysopt* &
     &        tmask(ji,jj,jk)
          END DO
        END DO
      END DO
      RETURN
      END
#endif
