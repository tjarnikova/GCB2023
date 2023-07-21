      SUBROUTINE bgcsnk

      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
#if defined key_iomput
      USE iom
#endif
      IMPLICIT NONE
#if defined key_planktom  && key_top
!CC---------------------------------------------------------------------
!CC
!CC                       ROUTINE bgcsnk
!CC                     *****************
!CC
!CC  PURPOSE :
!CC  ---------
!CC         Compute vertical flux of particulate matter due to
!CC         gravitational sinking
!CC
!C   METHOD :
!C   -------
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
 
! local declarations
      INTEGER ji, jj, jk,jl
      REAL bacfer,maxlig,omeara,omecal
      REAL remara,remco3,remik,sedflx,siremin,snkspd(jpdsi:jpgoc)
      REAL xagg1,xagg2,xagg3,xagg4,xdens,xfeequi,xkeq,xlam1b,xlamfc,xpack
      REAL zdenom,zrfe2
      REAL xlibad,xlibas,xphdms,zbldmd,xcldmd
#  include "domzgr_substitute.h90"
! Pelagic model
!
!     Computation of the vertical sinking speed
!     -----------------------------------------
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!    
!     Sinking flux of pFe, POC, GOC, BSi, and Cal
!     ------------------------------------------------------------------
!     WHEN YOU CHANGE THIS, CHANGE THE SEDIMENT MODEL BELOW AS WELL
          snkpoc(ji,jj,jk+1) =                       &
     &           rn_snkpoc/rjjss*trn(ji,jj,jk,jppoc) &
     &           *rfact*tmask(ji,jj,jk+1)
          snksfe(ji,jj,jk+1) =                       &
     &           rn_snkpoc/rjjss*trn(ji,jj,jk,jpsfe) &
     &           *rfact*tmask(ji,jj,jk+1)
!
! 24,200,60 are the molar weight of dry organic material,CaCO3,SiO2
! 1040380,1660000,1245000 are the density of organic material,CaCO3,SiO2
! rn_snkgoc and rn_singoc are derived by fitting sinking speeds 
! calculated with the sedimentation function of Buitenhuis et al. 2001
! to an exponential function 
! maximum density is 1091314 for maximum sinking speed of 10 m/timestep
!
          xpack = 240.
!         xpack = 24.
!
          xdens = (trn(ji,jj,jk,jpgoc)*xpack+trn(ji,jj,jk,jpcal)*100.    &
     &     +trn(ji,jj,jk,jpara)*100.                                    &
     &     +trn(ji,jj,jk,jpdsi)*60.1)*1e6/                              &
     &     MAX(trn(ji,jj,jk,jpgoc)*xpack/1.08                        &
     &     +trn(ji,jj,jk,jpcal)*100./1.335                          &
     &     +trn(ji,jj,jk,jpara)*100./1.335                              &
     &     +trn(ji,jj,jk,jpdsi)*60.1/1.2,rtrn)-rhop(ji,jj,jk)*1000.

          xvsink(ji,jj,jk) = MIN(rn_singoc*MAX(xdens,dnsmin)**rn_snkgoc,snkmax(jk))
          snkgoc(ji,jj,jk+1) =                    &
     &            xvsink(ji,jj,jk)/rjjss*trn(ji,jj,jk,jpgoc) &
     &            *rfact*tmask(ji,jj,jk+1)
          snkgon(ji,jj,jk+1) =                    &
     &            xvsink(ji,jj,jk)/rjjss*trn(ji,jj,jk,jpgon) &
     &            *rfact*tmask(ji,jj,jk+1)
!
          snkbfe(ji,jj,jk+1) =                      &
     &              xvsink(ji,jj,jk)/rjjss*trn(ji,jj,jk,jpbfe) &
     &              *rfact*tmask(ji,jj,jk+1)
 
!
          snkdsi(ji,jj,jk+1) = xvsink(ji,jj,jk)*trn(ji,jj,jk,jpdsi) &
     &           /rjjss*rfact*tmask(ji,jj,jk+1)
!
          snkcal(ji,jj,jk+1) = xvsink(ji,jj,jk)*trn(ji,jj,jk,jpcal) &
     &           /rjjss*rfact*tmask(ji,jj,jk+1)
          snkara(ji,jj,jk+1) = xvsink(ji,jj,jk)*trn(ji,jj,jk,jpara)  &
     &           /rjjss*rfact*tmask(ji,jj,jk+1)
          END DO
        END DO
      END DO
!
!  Exchange between organic matter compartments due to
!  coagulation/disaggregation
!  ---------------------------------------------------
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
!    Part I : Coagulation dependent on turbulence
!    ----------------------------------------------
!
       xagg1=rn_ag1poc/rjjss*rfact*min(avt(ji,jj,jk)/5.E-4,1.)* &
     &       trn(ji,jj,jk,jppoc)**2
       xagg2=rn_ag2poc/rjjss*rfact*min(avt(ji,jj,jk)/5.E-4,1.)* &
     &       trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jpgoc)
!
!    Aggregation of small into large particles
!    Part II : Differential settling
!    -------------------------------------------------------------------
!
       xagg3=rn_ag3poc/rjjss*rfact*                  &
     &       trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jpgoc)
       xagg4=rn_ag4poc/rjjss*rfact* &
     &       trn(ji,jj,jk,jppoc)**2
       xagg(ji,jj,jk)=xagg1+xagg2+xagg3+xagg4
       xaggfe(ji,jj,jk)=xagg(ji,jj,jk)*trn(ji,jj,jk,jpsfe)/        &
     &                  (trn(ji,jj,jk,jppoc)+rtrn)*tmask(ji,jj,jk)
!
!     Aggregation of DOC to small particles
!     --------------------------------------
!
        xaggdoc(ji,jj,jk)=(rn_ag5doc*trn(ji,jj,jk,jpdoc)                &
     & +rn_ag7doc*trn(ji,jj,jk,jppoc))/rjjss*rfact                      &
     & *min(avt(ji,jj,jk)/5.E-4,1.)*trn(ji,jj,jk,jpdoc)*tmask(ji,jj,jk)
        xaggdoc2(ji,jj,jk)=rn_ag6doc*trn(ji,jj,jk,jpgoc)*rfact          &
     & /rjjss*min(avt(ji,jj,jk)/5.E-4,1.)*trn(ji,jj,jk,jpdoc)           &
     & *tmask(ji,jj,jk)
          END DO
        END DO
      END DO
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
!     TOC remineralization. Depends on temperature and nutrients
!     ------------------------------------------------------------------
!
            zdenom=rn_kmobac+rn_gbadoc*trn(ji,jj,jk,jpdoc) &
     &        +rn_gbapoc*trn(ji,jj,jk,jppoc) &
     &        +rn_gbagoc*trn(ji,jj,jk,jpgoc)
            remik=rn_grabac/rjjss*rfact*tmask(ji,jj,jk) &
     &        *tgfunc(ji,jj,jk,jpbac) &
     &        *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy)) &
     &        *trn(ji,jj,jk,jpbac)
!
! remineralization of DOC
!
            remdoc(ji,jj,jk) =remik*rn_gbadoc*trn(ji,jj,jk,jpdoc)/zdenom
!
! remineralization of POC
! WHEN YOU CHANGE THIS, CHANGE THE SEDIMENT MODEL BELOW AS WELL
            remik=rn_grabac/rjjss*rfact*tmask(ji,jj,jk) &
     &        *tgfunc(ji,jj,jk,jpbac) &
     &        *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy)) &
     &        *max(trn(ji,jj,jk,jpbac),rn_rembac)
            rempoc(ji,jj,jk)  =remik*rn_gbapoc*trn(ji,jj,jk,jppoc)/zdenom
!
! remineralization of GOC
!
            remgoc(ji,jj,jk) =remik*rn_gbagoc*trn(ji,jj,jk,jpgoc)/zdenom
            remgon(ji,jj,jk) =remik*rn_gbagon*trn(ji,jj,jk,jpgon)/zdenom! &
!     &        *trn(ji,jj,jk,jpgon)/trn(ji,jj,jk,jpgoc)
!
! remineralization of Iron in POC and GOC
!
            remsfe(ji,jj,jk)  =remik*rn_gbapoc*trn(ji,jj,jk,jpsfe)/zdenom
            rembfe(ji,jj,jk) =remik*rn_gbagoc*trn(ji,jj,jk,jpbfe)/zdenom
!
            bactge(ji,jj,jk)=rn_ggebac-rn_ggtbac*tsn(ji,jj,jk,1)
!
! Fe that is taken up by bacteria minus what is available in
! DOFe, sPOFe and bPOFe
!
            bacfer=bactge(ji,jj,jk)*ferat3*(remdoc(ji,jj,jk)             &
     &        +rempoc(ji,jj,jk)+remgoc(ji,jj,jk))                          &
     &        -remsfe(ji,jj,jk)-rembfe(ji,jj,jk)
!
! If there is not enough Fe in DOFe, sPOFe and bPOFe, then take up dissolved Fe
!
            ubafer(ji,jj,jk)=max(bacfer*trn(ji,jj,jk,jpfer)/            &
     &        (trn(ji,jj,jk,jpfer)+rn_kmfbac),0.)
!
! If there is not enough dissolved Fe, then decrease ggebac
!
            bactge(ji,jj,jk)=min(bactge(ji,jj,jk),                      &
     &        (ubafer(ji,jj,jk)+remsfe(ji,jj,jk)+rembfe(ji,jj,jk))         &
     &        /max((remdoc(ji,jj,jk)+rempoc(ji,jj,jk)+remgoc(ji,jj,jk))     &
     &        *ferat3,minfer))
!
! Fe in excess of that taken up by bacteria
!
            rbafer(ji,jj,jk)=max(-bacfer,0.)
            resbac(ji,jj,jk) = rn_resbac*(rn_retbac**tsn(ji,jj,jk,1)) &
     &        /rjjss*rfact*max(trn(ji,jj,jk,jpbac)-1e-10,0.)*tmask(ji,jj,jk)
!
!     Remineralisation rate of BSi dependent on T and O2
!     --------------------------------------------------
!
!          siremin=min(1.32E16*exp(-11200./(273.15+tsn(ji,jj,jk,1)))       &
!     &      ,0.1)/rjjss*rfact*tmask(ji,jj,jk) &
            siremin=min(rn_remdsi*exp(rn_retdsi/(273.15+tsn(ji,jj,jk,1))) &
     &        ,rn_readsi)/rjjss*rfact*tmask(ji,jj,jk) &
     &        *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy))
            remdsi(ji,jj,jk)  =siremin*trn(ji,jj,jk,jpdsi)

!
!     scavenging rate of iron based on Parekh et al., GBC 2005 as implemented by 
!     Aumont & Bopp GBC 2006
!     -------------------------------------------------------------------------
!
         xkeq    = 10**(17.27 - 1565.7 / ( 273.15 + tsn(ji,jj,jk,1) +         &
     &              (1.-tmask(ji,jj,jk))*20. ) )

         xfeequi = (-(1.+ligfer(jj,jk)*xkeq-xkeq*trn(ji,jj,jk,jpfer))+            &
     &             ((1.+ligfer(jj,jk)*xkeq-xkeq*trn(ji,jj,jk,jpfer))**2           &
     &             +4.*trn(ji,jj,jk,jpfer)*xkeq)**0.5)/(2.*xkeq)
         xlam1b  = rn_scmfer+rn_scofer*(trn(ji,jj,jk,jppoc)+trn(ji,jj,jk,jpgoc) &
     &             +trn(ji,jj,jk,jpcal)+trn(ji,jj,jk,jpdsi))*1E6

         xscave(ji,jj,jk) = xfeequi*xlam1b/rjjss*rfact*tmask(ji,jj,jk)
         scabfe(ji,jj,jk) = 1./(trn(ji,jj,jk,jppoc)+trn(ji,jj,jk,jpgoc)    &
     &                + trn(ji,jj,jk,jpdsi) + trn(ji,jj,jk,jpcal) + rtrn )
         scasfe(ji,jj,jk) = xscave(ji,jj,jk)*trn(ji,jj,jk,jppoc)/scabfe(ji,jj,jk)
         scabfe(ji,jj,jk) = xscave(ji,jj,jk)*trn(ji,jj,jk,jpgoc)/scabfe(ji,jj,jk)
#    if defined key_trc_dms
!
!     Bacterial production and degradation of DMS(Pd)
!
	 zbldmd=min(1., &
     &         max(0.66,1.-(etot(ji,jj,jk)/rn_etomax)**6.+0.18))
	 xlibad=min(xlimbac(ji,jj,jk),trn(ji,jj,jk,jpdmd)/      &
     &         (trn(ji,jj,jk,jpdmd)+rn_xkdmd))                 &
     &         * zbldmd
!
! cleavage of dmspd
!
         xcldmd=rn_xcldmd*rfact/rjjss*trn(ji,jj,jk,jpdmd)

         degdmd(ji,jj,jk)= 2*rn_grabac*xlibad  &
     &    *1.12**(tsn(ji,jj,jk,1))*2.*trn(ji,jj,jk,jpbac) &
     &    * rfact/rjjss *trn(ji,jj,jk,jpdmd)

         dmddms(ji,jj,jk)=rn_dmsyld*degdmd(ji,jj,jk)+xcldmd
         degdmd(ji,jj,jk)=degdmd(ji,jj,jk)+xcldmd
         xlibas=min(xlimbac(ji,jj,jk), trn(ji,jj,jk,jpdms)/ &
     &       (trn(ji,jj,jk,jpdms)+rn_xkdms)) &
     &         * zbldmd
!
! photolysis of dms
!
         degdms(ji,jj,jk)=2*(rn_grabac*xlibas  &
     &        *1.12**(tsn(ji,jj,jk,1))*0.6*trn(ji,jj,jk,jpbac) &
     &        +rn_xpodms*etot(ji,jj,jk))*rfact/rjjss*trn(ji,jj,jk,jpdms)
#    endif
          END DO
        END DO
      END DO
! Sediment model
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!     Sinking fluxes out of water layer above sediment
!     WHEN YOU CHANGE THIS, CHANGE THE PELAGIC MODEL ABOVE AS WELL
          jk=max(2,mbathy(ji,jj))-1
!     No T, BAC and OXY in the sediment, use values in overlying water
          remik=rn_grabac/rjjss*tmask(ji,jj,jk)*tgfunc(ji,jj,jk,jpbac) &
     &      *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy)) &
     &      *max(trn(ji,jj,jk,jpbac),rn_rembac)
          siremin=min(rn_remdsi*exp(rn_retdsi/(273.15+tsn(ji,jj,jk,1))) &
     &      ,rn_readsi)/rjjss*tmask(ji,jj,jk) &
     &      *(trn(ji,jj,jk,jpoxy)+3.E-6)/(10.E-6+trn(ji,jj,jk,jpoxy))
!     POM remineralisation in the sediment. No DOC in the sediment
          zdenom=rn_kmobac                                            &
     &      +rn_gbapoc*trnsed(ji,jj,jppoc) &
     &      +rn_gbagoc*trnsed(ji,jj,jpgoc)
          remsed(ji,jj,jppoc) =remik*rn_gbapoc*trnsed(ji,jj,jppoc)/zdenom
          remsed(ji,jj,jpgoc) =remik*rn_gbagoc*trnsed(ji,jj,jpgoc)/zdenom
          remsed(ji,jj,jpgon) =remik*rn_gbagon*trnsed(ji,jj,jpgon)/zdenom
          remsed(ji,jj,jpsfe) = remik*rn_gbapoc*trnsed(ji,jj,jpsfe)/zdenom
          remsed(ji,jj,jpbfe) = remik*rn_gbagoc*trnsed(ji,jj,jpbfe)/zdenom
          remsed(ji,jj,jpdsi) = siremin*trnsed(ji,jj,jpdsi)
          omecal = co3(ji,jj,jk)*calcon/aksp(ji,jj,jk)
          remsed(ji,jj,jpcal) = max(trnsed(ji,jj,jpcal)*rn_lyscal/rjjss*(1.-omecal),0.)
          omeara = co3(ji,jj,jk)*calcon/aksara(ji,jj,jk)
          remsed(ji,jj,jpara) = max(trnsed(ji,jj,jpara)*rn_lysara/rjjss*(1.-omeara),0.)
          xpack = 240.
          xdens = (trn(ji,jj,jk,jpgoc)*xpack+trn(ji,jj,jk,jpcal)*100.    &
     &     +trn(ji,jj,jk,jpara)*100.                                    &
     &     +trn(ji,jj,jk,jpdsi)*60.1)*1e6/                              &
     &     MAX(trn(ji,jj,jk,jpgoc)*xpack/1.08                        &
     &     +trn(ji,jj,jk,jpcal)*100./1.335                          &
     &     +trn(ji,jj,jk,jpara)*100./1.335                              &
     &     +trn(ji,jj,jk,jpdsi)*60.1/1.2,rtrn)-rhop(ji,jj,jk)*1000.
          xvsink(ji,jj,jk) = MIN(rn_singoc*MAX(xdens,dnsmin)**rn_snkgoc,snkmax(jk))
          snkspd=xvsink(ji,jj,jk)/rjjss
          snkspd(jppoc)=rn_snkpoc/rjjss
          snkspd(jpsfe)=rn_snkpoc/rjjss
          DO jl = jpdsi, jpgoc
            sedflx = snkspd(jl)*trn(ji,jj,jk,jl)*tmask(ji,jj,jk)/fse3t(ji,jj,jk)
            remsed(ji,jj,jl) = max(min((trnsed(ji,jj,jl)-rn_gramin)*rfactr  &
     &        +sedflx,remsed(ji,jj,jl)),0.)
            trn(ji,jj,jk,jl) = trn(ji,jj,jk,jl)-sedflx*rfact
            trnsed(ji,jj,jl) = trnsed(ji,jj,jl)+(sedflx-remsed(ji,jj,jl))*rfact
          END DO
!!     Nutrients are instantaneously added to the overlying water
          trn(ji,jj,jk,jpsil) = trn(ji,jj,jk,jpsil)+remsed(ji,jj,jpdsi)*rfact
          trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)                     &
     &      +(remsed(ji,jj,jppoc)+remsed(ji,jj,jpgon)*ratc2n)*rfact
          remik = remsed(ji,jj,jppoc)+remsed(ji,jj,jpgoc)
          sedflx = min(remik*rato2c,trn(ji,jj,jk,jpoxy)*rfactr)
          trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)-sedflx*rfact
!! For now, in contrast to the water column, denitrification starts when O2=0.
          trn(ji,jj,jk,jpdin) = trn(ji,jj,jk,jpdin)+(0.8*(sedflx-remik*rato2c) &
     &      +remsed(ji,jj,jppoc)*ratn2c+remsed(ji,jj,jpgon))*rfact
          trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                     &
     &      +(remsed(ji,jj,jpsfe)+remsed(ji,jj,jpbfe))*rfact
          trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)+(remik              &
     &      +remsed(ji,jj,jpcal)+remsed(ji,jj,jpara))*rfact
#  if defined key_trc_piic
          trn(ji,jj,jk,jppiic) = trn(ji,jj,jk,jppiic)+(remik            &
     &      +remsed(ji,jj,jpcal)+remsed(ji,jj,jpara))*rfact
#  endif
          trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)+(-alknut*remik      &
     &      +2.*remsed(ji,jj,jpcal)+2.*remsed(ji,jj,jpara))*rfact
        END DO
      END DO
#if defined key_iomput
!      where (tmask(:,:,:) .eq. 0 ) 
!        out3d = ncf_fill
!      else where
        out3d = xvsink/rjjss
!      end where
      CALL iom_put("vsink", out3d )
      CALL iom_put("remdsised",remsed(:,:,jpdsi)*1e3)
#endif
 
#endif
      RETURN
      END SUBROUTINE bgcsnk
