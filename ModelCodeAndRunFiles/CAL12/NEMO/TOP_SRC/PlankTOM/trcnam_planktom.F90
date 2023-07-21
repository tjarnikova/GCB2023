MODULE trcnam_planktom
#if defined key_planktom
!!!
!!!                       trcnam_planktom
!!!                       ****************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     READs options for DGOM namelist
!!!     Variables that affect phytoplankton have been put in arrays
!!!     This is incompatible with the previous use of #ifdef key
!!!     Therefore these keys have been removed and this version is only
!!!     suited for the green ocean model.
!!!
!!   METHOD :                   : no
!!   -------
!!
!!   INPUT :
!!   -----
!!            &natgas           : gas exchange parameters
!!            &natagg           : particulate aggregation parameters 
!!            &natriv           : river elemental input parameters 
!!            &natfer           : iron model parameters 
!!            &natzoo           : biological parameters related to heterotrophs
!!            &natphy           : biological parameters related to phyto
!!            &natpft           : biological parameters related to all PFTs
!!            &natquo           : iron quota model parameters 
!!            &natlit           : light-iron model parameters 
!!            &natdms           : DMS parameters 
!!
!!   OUTPUT :
!!   ------
!!      COMMON
!!           /cotgas/
!!           /cotcon/
!!           /cotbio/
!!
!!   WORKSPACE :                : no
!!   ---------
!!
!!   MODIFICATIONS:
!!   --------------
!!      original  : 99-10 (M.A. Foujols, M. Levy) passive tracer
!!      addition  : 00-01 (L. Bopp) hamocc3,p3zd
!!      modification : 02 (E. Buitenhuis) dgom
!!---------------- !!----------------------------------------------------------------------
   !! trc_nam_planktom       : PlankTOM model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_planktom    ! sms trends
   USE iom  ! I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_planktom   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_planktom.F90 2567 2011-01-25 09:36:27Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_planktom
!!------------------------------------------------------
!!----------------------------------------------------------------------
!! local declarations
!! ==================

      INTEGER jl, jm, jn
!----------------------------------------------------------------------
!! statement functions
!! ===================
!!#include "stafun.h"
!!
!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (15/11/96)
!!!---------------------------------------------------------------------
!
! 0. initializations
! ------------------
!
       namelist/natagg/rn_ag1poc,rn_ag2poc,rn_ag3poc,rn_ag4poc,rn_ag5doc,rn_ag6doc,rn_ag7doc,&
     &                 rn_singoc,rn_snkgoc,rn_snkpoc,rn_snkspd
       namelist/natriv/rn_rivdic,rn_rivdoc,rn_rivfer,rn_rivnit,rn_rivpo4,rn_rivpoc,rn_rivsil
       namelist/natfer/rn_fersol,rn_silsol,rn_sedfer,rn_scofer,rn_scmfer,rn_ligdep,rn_ligfer,rn_liglat
       namelist/natpre/rn_gbadoc,rn_gbagoc,rn_gbagon,rn_gbapoc,rn_prfzoo,rn_readsi,rn_remdsi,rn_retdsi
       namelist/natzoo/nn_sizzoo,rn_ggebac,rn_ggtbac,rn_ggezoo,&
     &   rn_sigzoo,                                  &
     &   rn_unazoo,                                  &
     &   rn_grabac,rn_gramin,rn_grazoo,                        &
     &   rn_kmobac,rn_grkzoo,                        &
     &   rn_kmfbac,rn_kmpbac,                                            &
     &   rn_mormac,rn_motmac,                        &
#if ! defined key_trc_foram
     &   rn_morgel,rn_motgel,                        &
#endif
     &   rn_rembac,rn_resbac,rn_reszoo,              &
     &   rn_retbac,rn_retzoo,                        &
     &   rn_icemac,rn_trnmac
       namelist/natphy/rn_bsidia,rn_coccal,rn_discal,rn_ferbsi,rn_kmsbsi,rn_silbsi,&
     &   rn_lyoco3,rn_lyscal,rn_sildia,rn_munfix,                        &
     &   rn_nutthe,                                            &
     &   rn_resphy,                                            &
     &   rn_docphy,rn_domphy,                                            &
     &   rn_kmnphy,rn_kmpphy,rn_mumpft
       namelist/natpft/rn_mutpft,rn_mudpft
       namelist/natquo/rn_kmfphy,rn_rhfphy,rn_qmaphy,rn_qmiphy,rn_qopphy
       namelist/natlit/rn_ekwgrn,rn_ekwred,rn_alpphy,rn_kgrphy,rn_krdphy,rn_thmphy,rn_tliphy
       namelist/natzca/rn_disara,rn_disfor,rn_forcal,rn_pteara,rn_lysara
#    if defined key_trc_dms
       namelist/natdms/rn_dmsyld,rn_etomax,rn_rdddms,rn_xcldmd,rn_xpodms,&
     &                 rn_xkdms,rn_xkdmd,rn_assdms,rn_rphdmd,rn_xpldmd
#    endif
#    if defined key_trc_n2o
       namelist/natn2o/nn_deun2s,rn_aoun2s,rn_atmn2o,rn_betn2s,rn_decn2s,rn_omxn2s
#    endif
#    if defined key_trc_ch4
       namelist/natch4/rn_proch1,rn_proch2,rn_proch3,rn_proch4,rn_proch5,&
     &                 rn_botch1,rn_botch2,rn_botch3,rn_botch4,rn_botch5,&
     &                 rn_conch1,rn_conch2,rn_conch3,rn_conch4,rn_conch5,&
     &                 rn_decch4
#    endif
#  if defined key_trc_diaadd
      namelist/natadd/ctrc3d,ctrc3l,ctrc2d,ctrc2l, ctrc3u, ctrc2u,     &
         nn_writedia
#  endif

        CHARACTER (len=39) ::   clname
        integer :: numnatpl
! initialize the number of LOGICAL UNIT used
!
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' IN trclsm.dgom.h90'
          WRITE(numout,*) ' ******************'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' Read namelist for DGOM model'
          WRITE(numout,*) ' ***********************'
          WRITE(numout,*) ' '
          CALL FLUSH(numout)
      ENDIF
#  if defined key_trc_diaadd && ! key_iomput
      STOP 'key_trc_diaadd without key_iomput is no longer supported'
#  endif
    
      numnatpl=80
      clname='namelist.trc.sms'
      OPEN( numnatpl, FILE= clname, FORM='formatted', STATUS = 'old')
!
! 2 Namelist :
!
      READ(numnatpl,natagg)
      READ(numnatpl,natriv)
      READ(numnatpl,natfer)
      READ(numnatpl,natpre)
      READ(numnatpl,natzoo)
      READ(numnatpl,natphy)
      READ(numnatpl,natpft)
      READ(numnatpl,natquo)
      READ(numnatpl,natlit)
      READ(numnatpl,natzca)
#    if defined key_trc_dms
      READ(numnatpl,natdms)
#    endif
#    if defined key_trc_n2o
      READ(numnatpl,natn2o)
#    endif
#    if defined key_trc_ch4
      READ(numnatpl,natch4)
#    endif
      REWIND( numnatpl ) 
!
      dnsmin=(rn_snkpoc/rn_singoc)**(1/rn_snkgoc)
!     write(numout,*) "namelist input",rn_singoc,rn_snkpoc,rn_snkgoc,dnsmin
      jpfoo = 0
      DO jm = 1, jpzft
        DO jl = jppoc,jpdia+jppft-1
          IF (rn_prfzoo(jm,jl) .GT. rtrn) THEN
            jpfoo=jpfoo+1
            IF (jpfoo .GT. jpfmx) STOP 'increase jpfmx'
            grizoo(1,jpfoo) = jm
            grizoo(2,jpfoo) = jl
          ENDIF
        END DO
      END DO
      DO jl = jpdia, jpdia+jppft-1
        IF (rn_resphy(jl) .eq. 1. .and. rn_nutthe .ne. 0.) STOP 'set rn_nutthe cannot gt 0 when rn_resphy eq 1 '
        IF (rn_resphy(jl) .eq. 1. .and. rn_tliphy(jl) .ne. 0.) STOP 'set rn_tliphy cannot gt 0 when rn_resphy eq 1 '
      ENDDO
!
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natagg'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' aggregation term 1 (rn_ag1poc) =', rn_ag1poc
          WRITE(numout,*) ' aggregation term 2 (rn_ag2poc) =', rn_ag2poc
          WRITE(numout,*) ' aggregation term 3 (rn_ag3poc) =', rn_ag3poc
          WRITE(numout,*) ' aggregation term 4 (rn_ag4poc) =', rn_ag4poc
          WRITE(numout,*) ' aggregation term 5 (rn_ag5doc) =', rn_ag5doc
          WRITE(numout,*) ' aggregation term 6 (rn_ag6doc) =', rn_ag6doc
          WRITE(numout,*) ' aggregation term 7 (rn_ag7doc) =', rn_ag7doc
          WRITE(numout,*) ' POC sinking speed (rn_snkpoc)  =', rn_snkpoc
          WRITE(numout,*) ' Big particles sinking multiplication (rn_snkgoc) =', rn_snkgoc
          WRITE(numout,*) ' Big particles sinking division (rn_singoc) =', rn_singoc
          WRITE(numout,*) ' Maximum sink speed (rn_snkspd) =', rn_snkspd

!
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natriv'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' river conversion of dic (rn_rivdic) =', rn_rivdic
          WRITE(numout,*) ' river conversion of doc (rn_rivdoc) =', rn_rivdoc
          WRITE(numout,*) ' river conversion of fer (rn_rivfer) =', rn_rivfer
          WRITE(numout,*) ' river conversion of nit (rn_rivnit) =', rn_rivnit
          WRITE(numout,*) ' river conversion of po4 (rn_rivpo4) =', rn_rivpo4
          WRITE(numout,*) ' river conversion of poc (rn_rivpoc) =', rn_rivpoc
          WRITE(numout,*) ' river conversion of sil (rn_rivsil) =', rn_rivsil
!
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natfer'
          WRITE(numout,*) ' '
          WRITE(numout,*) ' solubility of iron in dust (rn_fersol) =', rn_fersol
          WRITE(numout,*) ' coastal release of iron (rn_sedfer) =', rn_sedfer
          WRITE(numout,*) ' scavenging rate of iron (rn_scofer) =', rn_scofer
          WRITE(numout,*) ' minimum scavenging rate of iron (rn_scmfer) =', rn_scmfer
          WRITE(numout,*) ' ligand concentration of iron (rn_ligfer) =', rn_ligfer
!
          WRITE(numout,*) ' '
          WRITE(numout,*) 'natbio'
          WRITE(numout,*) ' '
          WRITE(numout,*)&
     &        ' maximum silification:photosynthesis diatoms rn_bsidia =', rn_bsidia
          WRITE(numout,*)&
     &        ' calcification:photosynthesis coccolithoph. rn_coccal =', rn_coccal
          WRITE(numout,*)&
     &        ' calcite dissolution with coc loss rn_discal  =', rn_discal
          WRITE(numout,*)&
     &        ' aragonite dissolution with pte loss rn_disara=', rn_disara
#    if defined key_trc_foram
          WRITE(numout,*)&
     &        ' calcite dissolution with for loss rn_disfor  =', rn_disfor
          WRITE(numout,*)&
     &        ' calcification:grazing foraminifers rn_forcal =', rn_forcal
#    endif
          WRITE(numout,*)&
     &        ' calcification:grazing pteropods    rn_pteara =', rn_pteara
          WRITE(numout,*)&
     &        ' subsat. aragonite dissolution rate rn_lysara =', rn_lysara
          WRITE(numout,*)&
     &        ' subsat. calcite dissolution rate   rn_lyscal =', rn_lyscal
        DO jm = 1, jpzft
          WRITE(numout,*) ctrcnm(jpbac+jm)
          WRITE(numout,*)&
     &        ' maximal grazing rate               rn_grazoo =', rn_grazoo(jm)
          WRITE(numout,*)&
     &        ' half saturation grazing            rn_grkzoo =', rn_grkzoo(jm)
          WRITE(numout,*)&
     &        ' temperature dependance grazingrate rn_mutpft =', rn_mutpft(jpbac+jm)
          WRITE(numout,*)&
     &        ' temperature dependance grazingrate rn_mudpft =', rn_mudpft(jpbac+jm)
          WRITE(numout,*)&
     &        ' gross growth efficiency            rn_ggezoo =', rn_ggezoo(jm)
          WRITE(numout,*)&
     &        ' unassimilated fraction/fecal pell. rn_unazoo =', rn_unazoo(jm)
           WRITE(numout,*)&
     &        ' losses to POC=0 or GOC=1           nn_sizzoo =', nn_sizzoo(jm)
          WRITE(numout,*)&
     &        ' fraction of respiration as PO4     rn_sigzoo =', rn_sigzoo(jm)
          WRITE(numout,*)&
     &        ' respiration rate of zooplankton    rn_reszoo =', rn_reszoo(jm)
          WRITE(numout,*)&
     &        ' temp. dep. zooplankton respiration rn_retzoo =', rn_retzoo(jm)
          DO jl = jppoc,jpdia+jppft-1
            WRITE(numout,*) &
     &        ' preference for ',ctrcnm(jl),rn_prfzoo(jm,jl)
          END DO
        END DO
          WRITE(numout,*)&
     &        ' minimum food for grazing rn_gramin =', rn_gramin
          WRITE(numout,*)&
     &        ' macrozooplankton mortality rate            =', rn_mormac
          WRITE(numout,*)&
     &        ' temp. dep. macrozooplankton mortality rate =', rn_motmac
          WRITE(numout,*)&
     &        ' gelatinouszooplankton mortality rate       =', rn_morgel
          WRITE(numout,*)&
     &        ' temp. dep. gelatinouszooplankton mortality rate =', rn_motgel
          WRITE(numout,*) ' macrozoo growth enhancement under ice (rn_icemac)=', rn_icemac
          WRITE(numout,*) ' fraction of macrozoo transported (rn_trnmac)=', rn_trnmac

          WRITE(numout,*) ' maximum growth rate bacteria (rn_grabac)   =',rn_grabac
          WRITE(numout,*) ' temperature depend. growth bacteria (rn_mutpft) =',rn_mutpft(jpbac)
          WRITE(numout,*) ' temperature depend. growth bacteria (rn_mudpft) =',rn_mudpft(jpbac)
          WRITE(numout,*) ' Fe half-sat. of bacteria (rn_kmfbac)       =', rn_kmfbac
          WRITE(numout,*) ' DOC half-sat. of bacteria (rn_kmobac)      =', rn_kmobac
          WRITE(numout,*) ' PO4 half-sat. of bacteria (rn_kmpbac)      =', rn_kmpbac
          WRITE(numout,*) ' bacteria preference for DOC, GOC, GON, POC =', rn_gbadoc,rn_gbagoc,rn_gbagon,rn_gbapoc
          WRITE(numout,*) ' bacteria growth efficiency and T depend.   =', rn_ggebac,rn_ggtbac
          WRITE(numout,*) ' bacteria respiration rate and T depend.    =', rn_resbac,rn_retbac
          WRITE(numout,*) ' dnsmin                          =', dnsmin
          WRITE(numout,*) ' silicate half saturation constant diatoms  rn_sildia =', rn_sildia
          WRITE(numout,*) ' iron limited multiplier of Si:C            rn_ferbsi =', rn_ferbsi
          WRITE(numout,*) ' silicate saturated multiplier of Si:C      rn_silbsi =', rn_silbsi
          WRITE(numout,*) ' silicate half saturation constant for Si:C rn_kmsbsi =', rn_kmsbsi
!
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*) ' nitrogen half saturation concentrat.   =',rn_kmnphy(jl)
          WRITE(numout,*) ' phosphate half saturation concentration =',rn_kmpphy(jl)
          WRITE(numout,*) ' phyto respiration as a fraction of growth rn_resphy =', rn_resphy(jl)
          WRITE(numout,*) ' excretion ratio                     =',rn_docphy(jl)
          WRITE(numout,*) ' maximum DOC excretion ratio (rn_domphy) =',rn_domphy(jl)
          WRITE(numout,*) ' maximum growth rate                 =',rn_mumpft(jl)
          WRITE(numout,*) ' temperature dependance growth rate  =',rn_mutpft(jl)
          WRITE(numout,*) ' temperature width growth rate mudpft=',rn_mudpft(jl)
        END DO
        WRITE(numout,*) ' fraction of growth rate during N2fix =',rn_munfix
        WRITE(numout,*) ' exponent for calcite dissolution rate (rn_lyoco3) =',rn_lyoco3
        WRITE(numout,*) ' minimum concentration for bacteria in remineralisation (rn_rembac) =',rn_rembac
        WRITE(numout,*) ' threshold for nutrient uptake for P, N, and Si rn_nutthe =',rn_nutthe
!
! Fe quota model
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natquo'
        WRITE(numout,*) ' '
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*)' iron half saturation conc. (rn_kmfphy)        =',rn_kmfphy(jl)
          WRITE(numout,*)' threshold for nutrient uptake for Fe (rn_rhfphy) =', rn_rhfphy(jl)
          WRITE(numout,*)' maximum quota for Fe (rn_qmaphy)        =',rn_qmaphy(jl)
          WRITE(numout,*)' minimum quota for Fe (rn_qmiphy)        =',rn_qmiphy(jl)
          WRITE(numout,*)' optimal quota for Fe (rn_qopphy)        =',rn_qopphy(jl)
        ENDDO
!
! Light model
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natlit'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' red light absorption coeff. of water       =', rn_ekwred
        WRITE(numout,*) ' green light absorption coeff. of water     =', rn_ekwgrn
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*) ' initial slope PI curve rn_alpphy       =',rn_alpphy(jl)
          WRITE(numout,*) ' maximum chl:C ratio rn_thmphy          =',rn_thmphy(jl)
          WRITE(numout,*) ' light absorption in the red rn_krdphy  =',rn_krdphy(jl)
          WRITE(numout,*) ' light absorption in the blue-green rn_kgrphy =',rn_kgrphy(jl)
          WRITE(numout,*) ' threshold by phytopl. for light rn_tliphy =', rn_tliphy(jl)
        ENDDO
#    if defined key_trc_dms
!
! dms
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natdms'
        WRITE(numout,*) ' '
        WRITE(numout,*)&
     &       ' Ratio of DMS/DMSP released during grazing =',rn_rdddms
        WRITE(numout,*)&
     &       ' DMSP-lyase DMSPd cleavage rate =',rn_xcldmd
        WRITE(numout,*)&
     &       ' microbial yield for DMS production =',rn_dmsyld
        WRITE(numout,*)&
     &       ' half saturation constant for bacteria on DMSPd =',rn_xkdmd
        WRITE(numout,*)&
     &       ' half saturation constant for bacteria on DMS =',rn_xkdms
        DO jm = 1, jpzft
          WRITE(numout,*)&
     &        ' % of DMSP in fecal pellets lost    rn_assdms =', rn_assdms(jm)
        END DO
        WRITE(numout,*)&
     &       ' maximum surface insolation =',rn_etomax
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*)&
     &        ' mean DMSP/C cell ratio rn_rphdmd =', rn_rphdmd(jl)
          WRITE(numout,*)&
     &        ' DMSP leakage coefficient rn_xpldmd =', rn_rphdmd(jl)
        ENDDO
#    endif
#    if defined key_trc_n2o
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natn2o'
        WRITE(numout,*) ' '
        WRITE(numout,*)&
     &       'First depth with N2S production              nn_deun2s =',nn_deun2s
        WRITE(numout,*)&
     &       'oxic N2O/AOU production ratio                rn_aoun2s =',rn_aoun2s
        WRITE(numout,*)&
     &       'Atmospheric pN2O if the date is not in atmco2.dat rn_atmn2o =',rn_atmn2o
        WRITE(numout,*)&
     &       'hypoxic N2O/AOU production ratio             rn_betn2s =',rn_betn2s
        WRITE(numout,*)&
     &       'Exponent of decreasing N2O production with O2 rn_decn2s =',rn_decn2s
        WRITE(numout,*)&
     &       'O2 concentration below which prod. stops increasing rn_omxn2s =',rn_omxn2s
#    endif
#    if defined key_trc_ch4
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natch4'
        WRITE(numout,*) ' '
        WRITE(numout,*)&
     &       'CH4 production ratios               rn_proch? =',rn_proch1,rn_proch2,rn_proch3,rn_proch4,rn_proch5
#    endif
        CALL FLUSH(numout)
      ENDIF

   END SUBROUTINE trc_nam_planktom

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PlankTOM bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_planktom                      ! Empty routine
   END  SUBROUTINE  trc_nam_planktom
#endif  

   !!======================================================================
END MODULE trcnam_planktom
