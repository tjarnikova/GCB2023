MODULE p2zbio
   !!======================================================================
   !!                         ***  MODULE p2zbio  ***
   !! TOP :   LOBSTER
   !!======================================================================
   !! History :    -   !  1999-07  (M. Levy) Original code
   !!              -   !  2000-12  (E. Kestenare) assign a parameter to name individual tracers
   !!              -   !  2001-03  (M. Levy)  LNO3 + dia2d 
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces_reduced'                                     LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   p2z_bio        :  
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc             ! 
   USE sms_pisces
   USE p2zopt
   USE lbclnk          ! 
   USE prtctl_trc      ! Print control for debbuging
   USE trd_oce
   USE trdtrc
   USE iom
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_bio         ! called in ???
   PUBLIC   p2z_bio_init    ! called in ???

   REAL(wp) ::   tmumax     ! maximal phytoplankton growth rate            [s-1]
   REAL(wp) ::   rgamma     ! phytoplankton exudation fraction             [%]
   REAL(wp) ::   fphylab    ! NH4 fraction of phytoplankton exsudation
   REAL(wp) ::   tmminp     ! minimal phytoplancton mortality rate         [0.05/86400 s-1=20 days]
   REAL(wp) ::   aki        ! light photosynthesis half saturation constant[W/m2]
   !
   REAL(wp) ::   akno3      ! nitrate limitation half-saturation value     [mmol/m3]
   REAL(wp) ::   aknh4      ! ammonium limitation half-saturation value    [mmol/m3]
   REAL(wp) ::   taunn      ! nitrification rate                           [s-1]
   REAL(wp) ::   psinut     ! inhibition of nitrate uptake by ammonium
   !
   REAL(wp) ::   taudn      ! detritus breakdown rate                        [0.1/86400 s-1=10 days]
   REAL(wp) ::   fdetlab    ! NH4 fraction of detritus dissolution
   !
   REAL(wp) ::   taudomn    ! DOM breakdown rate                             [s-1]
   !                        ! slow remineralization rate of semi-labile dom to nh4 (1 month)
   !
   REAL(wp) ::   rppz       ! ivlev coeff for zoo mortality
   REAL(wp) ::   taus       ! specific zooplankton maximal grazing rate              [s-1]
   !                                     ! 0.75/86400 s-1=8.680555E-6    1/86400 = 1.15e-5
   REAL(wp) ::   aks        ! half-saturation constant for total zooplankton grazing [mmolN.m-3]
   REAL(wp) ::   rpnaz      ! non-assimilated phytoplankton by zooplancton           [%]
   REAL(wp) ::   rdnaz      ! non-assimilated detritus by zooplankton                [%]
   REAL(wp) ::   tauzn      ! zooplancton specific excretion rate                    [0.1/86400 s-1=10 days]
   REAL(wp) ::   tmminz     ! minimal zooplankton mortality rate                     [(mmolN/m3)-1 d-1]
   REAL(wp) ::   fzoolab    ! NH4 fraction of zooplankton excretion
   REAL(wp) ::   fdbod      ! zooplankton mortality fraction that goes to detritus

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p2zbio.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p2z_bio( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_bio  ***
      !!
      !! ** Purpose :   compute the now trend due to biogeochemical processes
      !!              and add it to the general trend of passive tracers equations
      !!
      !! ** Method  :   each now biological flux is calculated in function of now
      !!              concentrations of tracers.
      !!              depending on the tracer, these fluxes are sources or sinks.
      !!              the total of the sources and sinks for each tracer
      !!              is added to the general trend.
      !!        
      !!                      tra = tra + zf...tra - zftra...
      !!                                     |         |
      !!                                     |         |
      !!                                  source      sink
      !!        
      !!              IF 'key_diabio' defined , the biogeochemical trends
      !!              for passive tracers are saved for futher diagnostics.
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER  ::   ji, jj, jk, jl
      REAL(wp) ::   zdet, zzoo, zphy, zno3, znh4, zdom      ! now concentrations
      REAL(wp) ::   zlno3, zlnh4, zle, zlt                  ! limitation terms for phyto
      REAL(wp) ::   zno3phy, znh4phy, zphynh4, zphydom
      REAL(wp) ::   zphydet, zphyzoo, zdetzoo
      REAL(wp) ::   zzoonh4, zzoodom, zzoodet, zdetnh4, zdetdom
      REAL(wp) ::   znh4no3, zdomnh4, zppz, zpdz, zpppz, zppdz, zfood
      REAL(wp) ::   zfilpz, zfildz, zphya, zzooa, zno3a
      REAL(wp) ::   znh4a, zdeta, zdoma, zzoobod, zboddet, zdomaju
      REAL(wp) ::   ze3t
      REAL(wp), POINTER,   DIMENSION(:,:,:) :: zw2d
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: zw3d
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p2z_bio')
      !
      IF( ln_diatrc .OR. lk_iomput ) THEN
         CALL wrk_alloc( jpi, jpj,     17, zw2d )
         CALL wrk_alloc( jpi, jpj, jpk, 3, zw3d )
      ENDIF

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' p2z_bio: LOBSTER bio-model'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
      ENDIF

      xksi(:,:) = 0.e0        ! zooplakton closure ( fbod)
      IF( ln_diatrc .OR. lk_iomput ) THEN
         zw2d  (:,:,:) = 0.e0
         zw3d(:,:,:,:) = 0.e0
      ENDIF

      !                                      ! -------------------------- !
      DO jk = 1, jpkbm1                      !  Upper ocean (bio-layers)  !
         !                                   ! -------------------------- !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ! trophic variables( det, zoo, phy, no3, nh4, dom)
               ! ------------------------------------------------

               ! negative trophic variables DO not contribute to the fluxes
               zdet = MAX( 0.e0, trn(ji,jj,jk,jpdet) )
               zzoo = MAX( 0.e0, trn(ji,jj,jk,jpzoo) )
               zphy = MAX( 0.e0, trn(ji,jj,jk,jpphy) )
               zno3 = MAX( 0.e0, trn(ji,jj,jk,jpno3) )
               znh4 = MAX( 0.e0, trn(ji,jj,jk,jpnh4) )
               zdom = MAX( 0.e0, trn(ji,jj,jk,jpdom) )

               ! Limitations
               zlt   = 1.
               zle   = 1. - EXP( -etot(ji,jj,jk) / aki / zlt )
               ! psinut,akno3,aknh4 added by asklod AS Kremeur 2005-03
               zlno3 = zno3 * EXP( -psinut * znh4 ) / ( akno3 + zno3 )
               zlnh4 = znh4 / (znh4+aknh4)  

               ! sinks and sources
               !    phytoplankton production and exsudation
               zno3phy = tmumax * zle * zlt * zlno3 * zphy
               znh4phy = tmumax * zle * zlt * zlnh4 * zphy

               !    fphylab added by asklod AS Kremeur 2005-03
               zphydom = rgamma * (1 - fphylab) * (zno3phy + znh4phy)
               zphynh4 = rgamma * fphylab * (zno3phy + znh4phy)
               ! zooplankton production
               !    preferences
               zppz = rppz
               zpdz = 1. - rppz
               zpppz = ( zppz * zphy ) / ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
               zppdz = ( zpdz * zdet ) / ( ( zppz * zphy + zpdz * zdet ) + 1.e-13 )
               zfood = zpppz * zphy + zppdz * zdet
               !    filtration 
               zfilpz = taus * zpppz / (aks + zfood)
               zfildz = taus * zppdz / (aks + zfood)
               !    grazing
               zphyzoo = zfilpz * zphy * zzoo
               zdetzoo = zfildz * zdet * zzoo

               ! fecal pellets production
               zzoodet = rpnaz * zphyzoo + rdnaz * zdetzoo

               ! zooplankton liquide excretion
               zzoonh4 = tauzn * fzoolab * zzoo  
               zzoodom = tauzn * (1 - fzoolab) * zzoo

               ! mortality
               !    phytoplankton mortality
               zphydet = tmminp * zphy

               !    zooplankton mortality
               !    closure : flux grazing is redistributed below level jpkbio
               zzoobod = tmminz * zzoo * zzoo
               xksi(ji,jj) = xksi(ji,jj) + (1-fdbod) * zzoobod * fse3t(ji,jj,jk)
               zboddet = fdbod * zzoobod

               ! detritus and dom breakdown
               zdetnh4 = taudn * fdetlab * zdet
               zdetdom = taudn * (1 - fdetlab) * zdet

               zdomnh4 = taudomn * zdom

               ! flux added to express how the excess of nitrogen from 
               ! PHY, ZOO and DET to DOM goes directly to NH4 (flux of ajustment)
               zdomaju = (1 - redf/reddom) * (zphydom + zzoodom + zdetdom)

               ! Nitrification 
               znh4no3 = taunn * znh4

               ! determination of trends
               !    total trend for each biological tracer
               zphya =   zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet
               zzooa =   zphyzoo + zdetzoo - zzoodet - zzoodom - zzoonh4 - zzoobod
               zno3a = - zno3phy + znh4no3
               znh4a = - znh4phy - znh4no3 + zphynh4 + zzoonh4 + zdomnh4 + zdetnh4 + zdomaju
               zdeta =   zphydet + zzoodet - zdetzoo - zdetnh4 - zdetdom + zboddet
               zdoma =   zphydom + zzoodom + zdetdom - zdomnh4 - zdomaju

               ! tracer flux at totox-point added to the general trend
               tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + zdeta
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zzooa
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zphya
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zno3a
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + znh4a
               tra(ji,jj,jk,jpdom) = tra(ji,jj,jk,jpdom) + zdoma


               IF( ( ln_diabio .AND. .NOT. lk_iomput ) .OR. l_trdtrc ) THEN
                  trbio(ji,jj,jk,jp_pcs0_trd     ) = zno3phy
                  trbio(ji,jj,jk,jp_pcs0_trd +  1) = znh4phy
                  trbio(ji,jj,jk,jp_pcs0_trd +  2) = zphynh4
                  trbio(ji,jj,jk,jp_pcs0_trd +  3) = zphydom
                  trbio(ji,jj,jk,jp_pcs0_trd +  4) = zphyzoo
                  trbio(ji,jj,jk,jp_pcs0_trd +  5) = zphydet
                  trbio(ji,jj,jk,jp_pcs0_trd +  6) = zdetzoo
                  !  trend number 8 in p2zsed
                  trbio(ji,jj,jk,jp_pcs0_trd +  8) = zzoodet
                  trbio(ji,jj,jk,jp_pcs0_trd +  9) = zzoobod
                  trbio(ji,jj,jk,jp_pcs0_trd + 10) = zzoonh4
                  trbio(ji,jj,jk,jp_pcs0_trd + 11) = zzoodom
                  trbio(ji,jj,jk,jp_pcs0_trd + 12) = znh4no3
                  trbio(ji,jj,jk,jp_pcs0_trd + 13) = zdomnh4
                  trbio(ji,jj,jk,jp_pcs0_trd + 14) = zdetnh4
                  trbio(ji,jj,jk,jp_pcs0_trd + 15) = zdetdom
                  !  trend number 17 in p2zexp
                ENDIF
                IF( ln_diatrc .OR. lk_iomput ) THEN
                  ! convert fluxes in per day
                  ze3t = fse3t(ji,jj,jk) * 86400.
                  zw2d(ji,jj,1)  = zw2d(ji,jj,1)  + zno3phy * ze3t
                  zw2d(ji,jj,2)  = zw2d(ji,jj,2)  + znh4phy * ze3t
                  zw2d(ji,jj,3)  = zw2d(ji,jj,3)  + zphydom * ze3t
                  zw2d(ji,jj,4)  = zw2d(ji,jj,4)  + zphynh4 * ze3t
                  zw2d(ji,jj,5)  = zw2d(ji,jj,5)  + zphyzoo * ze3t
                  zw2d(ji,jj,6)  = zw2d(ji,jj,6)  + zphydet * ze3t
                  zw2d(ji,jj,7)  = zw2d(ji,jj,7)  + zdetzoo * ze3t
                  zw2d(ji,jj,8)  = zw2d(ji,jj,8)  + zzoodet * ze3t
                  zw2d(ji,jj,9)  = zw2d(ji,jj,9)  + zzoobod * ze3t
                  zw2d(ji,jj,10) = zw2d(ji,jj,10) + zzoonh4 * ze3t
                  zw2d(ji,jj,11) = zw2d(ji,jj,11) + zzoodom * ze3t
                  zw2d(ji,jj,12) = zw2d(ji,jj,12) + znh4no3 * ze3t
                  zw2d(ji,jj,13) = zw2d(ji,jj,13) + zdomnh4 * ze3t
                  zw2d(ji,jj,14) = zw2d(ji,jj,14) + zdetnh4 * ze3t
                  zw2d(ji,jj,15) = zw2d(ji,jj,15) + ( zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet ) * ze3t
                  zw2d(ji,jj,16) = zw2d(ji,jj,16) + ( zphyzoo + zdetzoo - zzoodet - zzoobod - zzoonh4 - zzoodom ) * ze3t
                  zw2d(ji,jj,17) = zw2d(ji,jj,17) + zdetdom * ze3t
                  !   
                  zw3d(ji,jj,jk,1) = zno3phy * 86400
                  zw3d(ji,jj,jk,2) = znh4phy * 86400     
                  zw3d(ji,jj,jk,3) = znh4no3 * 86400   
                   ! 
                ENDIF
            END DO
         END DO
      END DO

      !                                      ! -------------------------- !
      DO jk = jpkb, jpkm1                    !  Upper ocean (bio-layers)  !
         !                                   ! -------------------------- !
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1 
               ! remineralisation of all quantities towards nitrate 

               !    trophic variables( det, zoo, phy, no3, nh4, dom)
               !       negative trophic variables DO not contribute to the fluxes
               zdet = MAX( 0.e0, trn(ji,jj,jk,jpdet) )
               zzoo = MAX( 0.e0, trn(ji,jj,jk,jpzoo) )
               zphy = MAX( 0.e0, trn(ji,jj,jk,jpphy) )
               zno3 = MAX( 0.e0, trn(ji,jj,jk,jpno3) )
               znh4 = MAX( 0.e0, trn(ji,jj,jk,jpnh4) )
               zdom = MAX( 0.e0, trn(ji,jj,jk,jpdom) )

               !    Limitations
               zlt   = 0.e0
               zle   = 0.e0
               zlno3 = 0.e0
               zlnh4 = 0.e0

               !    sinks and sources
               !       phytoplankton production and exsudation
               zno3phy = 0.e0
               znh4phy = 0.e0
               zphydom = 0.e0
               zphynh4 = 0.e0

               !    zooplankton production
               zphyzoo = 0.e0      ! grazing
               zdetzoo = 0.e0

               zzoodet = 0.e0      ! fecal pellets production

               zzoonh4 = tauzn * fzoolab * zzoo         ! zooplankton liquide excretion
               zzoodom = tauzn * (1 - fzoolab) * zzoo

               !    mortality
               zphydet = tmminp * zphy      ! phytoplankton mortality

               zzoobod = 0.e0               ! zooplankton mortality
               zboddet = 0.e0               ! closure : flux fbod is redistributed below level jpkbio

               !    detritus and dom breakdown
               zdetnh4 = taudn * fdetlab * zdet
               zdetdom = taudn * (1 - fdetlab) * zdet

               zdomnh4 = taudomn * zdom
               zdomaju = (1 - redf/reddom) * (zphydom + zzoodom + zdetdom)

               !    Nitrification
               znh4no3 = taunn * znh4


               ! determination of trends
               !     total trend for each biological tracer
               zphya =   zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet
               zzooa =   zphyzoo + zdetzoo - zzoodet - zzoodom - zzoonh4 - zzoobod
               zno3a = - zno3phy + znh4no3 
               znh4a = - znh4phy - znh4no3 + zphynh4 + zzoonh4 + zdomnh4 + zdetnh4 + zdomaju
               zdeta = zphydet + zzoodet  - zdetzoo - zdetnh4 - zdetdom + zboddet
               zdoma = zphydom + zzoodom + zdetdom - zdomnh4 - zdomaju

               ! tracer flux at totox-point added to the general trend
               tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + zdeta
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zzooa
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zphya
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zno3a
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + znh4a
               tra(ji,jj,jk,jpdom) = tra(ji,jj,jk,jpdom) + zdoma
               !
               IF( ( ln_diabio .AND. .NOT. lk_iomput ) .OR. l_trdtrc ) THEN
                  trbio(ji,jj,jk,jp_pcs0_trd     ) = zno3phy
                  trbio(ji,jj,jk,jp_pcs0_trd +  1) = znh4phy
                  trbio(ji,jj,jk,jp_pcs0_trd +  2) = zphynh4
                  trbio(ji,jj,jk,jp_pcs0_trd +  3) = zphydom
                  trbio(ji,jj,jk,jp_pcs0_trd +  4) = zphyzoo
                  trbio(ji,jj,jk,jp_pcs0_trd +  5) = zphydet
                  trbio(ji,jj,jk,jp_pcs0_trd +  6) = zdetzoo
                  !  trend number 8 in p2zsed
                  trbio(ji,jj,jk,jp_pcs0_trd +  8) = zzoodet
                  trbio(ji,jj,jk,jp_pcs0_trd +  9) = zzoobod
                  trbio(ji,jj,jk,jp_pcs0_trd + 10) = zzoonh4
                  trbio(ji,jj,jk,jp_pcs0_trd + 11) = zzoodom
                  trbio(ji,jj,jk,jp_pcs0_trd + 12) = znh4no3
                  trbio(ji,jj,jk,jp_pcs0_trd + 13) = zdomnh4
                  trbio(ji,jj,jk,jp_pcs0_trd + 14) = zdetnh4
                  trbio(ji,jj,jk,jp_pcs0_trd + 15) = zdetdom
                  !  trend number 17 in p2zexp 
                ENDIF
                IF( ln_diatrc .OR. lk_iomput ) THEN
                  ! convert fluxes in per day
                  ze3t = fse3t(ji,jj,jk) * 86400.
                  zw2d(ji,jj,1)  = zw2d(ji,jj,1)  + zno3phy * ze3t
                  zw2d(ji,jj,2)  = zw2d(ji,jj,2)  + znh4phy * ze3t
                  zw2d(ji,jj,3)  = zw2d(ji,jj,3)  + zphydom * ze3t
                  zw2d(ji,jj,4)  = zw2d(ji,jj,4)  + zphynh4 * ze3t
                  zw2d(ji,jj,5)  = zw2d(ji,jj,5)  + zphyzoo * ze3t
                  zw2d(ji,jj,6)  = zw2d(ji,jj,6)  + zphydet * ze3t
                  zw2d(ji,jj,7)  = zw2d(ji,jj,7)  + zdetzoo * ze3t
                  zw2d(ji,jj,8)  = zw2d(ji,jj,8)  + zzoodet * ze3t
                  zw2d(ji,jj,9)  = zw2d(ji,jj,9)  + zzoobod * ze3t
                  zw2d(ji,jj,10) = zw2d(ji,jj,10) + zzoonh4 * ze3t
                  zw2d(ji,jj,11) = zw2d(ji,jj,11) + zzoodom * ze3t
                  zw2d(ji,jj,12) = zw2d(ji,jj,12) + znh4no3 * ze3t
                  zw2d(ji,jj,13) = zw2d(ji,jj,13) + zdomnh4 * ze3t
                  zw2d(ji,jj,14) = zw2d(ji,jj,14) + zdetnh4 * ze3t
                  zw2d(ji,jj,15) = zw2d(ji,jj,15) + ( zno3phy + znh4phy - zphynh4 - zphydom - zphyzoo - zphydet ) * ze3t
                  zw2d(ji,jj,16) = zw2d(ji,jj,16) + ( zphyzoo + zdetzoo - zzoodet - zzoobod - zzoonh4 - zzoodom ) * ze3t
                  zw2d(ji,jj,17) = zw2d(ji,jj,17) + zdetdom * ze3t
                  !   
                  zw3d(ji,jj,jk,1) = zno3phy * 86400
                  zw3d(ji,jj,jk,2) = znh4phy * 86400
                  zw3d(ji,jj,jk,3) = znh4no3 * 86400
                   !
                ENDIF
            END DO
         END DO
      END DO

      IF( ln_diatrc .OR. lk_iomput ) THEN
         DO jl = 1, 17 
            CALL lbc_lnk( zw2d(:,:,jl),'T', 1. )
         END DO
         DO jl = 1, 3
            CALL lbc_lnk( zw3d(:,:,:,jl),'T', 1. )
         END DO
      ENDIF
      IF( lk_iomput ) THEN
            ! Save diagnostics
        CALL iom_put( "TNO3PHY", zw2d(:,:,1) )
        CALL iom_put( "TNH4PHY", zw2d(:,:,2) )
        CALL iom_put( "TPHYDOM", zw2d(:,:,3) )
        CALL iom_put( "TPHYNH4", zw2d(:,:,4) )
        CALL iom_put( "TPHYZOO", zw2d(:,:,5) )
        CALL iom_put( "TPHYDET", zw2d(:,:,6) )
        CALL iom_put( "TDETZOO", zw2d(:,:,7) )
        CALL iom_put( "TZOODET", zw2d(:,:,8) )
        CALL iom_put( "TZOOBOD", zw2d(:,:,9) )
        CALL iom_put( "TZOONH4", zw2d(:,:,10) )
        CALL iom_put( "TZOODOM", zw2d(:,:,11) )
        CALL iom_put( "TNH4NO3", zw2d(:,:,12) )
        CALL iom_put( "TDOMNH4", zw2d(:,:,13) )
        CALL iom_put( "TDETNH4", zw2d(:,:,14) )
        CALL iom_put( "TPHYTOT", zw2d(:,:,15) )
        CALL iom_put( "TZOOTOT", zw2d(:,:,16) )
         ! 
        CALL iom_put( "FNO3PHY", zw3d(:,:,:,1) )
        CALL iom_put( "FNH4PHY", zw3d(:,:,:,2) )
        CALL iom_put( "FNH4NO3", zw3d(:,:,:,3) )
         !
       ELSE
          IF( ln_diatrc ) THEN
            !
            trc2d(:,:,jp_pcs0_2d    ) = zw2d(:,:,1) 
            trc2d(:,:,jp_pcs0_2d + 1) = zw2d(:,:,2) 
            trc2d(:,:,jp_pcs0_2d + 2) = zw2d(:,:,3) 
            trc2d(:,:,jp_pcs0_2d + 3) = zw2d(:,:,4) 
            trc2d(:,:,jp_pcs0_2d + 4) = zw2d(:,:,5) 
            trc2d(:,:,jp_pcs0_2d + 5) = zw2d(:,:,6) 
            trc2d(:,:,jp_pcs0_2d + 6) = zw2d(:,:,7) 
                     ! trend number 8 is in p2zsed.F
            trc2d(:,:,jp_pcs0_2d +  8) = zw2d(:,:,8) 
            trc2d(:,:,jp_pcs0_2d +  9) = zw2d(:,:,9) 
            trc2d(:,:,jp_pcs0_2d + 10) = zw2d(:,:,10) 
            trc2d(:,:,jp_pcs0_2d + 11) = zw2d(:,:,11) 
            trc2d(:,:,jp_pcs0_2d + 12) = zw2d(:,:,12) 
            trc2d(:,:,jp_pcs0_2d + 13) = zw2d(:,:,13) 
            trc2d(:,:,jp_pcs0_2d + 14) = zw2d(:,:,14) 
            trc2d(:,:,jp_pcs0_2d + 15) = zw2d(:,:,15) 
            trc2d(:,:,jp_pcs0_2d + 16) = zw2d(:,:,16) 
            trc2d(:,:,jp_pcs0_2d + 17) = zw2d(:,:,17) 
            ! trend number 19 is in p2zexp.F
            trc3d(:,:,:,jp_pcs0_3d    ) = zw3d(:,:,:,1) 
            trc3d(:,:,:,jp_pcs0_3d + 1) = zw3d(:,:,:,2) 
            trc3d(:,:,:,jp_pcs0_3d + 2) = zw3d(:,:,:,3) 
         ENDIF
        !
      ENDIF

      IF( ln_diabio .AND. .NOT. lk_iomput )  THEN
         DO jl = jp_pcs0_trd, jp_pcs1_trd
            CALL lbc_lnk( trbio(:,:,1,jl),'T', 1. )
         END DO 
      ENDIF
      !
      IF( l_trdtrc ) THEN
         DO jl = jp_pcs0_trd, jp_pcs1_trd
            CALL trd_trc( trbio(:,:,:,jl), jl, kt )   ! handle the trend
         END DO
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_diatrc .OR. lk_iomput ) THEN
         CALL wrk_dealloc( jpi, jpj,     17, zw2d )
         CALL wrk_dealloc( jpi, jpj, jpk, 3, zw3d )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p2z_bio')
      !
   END SUBROUTINE p2z_bio

   SUBROUTINE p2z_bio_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_bio_init  ***
      !!
      !! ** Purpose :  biological parameters 
      !!
      !! ** Method  :   Read namelist and check the parameters
      !!
      !!----------------------------------------------------------------------
      NAMELIST/namlobphy/ tmumax, rgamma, fphylab, tmminp, aki
      NAMELIST/namlobnut/ akno3, aknh4, taunn, psinut
      NAMELIST/namlobzoo/ rppz, taus, aks, rpnaz, rdnaz, tauzn, fzoolab, fdbod, tmminz
      NAMELIST/namlobdet/  taudn, fdetlab
      NAMELIST/namlobdom/ taudomn
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist namlobphy in reference namelist : Lobster biological parameters
      READ  ( numnatp_ref, namlobphy, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobphy in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namlobphy in configuration namelist : Lobster biological parameters
      READ  ( numnatp_cfg, namlobphy, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobphy in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namlobphy )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobphy'
          WRITE(numout,*) '    phyto max growth rate                                tmumax    =', 86400 * tmumax, ' d'
          WRITE(numout,*) '    phytoplankton exudation fraction                     rgamma    =', rgamma
          WRITE(numout,*) '    NH4 fraction of phytoplankton exsudation             fphylab   =', fphylab
          WRITE(numout,*) '    minimal phyto mortality rate                         tmminp    =', 86400 * tmminp
          WRITE(numout,*) '    light hlaf saturation constant                       aki       =', aki
          WRITE(numout,*) ' '
      ENDIF

      REWIND( numnatp_ref )              ! Namelist namlobnut in reference namelist : Lobster nutriments parameters
      READ  ( numnatp_ref, namlobnut, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobnut in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namlobnut in configuration namelist : Lobster nutriments parameters
      READ  ( numnatp_cfg, namlobnut, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobnut in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namlobnut )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobnut'
          WRITE(numout,*) '    half-saturation nutrient for no3 uptake              akno3     =', akno3
          WRITE(numout,*) '    half-saturation nutrient for nh4 uptake              aknh4     =', aknh4
          WRITE(numout,*) '    nitrification rate                                   taunn     =', taunn
          WRITE(numout,*) '    inhibition of no3 uptake by nh4                      psinut    =', psinut
          WRITE(numout,*) ' '
      ENDIF

      REWIND( numnatp_ref )              ! Namelist namlobzoo in reference namelist : Lobster zooplankton parameters
      READ  ( numnatp_ref, namlobzoo, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobzoo in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namlobzoo in configuration namelist : Lobster zooplankton parameters
      READ  ( numnatp_cfg, namlobzoo, IOSTAT = ios, ERR = 906 )
906   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobzoo in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namlobzoo )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobzoo'
          WRITE(numout,*) '    zoo preference for phyto                             rppz      =', rppz
          WRITE(numout,*) '    maximal zoo grazing rate                             taus      =', 86400 * taus, ' d'
          WRITE(numout,*) '    half saturation constant for zoo food                aks       =', aks
          WRITE(numout,*) '    non-assimilated phyto by zoo                         rpnaz     =', rpnaz
          WRITE(numout,*) '    non-assimilated detritus by zoo                      rdnaz     =', rdnaz
          WRITE(numout,*) '    zoo specific excretion rate                          tauzn     =', 86400 * tauzn
          WRITE(numout,*) '    minimal zoo mortality rate                           tmminz    =', 86400 * tmminz
          WRITE(numout,*) '    NH4 fraction of zooplankton excretion                fzoolab   =', fzoolab
          WRITE(numout,*) '    Zooplankton mortality fraction that goes to detritus fdbod     =', fdbod
          WRITE(numout,*) ' '
      ENDIF

      REWIND( numnatp_ref )              ! Namelist namlobdet in reference namelist : Lobster detritus parameters
      READ  ( numnatp_ref, namlobdet, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobdet in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namlobdet in configuration namelist : Lobster detritus parameters
      READ  ( numnatp_cfg, namlobdet, IOSTAT = ios, ERR = 908 )
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobdet in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namlobdet )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobdet'
          WRITE(numout,*) '    detrital breakdown rate                              taudn     =', 86400 * taudn , ' d'
          WRITE(numout,*) '    NH4 fraction of detritus dissolution                 fdetlab   =', fdetlab
          WRITE(numout,*) ' '
      ENDIF

      REWIND( numnatp_ref )              ! Namelist namlobdom in reference namelist : Lobster DOM breakdown rate
      READ  ( numnatp_ref, namlobdom, IOSTAT = ios, ERR = 909)
909   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobdom in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist namlobdom in configuration namelist : Lobster DOM breakdown rate
      READ  ( numnatp_cfg, namlobdom, IOSTAT = ios, ERR = 910 )
910   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlobdom in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, namlobdom )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobdom'
          WRITE(numout,*) '    DOM breakdown rate                                 taudomn     =', 86400 * taudn , ' d'
          WRITE(numout,*) ' '
      ENDIF
      !
   END SUBROUTINE p2z_bio_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_bio( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_bio: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_bio
#endif 

   !!======================================================================
END MODULE  p2zbio
