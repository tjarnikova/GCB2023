MODULE par_planktom
!!======================================================================
!!
!!                         PARAMETER SMS
!!                       *******************************
!!
!!  purpose :
!!  ---------
!!     PARAMETER FILE for biogeochemical model
!!
!!======================================================================
!!  TOP 1.0,  LOCEAN-IPSL (2005)
!! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
!!----------------------------------------------------------------------
#if defined key_planktom
   IMPLICIT NONE
   PUBLIC
!!
!!    ASSIGN AN INTEGER TO NAME INDIVIDUAL TRACERS.
!!    FOR THE DGOM, THE TRACERS SHOULD BE IN THE ORDER:
!!    chemical tracers, detritus, zooplankton,
!!    nutrients, phytoplankton (C,Fe,Chl,Si), isotopes
!!    For documentation of variables see:
!!    http://lgmacweb.env.uea.ac.uk/green_ocean/model/code_description/var_dgom.html
!!----------------------------------------------------------------------
      LOGICAL, PUBLIC, PARAMETER ::   lk_planktom     = .TRUE.  !: PlankTOM flag 
      LOGICAL, PUBLIC, PARAMETER ::   lk_kriest       = .FALSE. !: Kriest flag 
      INTEGER, PUBLIC, PARAMETER :: jpzft=5
      INTEGER, PUBLIC, PARAMETER :: jppft=6,jpfmx=jpzft*(3+jpzft+jppft)
      INTEGER, PUBLIC, PARAMETER :: jptal=1,jpoxy=2,jpdic=3
#  if defined key_trc_piic
      INTEGER, PUBLIC, PARAMETER :: jppiic=4,jpdin=5
#  else
      INTEGER, PUBLIC, PARAMETER :: jpdin=4
#  endif
      INTEGER, PUBLIC, PARAMETER :: jpsil=jpdin+1,jppo4=jpdin+2
      INTEGER, PUBLIC, PARAMETER :: jpfer=jpdin+3,jpdoc=jpdin+4,jpdsi=jpdin+5,jpcal=jpdin+6
      INTEGER, PUBLIC, PARAMETER :: jpara=jpdin+7,jpgon=jpdin+8,jpsfe=jpdin+9
      INTEGER, PUBLIC, PARAMETER :: jpbfe=jpsfe+1,jppoc=jpsfe+2,jpgoc=jpsfe+3,jpbac=jpsfe+4
      INTEGER, PUBLIC, PARAMETER :: jpmic=jpsfe+5
#    if defined key_trc_foram
      INTEGER, PUBLIC, PARAMETER :: jpfor=jpsfe+6,jppte=jpsfe+7,jpmes=jpsfe+8,jpmac=jpsfe+9
#    else
      INTEGER, PUBLIC, PARAMETER :: jppte=jpsfe+6,jpmes=jpsfe+7,jpgel=jpsfe+8,jpmac=jpsfe+9
#    endif
      INTEGER, PUBLIC, PARAMETER :: jpdia=jpmac+1 ,jpmix=jpmac+2 ,jpcoc=jpmac+3 ,jppic=jpmac+4 ,jppha=jpmac+5 ,jpfix=jpmac+6
      INTEGER, PUBLIC, PARAMETER :: jpdfe=jpmac+7 ,jpnfe=jpmac+8 ,jpcfe=jpmac+9 ,jppfe=jpmac+10,jphfe=jpmac+11,jpffe=jpmac+12
      INTEGER, PUBLIC, PARAMETER :: jpdch=jpmac+13,jpnch=jpmac+14,jpcch=jpmac+15,jppch=jpmac+16,jphch=jpmac+17,jpfch=jpmac+18
      INTEGER, PUBLIC, PARAMETER :: jpbsi=jpmac+19
#  if defined key_trc_cfc11
      INTEGER, PUBLIC, PARAMETER :: jpc11=jpmac+20
#  else
      INTEGER, PUBLIC, PARAMETER :: jpc11=jpmac+19
#  endif
#  if defined key_trc_dms
      INTEGER, PUBLIC, PARAMETER :: jpdms=jpc11+1,jpdmd=jpc11+2
#  else
      INTEGER, PUBLIC, PARAMETER :: jpdmd=jpc11
#  endif
#  if defined key_trc_n2o
!      INTEGER, PUBLIC, PARAMETER :: jpn2o=jpdmd+1,jpn2s=jpdmd+2
      INTEGER, PUBLIC, PARAMETER :: jpn2s=jpdmd+1
#  else
      INTEGER, PUBLIC, PARAMETER :: jpn2s=jpdmd
#  endif
#  if defined key_trc_ch4
      INTEGER, PUBLIC, PARAMETER :: jpch1 =jpn2s+1 ,jpch2 =jpn2s+2 ,jpch3 =jpn2s+3 ,jpch4 =jpn2s+4 ,jpch5 =jpn2s+5
      INTEGER, PUBLIC, PARAMETER :: jpch6 =jpn2s+6 ,jpch7 =jpn2s+7 ,jpch8 =jpn2s+8 ,jpch9 =jpn2s+9 ,jpch10=jpn2s+10
      INTEGER, PUBLIC, PARAMETER :: jpch11=jpn2s+11,jpch12=jpn2s+12,jpch13=jpn2s+13,jpch14=jpn2s+14,jpch15=jpn2s+15
      INTEGER, PUBLIC, PARAMETER :: jpch16=jpn2s+16,jpch17=jpn2s+17,jpch18=jpn2s+18,jpch19=jpn2s+19,jpch20=jpn2s+20
      INTEGER, PUBLIC, PARAMETER :: jpch21=jpn2s+21,jpch22=jpn2s+22,jpch23=jpn2s+23,jpch24=jpn2s+24,jpch25=jpn2s+25
      INTEGER, PUBLIC, PARAMETER :: jp_planktom=jpn2s+25
#  else
      INTEGER, PUBLIC, PARAMETER :: jp_planktom=jpn2s
#  endif
#else
      LOGICAL, PUBLIC, PARAMETER ::   lk_planktom     = .FALSE.
      INTEGER, PUBLIC, PARAMETER ::   jp_planktom     =  0
#endif
END MODULE par_planktom
