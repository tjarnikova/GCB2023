#!/bin/sh
#SBATCH --mail-type=NONE
#SBATCH -q short-eth
#SBATCH -o ptom-totals.log
#SBATCH -e ptom-totals_err.log
. /etc/profile
module list
set -x
cd ~/scratch/$ENAM
ln -fs ~/Programs/basin_mask_v3.6.nc basin_mask.nc
ln -fs ~/Programs/bottom_volume_v3.6.nc bottom_volume.nc
echo $VERSION > total.arg
echo ${SIM} >> total.arg
echo ${YEARN} >> total.arg
if [ -f total.arg ]
 then ~/nemo_v3.6/${VERSION}/CONFIG/ORCA2_LIM_PlankTOM/EXP00/totalTOM12
fi
sed -e s/1845/$YEARN/g ~/Programs/total.DOCTRP > total.inp
~/Programs/totalDOCTRP
sed -e s/1750/$YEARN/g ~/Programs/total.PICflx > total.inp
if [ -f totals.PICflx ]
 then ~/Programs/totalanything_units |tail -n 1 >> totals.PICflx
 else ~/Programs/totalanything_units >> totals.PICflx
fi
sed -e s/1948/$YEARN/g ~/Programs/total.Siflx > total.inp
if [ -f totals.Siflx ]
 then ~/Programs/totalanything_units |tail -n 1 >> totals.Siflx
 else ~/Programs/totalanything_units >> totals.Siflx
fi
sed -e s/00021900/$timestep/g ~/Programs/total.Sirestart > total.inp
if [ -f totals.Sirestart ]
 then ~/Programs/totalanything_units |tail -n 1 >> totals.Sirestart
 else ~/Programs/totalanything_units >> totals.Sirestart
fi
sed -e s/2020/$YEARN/g ~/Programs/total.Siinv > total.inp
if [ -f totals.Siinv ]
 then ~/Programs/totalanything_units |tail -n 1 >> totals.Siinv
 else ~/Programs/totalanything_units >> totals.Siinv
fi
if [ $YEARN -ge 1970 ];then
  sed -e s/1967/$YEARN/g ../HPC_ET_HI04/totald.inp > total.inp
  ~/Programs/totalanymonth |tail -n 12 >> total_Cflu.monthly
#  ncra -v Cflx ORCA2_1m_${YEARP}*_diad_T.nc ORCA2_1y_${YEARP}_Cflx.nc
fi

if [ $YEARN -eq 2012 ];then
module add gcc hdf5/1.10.6/gcc netcdf/4.7.4/gcc nco/4.9.3/intel
  ncrcat -v nav_lon,nav_lat,C11 ORCA2_1m_198[123456789]0101_198?1231_ptrc_T.nc ORCA2_1m_199?0101_199?1231_ptrc_T.nc ORCA2_1m_200?0101_2???1231_ptrc_T.nc ORCA2_1m_201[012]0101_2???1231_ptrc_T.nc ${SIM}_1981_2012_C11.nc
module purge
module add gcc/9.2.0 netcdf/4.7.4/parallel/gcc-openmpi hdf5/1.10.6/gcc-openmpi java/jdk1.8.0_231 mpi/openmpi/4.0.3/gcc/ib perl/5.30.2 ferret/7.5.0
  sed -e s/TOM12_ET_PID4/$ENAM/ -e s/PID4/$SIM/ ../TOM12_ET_PID4/opa2reg_C11 > opa2reg_C11
  ln -fs opa2reg_C11 opa2reg_3d
  ~/Programs/orca2woa
  ln -fs ${ENAM}_1981_2012_CFC11.nc CFC11.nc
  ~/Programs/aveCFC11
  \rm opa2reg_3d
fi
if [ $YEARP -eq 2019 ];then
  sed s/\ \ \ /\ /g totals3.output > ${ENAM}_totals3.csv
  sed s/\ \ /\ /g ${ENAM}_totals3.csv > tmp
  sed s/\ /,/g tmp > ${ENAM}_totals3.csv
  module add gcc hdf5/1.10.6/gcc netcdf/4.7.4/gcc nco/4.9.3/intel
  ncrcat -v nav_lon,nav_lat,pCO2 ORCA2_1m_19[789]?0101_19??1231_diad_T.nc ORCA2_1m_2???0101_2???1231_diad_T.nc ${SIM}_1970_2019_pCO2.nc
  module purge
  module add gcc/9.2.0 netcdf/4.7.4/parallel/gcc-openmpi hdf5/1.10.6/gcc-openmpi java/jdk1.8.0_231 mpi/openmpi/4.0.3/gcc/ib perl/5.30.2 ferret/7.5.0
  sed -e s/TOM12_ET_PIGT/${ENAM}/ -e s/PIGT/$SIM/ ~/Programs/ORCA2WOA/opa2reg_spco2 > opa2reg_spco2
  ln -fs opa2reg_spco2 opa2reg_2d
  ~/Programs/orca2woa_ada
  sed  -e s/TOM12_ET_PIGT/$ENAM/ -e s/PIGT/$SIM/ ~/Programs/ORCA2WOA/opa2reg_transpose > opa2reg_transpose
  ln -s opa2reg_transpose opa2reg_in
  /gpfs/env/e031/Programs/ORCA2WOA/Surface/transpose_ada
  ln -fs ${ENAM}_1970_2019_spco2_dateline.nc spco2.nc
  ../SOCOM/rss_monthly_ada_2020 > ${ENAM}_RSS_vs_Cflu
  \rm opa2reg_2d
  sed -e s/\ \ \ \ \ /\ /g ${ENAM}_RSS_vs_Cflu > ${ENAM}_RSS_vs_Cflu_land2020.csv
  sed -e s/\ \ \ \ /\ /g ${ENAM}_RSS_vs_Cflu_land2020.csv > ${ENAM}_RSS_vs_Cflu
  sed -e s/\ \ \ /\ /g ${ENAM}_RSS_vs_Cflu > ${ENAM}_RSS_vs_Cflu_land2020.csv
  sed -e s/\ \ /\ /g ${ENAM}_RSS_vs_Cflu_land2020.csv > ${ENAM}_RSS_vs_Cflu
  sed -e s/\ /,/g ${ENAM}_RSS_vs_Cflu > ${ENAM}_RSS_vs_Cflu_land2020.csv
  array=( $(tail -n 1 ${ENAM}_RSS_vs_Cflu) )
  echo ${array[3]}
  ../SOCOM/rss_lat_2020 >> ../SOCOM/rss_lat.row_2020
  module add gcc hdf5/1.10.6/gcc netcdf/4.7.4/gcc nco/4.9.3/intel
#  ncrcat ORCA2_1y_????_Cflx.nc ORCA2_1y_1970_2019_Cflx.nc && \rm ORCA2_1y_????_Cflx.nc
#  ln -fs ../TOM12_ET_PIGT/opa2reg_Cfl_ann opa2reg_2d
#  ~/Programs/orca2woa_ada
#  ln -fs ../TOM12_ET_PIGT/opa2reg_Cfl_transp opa2reg_in
#  /gpfs/env/e031/Programs/ORCA2WOA/Surface/transpose_ada
#  ~/Programs/totprovint_2020 >> ../SOCOM/RMSE_SO_2020
#  ~/Programs/totlatCfl_MSE_2020 >> ../SOCOM/RMSE_lat_2020
#  ~/Programs/totallatinteran > totals.Cflx_lat_2020

  ncrcat -v nav_lon,nav_lat,votemper ORCA2_1m_19[6789]?0101_19??1231_grid_T.nc ORCA2_1m_2???0101_2???1231_grid_T.nc ${SIM}_1955_2019_votemper.nc
fi
