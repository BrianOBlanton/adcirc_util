function AdcUtil=adcirc_util_init
echo off

AdcUtil.MaxVarMapping.maxele='zeta_max';
AdcUtil.MaxVarMapping.maxvel='vel_max';
AdcUtil.MaxVarMapping.maxwvel={'u-wind','v-wind'};
AdcUtil.MaxVarMapping.initiallydry='initiallydry';
AdcUtil.MaxVarMapping.maxinundepth='inun_max';

AdcUtil.MaxVarNaN.maxele=-1000;
AdcUtil.MaxVarNaN.maxvel=eps*1e6;
AdcUtil.MaxVarNaN.maxwvel=eps*1e6;
AdcUtil.MaxVarNaN.initiallydry=-1000;
AdcUtil.MaxVarNaN.maxinundepth=-1000;

