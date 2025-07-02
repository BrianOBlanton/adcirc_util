function AdcUtil=adcirc_util_init
echo off

AdcUtil.MaxVarMapping.maxele='zeta_max';
AdcUtil.MaxVarMapping.maxvel='vel_max';
AdcUtil.MaxVarMapping.maxwvel='wind_max';
AdcUtil.MaxVarMapping.initiallydry='initiallydry';
AdcUtil.MaxVarMapping.maxinundepth='inun_max';
AdcUtil.MaxVarMapping.swan_HS_max='swan_HS_max';
AdcUtil.MaxVarMapping.minpr='pressure_min';

AdcUtil.MaxVarNaN.maxele=-1000;
AdcUtil.MaxVarNaN.maxvel=eps*1e6;
AdcUtil.MaxVarNaN.maxwvel=eps*1e6;
AdcUtil.MaxVarNaN.initiallydry=-1000;
AdcUtil.MaxVarNaN.maxinundepth=-1000;
AdcUtil.MaxVarNaN.swan_HS_max=-1;
AdcUtil.MaxVarNaN.minpr=-1;

