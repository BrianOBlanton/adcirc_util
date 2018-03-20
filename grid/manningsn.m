function cd=manningsn(z,n,lowerboundonCd)
%  cd=manningsn(z,n,lowerboundonCd)

temp=z;
temp(abs(z)<.1)=.1;
cd=(9.81*n.^2)./(abs(temp).^(1/3));

if exist('lowerboundonCd')
    cd(cd<lowerboundonCd)=lowerboundonCd;
end


%{
      if (nlcd_class.eq.11) mannings_n=0.02    !Open Water
      if (nlcd_class.eq.12) mannings_n=0.010   !Perennial Ice/Snow
      if (nlcd_class.eq.21) mannings_n=0.020   !Developed - Open Space
      if (nlcd_class.eq.22) mannings_n=0.050   !Developed - Low Intensity
      if (nlcd_class.eq.23) mannings_n=0.100   !Developed - Medium Intensity
      if (nlcd_class.eq.24) mannings_n=0.150   !Developed - High Intensity
      if (nlcd_class.eq.31) mannings_n=0.090   !Barren Land (Rock/Sand/Clay)
      if (nlcd_class.eq.32) mannings_n=0.040   !Unconsolidated Shore
      if (nlcd_class.eq.41) mannings_n=0.100   !Deciduous Forest
      if (nlcd_class.eq.42) mannings_n=0.110   !Evergreen Forest
      if (nlcd_class.eq.43) mannings_n=0.100   !Mixed Forest
      if (nlcd_class.eq.51) mannings_n=0.040   !Dwarf Scrub
      if (nlcd_class.eq.52) mannings_n=0.050   !Shrub/Scrub
      if (nlcd_class.eq.71) mannings_n=0.034   !Grassland/Herbaceous
      if (nlcd_class.eq.72) mannings_n=0.030   !Sedge/Herbaceous
      if (nlcd_class.eq.73) mannings_n=0.027   !Lichens
      if (nlcd_class.eq.74) mannings_n=0.025   !Moss
      if (nlcd_class.eq.81) mannings_n=0.033   !Pasture/Hay
      if (nlcd_class.eq.82) mannings_n=0.037   !Cultivated Crops
      if (nlcd_class.eq.90) mannings_n=0.100   !Woody Wetlands
      if (nlcd_class.eq.91) mannings_n=0.100   !Palustrine Forested Wetland
      if (nlcd_class.eq.92) mannings_n=0.048   !Palustrine Scrub/Shrib Wetland
      if (nlcd_class.eq.93) mannings_n=0.100   !Estuarine Forested Wetland
      if (nlcd_class.eq.94) mannings_n=0.048   !Estuarine Scrub/Shrub Wetland
      if (nlcd_class.eq.95) mannings_n=0.045   !Emergent Herbaceous Wetlands
      if (nlcd_class.eq.96) mannings_n=0.045   !Palustrine Emergent Wetland (Persistant)
      if (nlcd_class.eq.97) mannings_n=0.045   !Estuarine Emergent Wetland
      if (nlcd_class.eq.98) mannings_n=0.015   !Palustrine Aquatic Bed
      if (nlcd_class.eq.99) mannings_n=0.015   !Estuarine Aquatic Bed
%}