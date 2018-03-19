function newzeta=AdcircExtendSurface(fgs,zeta)

ze=zeta(fgs.e);
i1=find(sum(isfinite(ze'))==1);                                                            
i2=find(sum(isfinite(ze'))==2);    

newzeta=zeta;

for i=1:length(i1)
    thise=fgs.e(i1(i),:);
    igood=find(isfinite(zeta(thise)));
    inan=find(isnan(zeta(thise)));
    newzeta(thise(inan))=newzeta(thise(igood));
end
for i=1:length(i2)
    thise=fgs.e(i2(i),:);
    igood=find(isfinite(zeta(thise)));
    inan=find(isnan(zeta(thise)));
    newzeta(thise(inan))=mean(newzeta(thise(igood)));
end
