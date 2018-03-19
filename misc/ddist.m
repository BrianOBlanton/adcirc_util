function d=ddist(Lat1,Lon1,Lat2,Lon2,R)

Lat1=Lat1*pi/180;
Lat2=Lat2*pi/180;
Lon1=Lon1*pi/180;
Lon2=Lon2*pi/180;

a=sin((Lat2-Lat1)/2).^2 + cos(Lat1) .* cos(Lat2) .* sin((Lon2-Lon1)/2).^2;
a(a<0)=0;
a(a>1)=1;
d = R * 2 * atan2(sqrt(a),sqrt(1 - a));

end


