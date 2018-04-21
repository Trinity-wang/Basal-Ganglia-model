function [wsg,wgs,wlatstn,wlatgpe]= weightcal(D2)
% -------------------parameters for version 5 -------------%
ssmax = 0.15;  %Strength of lateral connections in stn
sgmax = 1;   %Strength of lateral connections in gpe
wsg = 1; %Strength of  connections from stn to gpe
wgs =20; %Strength of  connections from gpe to stn
rs = 5; %Radius of lateral connections in stn
rg = 25; %Radius of lateral connections in gpe
nlatstn = 11; % number of laterals in stn
nlatgpe = 15;% number of laterals in gpe

CD21=1;
CD2=0.1;

rs=rs/(CD21*D2);
rg=rg/(1-D2);
wsg=((1-CD2*D2))*wsg;
wgs=((1-CD2*D2))*wgs;
wlatstn = calclatwts(nlatstn,ssmax, rs);
wlatgpe = calclatwts(nlatgpe,sgmax, rg);

end