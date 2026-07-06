function [NArray, EArray, ZArray] = sweref(latArray,longArray,hArray)
%Converts coordinates in lat/long WGS 84 to SWEREF99TM
%
%[NArray, EArray, ZArray] = sweref(north,east,hight)
%
% accepts both decimal-degrees and "NMEA-format" (ddmm.mmmmmm)
if latArray(1)>90
    latArray=latArray/100;
    latArray=(floor(latArray)+(latArray-floor(latArray))*100/60) ;
end
latArray = latArray *pi/180;
if longArray(1)>90
    longArray=longArray/100;
    longArray=(floor(longArray)+(longArray-floor(longArray))*100/60) ;
end
longArray = longArray *pi/180;
%Calculations from www.lantmateriet.se.
NArray=zeros(length(latArray),1);  %preallocate to increase speed
EArray=zeros(length(latArray),1);
ZArray=zeros(length(latArray),1);
%Constants
a=6378137; %[m]
f=1/298.257222101;
e2=f*(2-f);
n=f/(2-f);
at=a/(1+n) * (1+n^2/4+n^4/64); %serie
A=e2;
B=1/6 * (5*e2^2-e2^3);
C=1/120 * (104*e2^3 -45*e2^4); %serie
D=1/1260 * (1237*e2^4);       %serie
long_av=pi/180*15.0;  %15� %Medelmeridianen i SWEREF99
k=0.9996;                  %skalfaktor p� medelmeridianen
FN=0;                      %x-till�gg (false northing)
FE=500000;                 %y-till�gg (false easting)
b1=1/2*n - 2/3*n^2 + 5/16*n^3 + 41/180*n^4;    %serie
b2=13/48*n^2 - 3/5*n^3 + 557/1440*n^4;         %serie
b3=61/240*n^3 - 103/140*n^4;                    %serie
b4=49561/161280*n^4;                            %serie
%Converstion from Lat,Long to X,Y in GRS80 (SWEREF).
%(The difference between WGS84 and GRS80 is neglected here.)
for i =1:length(latArray)
    lat=latArray(i); %Allways North! 
    long=longArray(i); %Allways East!
    h=hArray(i);
    dlong=long-long_av;
    %Calculating 'conformal lattitude'
    lat1=lat - sin(lat)*cos(lat)*(A + B*sin(lat)^2 + C*sin(lat)^4 + D*sin(lat)^6);
    es=atan(tan(lat1)/cos(dlong));
    ns=atanh(cos(lat1)*sin(dlong));
    NArray(i)=k*at*(es + b1*(sin(2*es)*cosh(2*ns)) + b2*(sin(4*es)*cosh(4*ns)) + b3*(sin(6*es)*cosh(6*ns)) + b4*(sin(8*es)*cosh(8*ns))) + FN;
    EArray(i)=k*at*(ns + b1*(cos(2*es)*sinh(2*ns)) + b2*(cos(4*es)*sinh(4*ns)) + b3*(cos(6*es)*sinh(6*ns)) + b4*(cos(8*es)*sinh(8*ns))) + FE;
    ZArray(i)=h;
end
%ZArray=ZArray';
%Calculating geoidheight for each point
%Loading geiod-grid for RH2000
load swen08
swen_H = interp2(swen08_X,swen08_Y,swen08_Z,latArray*180/pi,longArray*180/pi);
if size(ZArray,1)==size(swen_H,1)
    ZArray=ZArray-swen_H;
else
    ZArray=ZArray-swen_H';
end
