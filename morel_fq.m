function fq = morel_fq(chl,solz,theta,relaz,wvl,foq_data)
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 201
% fqa = -999;
% fqc = -999;

nw = 1.34;
w = [412.5,442.5,490,510,560,620,660];
thetap = asin(sin(theta *(pi/180))/nw)/(pi/180);

fq = [];
f0 = [];

for i = 1:7
	fq = [fq,get_fq(w(i),solz,chl,thetap,relaz,foq_data)];
% 	f0 = [f0,get_fq(w(i),0,chl,0,0,foq_data)];
end
% 
% fc = f0 ./ fq;

if exist('wvl') == 0; wvl = 380:1:700; end

fq = interp1(w,fq,wvl,'pchip');
% fc = interp1(w,fc,wvl,'pchip');
end