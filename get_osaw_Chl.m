function [osaw,a,bb,chl_derived]=get_osaw_Chl(wl,Rrs,solz,ws,opt)
%% Xiaolong Yu @ xmu.edu.cn
% calculate the water leaving albedo using Morel LUT (chl-based);
% revison history:
% draft on Jan 20,2021
% references: Xiaolong Yu et al., 2021, Comments on chlorophyll-a-based schemes for the estimation of water-leaving albedo

%% input
% wl: wavelength (nm)
% Rrs: remote senisng reflectance (sr-1)
% solz: solar zenith angle (deg)
% ws  : wind speed (m/s)

% build the parameterization control structure 'opt'
% opt.yrat: the faction of CDOM to particle absorption, optional, default 0.2
% opt.chla: optional chlorophyll concentration (mg/m3); if not present, use the OCI
%           algroithm (Hu et al., 2002) to derive chla
% opt.inverse: IOPs inverse algroithm, either 'QAA'  (Lee et al., 2002) or
%             chla-based, optional; 'GIOP' maybe included later
% opt.aw: pure seawater absortion, either 'Lee'  (Lee et al., 2015) or
%             default (Pope & Fry, 1997), optional 
% opt.bbw: pure seawater backscattering, either 'Zhang'  (Zhang et al., 2009) or
%             default (Morel, 1974), optional 
%% output
% Rw: water leaving albedo, defined as water leaving irradiance/downwelling
%     irradaince
% a:  estimated total absorption coefficients from chla , optional
% bb: estimated total backscattering coefficients from chla , optional 
% chl_oci: derived chla using OCI algorithm, optional 

%% main 
nwl=length(wl);
wl(1:nwl)=wl;
Rrs(1:nwl)=Rrs;  

%% select input chla
chl=chl_oci(wl,Rrs);
chl_derived=chl; 
chl(chl<0.02) = 0.02;  % force minimal chl = 0.02

%% select aw and bbw 
aw=opt.aw;
bbw=opt.bbw;

%% derive IOPs from Rrs 
yrat=opt.yrat;
[a, bb]= get_iops(wl,chl,opt);  % use chla estimated IOPs, either from input Chla or from OCI



%% angular quad--> match Hydrolight simulation
senz=[0,10,20,30,40,50,60,70,80,87.5];
phi=[0,15,30,45,60,75,90,105,120,135,150,165,180]; 

n_s=length(senz);
n_p=length(phi);
%% reflection–refraction term (see Morel & Gentili, 1996, AO)
RR=gothic_R(wl,solz,senz,ws);
%% get angular Rrs(senz,phi)
foq_data=opt.fq;

for i=1:n_s   
    for j=1:n_p 
      fq= morel_fq(chl,solz,senz(i),phi(j),wl,foq_data);   % get  f/Q for given geometry         
      Rrs_temp=RR(i,:).*fq.*(bb./a);       
      Rrs_simu(i,j,1:nwl)=Rrs_temp;
    end
end

for k = 1:nwl
    A_Rrs=Rrs_simu(:,:,k);    % angular Rrs
    osaw(k)=getEw_interp2(A_Rrs,senz,phi);
end


end
