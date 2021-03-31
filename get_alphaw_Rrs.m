function [osaw,chl_derived]=get_alphaw_Rrs(wl,Rrs,solz,ws,foq_data)
%% Xiaolong Yu @ xmu.edu.cn
% calculate the water leaving albedo using Morel LUT (chl-based);
% revison history:
% draft on Jan 20,2021
% 2020-3-23; skip Chl-inferred IOPs, use LUT to get angular Rrs
% references: Xiaolong Yu et al., 2021, Comments on chlorophyll-a-based schemes for the estimation of water-leaving albedo

%% input
% wl: wavelength (nm)
% Rrs: remote senisng reflectance (sr-1)
% solz: solar zenith angle (deg)
% ws  : wind speed (m/s)

% build the parameterization control structure 'opt'
% opt.chla: optional chlorophyll concentration (mg/m3); if not present, use the OCI
%           algroithm (Hu et al., 2002) to derive chla
% foq data: LUT for f/Q 
%% output
% osaw: water leaving albedo, defined as water leaving irradiance/downwelling
%     irradaince
% chl_oci: derived chla using OCI algorithm, optional 

%% main  

nwl=length(wl);
wl(1:nwl)=wl;
Rrs(1:nwl)=Rrs; 

%% select input chla
chl=chl_oci(wl,Rrs);
chl_derived=chl; 
chl(chl<0.02) = 0.02;  % force minimal chl = 0.02


%% angular quad--> match Hydrolight simulation
senz=[0,10,20,30,40,50,60,70,80,87.5];
phi=[0,15,30,45,60,75,90,105,120,135,150,165,180]; 

n_s=length(senz);
n_p=length(phi);
%% reflection–refraction term (see Morel & Gentili, 1996, AO)
RR=gothic_R(wl,solz,senz,ws);

%% get angular Rrs(senz,phi)
% fq= get_fq_VIIRS(solz,chl,foq);   % speed up computation time for VIIRS imagery
% f0(1:nwl)=fq(1,1,:);

f0= morel_fq(chl,solz,0,0,wl,foq_data);

for i=1:n_s  
    R_ratio=RR(i,:)/RR(1,:);     
    for j=1:n_p 
%         fq_temp(1:nwl)=fq(i,j,:);
      fq_temp = morel_fq(chl,solz,senz(i),phi(j),wl,foq_data);   % get  f/Q for given geometry   
      Rrs_temp=R_ratio*fq_temp./f0.*Rrs;       
      Rrs_simu(i,j,1:nwl)=Rrs_temp;
    end
end

for k = 1:nwl
    A_Rrs=Rrs_simu(:,:,k);  % angular Rrs
    osaw(k)=getEw_interp2(A_Rrs,senz,phi);
end

end
