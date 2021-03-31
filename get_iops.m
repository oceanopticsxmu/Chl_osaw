function [a,bb]= get_iops(wl,chl,opt)

%% Xiaolong Yu @ xmu.edu.cn
% calculate the water leaving albedo using Morel LUT (chl-based);
% revison history:
% draft on Jan 20,2021


%% input
% wl: wavelength (nm)
% yrat: the faction of CDOM to particle absorption, optional, default 0.2
% chla: optional chlorophyll concentration (mg/m3); if not present, use the OCI
%       algroithm (Hu et al., 2002) to derive chla
% opt.aw: pure seawater absortion, either 'Lee'  (Lee et al., 2015) or
%             default (Pope & Fry, 1997), optional 
% opt.bbw: pure seawater backscattering, either 'Zhang'  (Zhang et al., 2009) or
%             default (Morel, 1974), optional 
% opt.dat: Coefficients in Bricaud 1998, spectra with aph*(443) normalized to 0.055 m2/mg

%% output
% a:  estimated total absorption coefficients from chla , optional
% bb: estimated total backscattering coefficients from chla , optional 


n_wl=length(wl);
wl(1:n_wl)=wl;
[~,i440] =min(abs(wl-440));
%% get IOPs from Chla
% pure seawater IOPs
aw=opt.aw;
bbw=opt.bbw;
yrat=opt.yrat;

% total absorption coefficient
% [ap_temp,~,~,~] = get_bricaud_aph(chl,wl,'norm');
dat = opt.dat;
ap0   = dat(:,2).* chl.^dat(:,3);
ap_temp = interp1(dat(:,1),ap0,wl,'cubic');
nanidx=find(isnan(ap_temp));

if isnan(nanidx)
    ap_temp(nanidx)=ap_temp(nanidx(1)-1);
end

ap(1:n_wl)=ap_temp;

ay440=yrat*(aw(i440)+ap(i440));   % Eq. 18 in MOREL AND MARITORENA(2001);  constant fraction (0.2) is adopted from Eq. 18 of Prieur & S. 1981 
 
ay(1:n_wl)=ay440*exp(-0.014*(wl-440));    % Eq. 16 in MOREL AND MARITORENA(2001) 

a(1:n_wl)=ap+ay+aw;

% total backscattering coefficient
% bp550=0.3*chl^0.62-bbw(i550)*2;   %% Eq. 9 and 9' in MOREL AND MARITORENA(2001), source Gordon & Morel, 1983 
bp550=0.416*chl^0.766;   %% Eq. 12 in MOREL AND MARITORENA(2001), based on a larger dataset (Loisel & Morel, 1998)
if  chl >=0.02 && chl < 2 
    v=0.5*(log10(chl)-0.3);
else
    v=0;
end
bbp(1:n_wl)=(0.002+0.01*(0.5-0.25*log10(chl))*(wl/550).^v)*bp550; % Eq. 13 in MOREL AND MARITORENA(2001) 
bb(1:n_wl)=bbp+bbw;

end
