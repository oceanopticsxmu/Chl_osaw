%% 
load('HLosa_morel_final.mat'); 
wl=HLosa.wl;
id700=find(wl==700);
wl=wl(1:id700);
wlno=length(wl);
Rrs=HLosa.nRrs;
Rrs=Rrs(1:id700,:);

solz=HLosa.sza;
ws=HLosa.ws;
Sno=length(solz);
IOCCG_chl=[0.03,0.05,0.07,0.1,0.15,0.2,0.3,0.5,0.7,1,1.5,2,3,5,7,10,15,20,25,30];  % every 25 samples

%% read f/Q data
foq_data = read_fq;

%% Bricaud 1998 spectra with aph*(443) normalized to 0.055 m2/mg
% [ap_temp,~,~,~] = get_bricaud_aph(chl,wl,'norm');
aphDat = load('bricaud_1998_aph.txt');
opt.dat=aphDat;
%% structure for the input, presented in opt
opt.fq=foq_data;
opt.yrat=0.2;       % the ratio of yellow substance to particle absorption, default 0.2 (Morel et al, 2001)
opt.bbw=0.5*h2o_iops_Zhh_lee(wl,'b');  % if not presented, default Morel, 1974
opt.aw= h2o_iops_Zhh_lee(wl,'a');     % if not presented, default Pope and Fry, 1997

for fcnt=1:Sno
    isolz=solz(fcnt);
    iws=ws(fcnt); 
    iRrs(1:wlno)=Rrs(:,fcnt);
    [osaw_temp,a_temp,bb_temp,chl_temp]=get_osaw_Chl(wl,iRrs,isolz,iws,opt); 

    chl_oci(fcnt)=chl_temp;
    osaw(:,fcnt)=osaw_temp;  
    a(:,fcnt)=a_temp;    
    bb(:,fcnt)=bb_temp;    
    disp(['processing the No. ' num2str(fcnt) ' file;']);
   
end

HLosa.osaw2=osaw;
HLosa.chl_a=a;
HLosa.chl_bb=bb;
HLosa.chl_oci2=chl_oci;
HLosa.info={'osaw1:using fq ratio, solz corrected';'osaw2: using CHl, solz corrected,Gordon05-RR'};
save('HLosa_morel_final.mat','HLosa'); 

toc;
