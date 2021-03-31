clear;
clc;
tic;

load('HLosa_morel_final.mat'); 
wl=HLosa.wl;
wlno=length(wl);
Rrs=HLosa.nRrs;
solz=HLosa.sza;
ws=HLosa.ws;
Sno=length(solz);
IOCCG_chl=[0.03,0.05,0.07,0.1,0.15,0.2,0.3,0.5,0.7,1,1.5,2,3,5,7,10,15,20,25,30];  % every 25 samples

%% read f/Q data
foq_data = read_fq;

for fcnt=1:Sno
    isolz=solz(fcnt);
    iws=ws(fcnt); 
    iRrs(1:wlno)=Rrs(:,fcnt);
    [osaw_temp,chl_temp]=get_alphaw_ratio(wl,iRrs,isolz,iws,foq_data); 
    chl_oci(fcnt)=chl_temp;
    osaw(:,fcnt)=osaw_temp;    
    disp(['processing the No. ' num2str(fcnt) ' file;']);
   
end

HLosa.osaw4=osaw;
HLosa.chl_oci=chl_oci;
HLosa.info={'osaw:using fq ratio';'osaw1:using fq ratio, solz corrected';'osaw2: using CHl, solz corrected';'osaw3: using CHl, solz corrected,Gordon05-RR';'osaw3: using Ratio, solz corrected,Gordon05-RR'};
save('HLosa_morel_final.mat','HLosa'); 

toc;
