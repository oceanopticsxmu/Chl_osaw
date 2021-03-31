function chl=chl_oc3c(wave, Rrs) 

a= [0.2515, -2.3798, 1.5823, -0.6372, -0.5692];
wl = [443, 490, 555];
[~,ib1] = min(abs(wave-wl(1)));
[~,ib2] = min(abs(wave-wl(2)));
[~,ib3] = min(abs(wave-wl(3)));
 
 Rrs1 = Rrs(ib1);
 Rrs2 = Rrs(ib2);
 Rrs3 = Rrs(ib3); 
 
% minRrs = MIN(Rrs1, Rrs2);  
% if (Rrs3 > 0.0 && Rrs2 > 0.0 && minRrs > -0.001) 
% Rrs3 = conv_rrs_to_555(wl(3),Rrs3);
% end

rat = max(Rrs1, Rrs2) / Rrs3;       
rat = log10(rat); 
chl=10.0^(a(1) + rat * (a(2) + rat * (a(3)+rat*(a(4) + rat * a(5)))));
       
end