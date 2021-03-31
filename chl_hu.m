function chl= chl_hu(wave, Rrs)  

 wl = [443, 555, 670];
 a = [-0.4909, 191.6590];  % Eq. (4) in Hu et al., 2012
 [~,ib1] = min(abs(wave-wl(1)));
 [~,ib2] = min(abs(wave-wl(2)));
 [~,ib3] = min(abs(wave-wl(3)));
 
 Rrs1 = Rrs(ib1);
 Rrs2 = Rrs(ib2);
 Rrs3 = Rrs(ib3); 
% Eq. (3) in Hu et al., 2012
 ci = min(Rrs2 - (Rrs1 + (wl(2) - wl(1)) / (wl(3) - wl(1))*(Rrs3 - Rrs1)), 0.0);
 chl = 10.0^(a(1) + a(2)* ci);

end

 
 