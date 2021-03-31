function chl=chl_oci(wave,Rrs) 
 t1 = 0.15;
 t2 = 0.20;

chl1 = chl_hu(wave, Rrs);

if   chl1 <= t1      
     chl = chl1;
else 
     chl2 = chl_oc3c(wave, Rrs); 
end

if chl1 >= t2    
    chl = chl2;
end

if chl1>t1 && chl1<t2
     chl = chl1 * (t2 - chl1) / (t2 - t1)+ chl2 * (chl1 - t1) / (t2 - t1);
end

end

  