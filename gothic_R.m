% /* ---------------------------------------------------------------------------- */
% /* gothic_R - returns total effect of transmission through air-sea interface    */
% Xiaolong Yu, xlyu@xmu.edu. Mar 31, 2021
%   ---------------------------------------------------------------------------- */
function R=gothic_R(wl,solz,senz,ws) 
  nw=1.34; 
  
  %% correction for solar-zenith angle effects, Wang et al., 2006
  brdf_sol = fresnel_sol(wl,solz,ws);
  %% use l2gen code
%   brdf_sen=brdf_sen';  
%   R = brdf_sol.*brdf_sen/nw/nw; 

%% use Gordon 2005 Fig. 5
  R_senz=get_RR(senz,ws);
  
  %% final Gothic_R
  R=R_senz'.*brdf_sol;
  
end

