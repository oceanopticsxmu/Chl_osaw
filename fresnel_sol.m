%/* ----------------------------------------------------------------------------  */
%/* fresnel_sol() - effects of the air-sea transmittance for solar path           */
%/* source code: https://oceancolor.gsfc.nasa.gov/docs/ocssw/brdf_8c_source.html#l00314  */                                                                          */
%/* Description:                                                                  */
%/*   This computes the correction factor on normalized water-leaving radiance    */
%/*   to account for the solar zenith angle effects on the downward irradiance    */
%/*   from above the ocean surface to below the surface.                          */
%/*   Menghua Wang 9/27/04.                                                       */                                                                             */
%/*   Added windspeed dependence, December 2004, BAF                              */
%/*   Modified to return air-sea transmittance as option, December 2008, BAF      */  

% Inputs:
% c!        wl  --- wavelength.
% c!        solz  --- solar zenith angle.
% c!        mu0  --- cosine of solar zenith angle.optional
% c!        wspeed  --- wind speed in (m/s).
% c!
% c! Outputs:
% c!        corr_f(n)  --- correction factors for the SeaWiFS and MODIS
% c!                ocean bands 1-6.  The correction factor for the MODIS 
% c!                bands 667-nm and 678-nm above are the same, i.e., using 
% c!                the same correction factor for the MODIS spectral bands 
% c!                for wavelength > 667 nm. 
% c!        The correction factor is applied as:
% c!                nLw(new) = corr_f * nLw(old) 
% c!
% c! Reference:
% c!        Wang, M.,(2006) "Effects of ocean surface reflectance variation with
% c!        solar elevation on normalized water-leaving radiance", Appl. Opt.
 
%/*   Converted to matlab code by Xiaolong Yu(xlyu@xmu.edu.cn),Mar 21,2021        */
%----------------------------------------------------------------------------     */ 

function brdf_sol =fresnel_sol(wl,solz,ws)
    
  radeg=pi/180;     
  twave= [412, 443, 490, 510, 555, 670];  
  
    
% /* M Wang, personal communication, red-nir iterpolated */
%   tf0_w =[412, 443, 490, 510, 555, 670, 765, 865];
%   tf0_v= [0.965980, 0.968320, 0.971040, 0.971860, 0.973450, 0.977513, 0.980870, 0.984403];
%   tf0 =interp1(tf0_w,tf0_v,wl);
 tf0= [0.965980, 0.968320, 0.971040, 0.971860, 0.973450, 0.977513]; % match twave
  
% Coefficients in Eq.(16) of Wang 2006. c=size(coeff.,ws,wave) %  ws=[0,1.9,7.5,16.9,30];   
% { /* ws=0.0 */
c(:,1,:)=[-0.0087, -0.0122, -0.0156, -0.0163, -0.0172, -0.0172;
   0.0638, 0.0415, 0.0188, 0.0133, 0.0048, -0.0003;
   -0.0379, -0.0780, -0.1156, -0.1244, -0.1368, -0.1430;
   -0.0311, -0.0427, -0.0511, -0.0523, -0.0526, -0.0478];
% { /* ws=1.9 */
c(:,2,:)=[-0.0011, -0.0037, -0.0068, -0.0077, -0.0090, -0.0106;
   0.0926, 0.0746, 0.0534, 0.0473, 0.0368, 0.0237;
   -5.3E-4, -0.0371, -0.0762, -0.0869, -0.1048, -0.1260;
   -0.0205, -0.0325, -0.0438, -0.0465, -0.0506, -0.0541];
% { /* ws=7.5 */
c(:,3,:)=[6.8E-5, -0.0018, -0.0011, -0.0012, -0.0015, -0.0013;
   0.1150, 0.1115, 0.1075, 0.1064, 0.1044, 0.1029;
   0.0649, 0.0379, 0.0342, 0.0301, 0.0232, 0.0158;
   0.0065, -0.0039, -0.0036, -0.0047, -0.0062, -0.0072];
% { /* ws=16.9 */
c(:,4,:)=[-0.0088, -0.0097, -0.0104, -0.0106, -0.0110, -0.0111;
  0.0697, 0.0678, 0.0657, 0.0651, 0.0640, 0.0637;
  0.0424, 0.0328, 0.0233, 0.0208, 0.0166, 0.0125;
  0.0047, 0.0013, -0.0016, -0.0022, -0.0031, -0.0036];  
% { /* ws=30 */
c(:,5,:)=[-0.0081, -0.0089, -0.0096, -0.0098, -0.0101, -0.0104;
 0.0482, 0.0466, 0.0450, 0.0444, 0.0439, 0.0434;
 0.0290, 0.0220, 0.0150, 0.0131, 0.0103, 0.0070;
 0.0029, 0.0004, -0.0017, -0.0022, -0.0029, -0.0033]; 

tsigma= [0.0, 0.1, 0.2, 0.3, 0.4];    %   tsigma = 0.0731 * sqrt(ws0);
sigma = 0.0731 * sqrt(ws);

x = log(cos(min(solz, 80.0)*radeg)); 
x2 = x*x;
x3 = x*x2;
x4 = x*x3;

%    /* find bracketing table winds */
%  Figure out the wind speed index to locate proper pair of sigmas
sidx=find(tsigma-sigma>0);
is2=sidx(1); 
is1=is2-1; 
slp=(sigma-tsigma(is1))/(tsigma(is2)-tsigma(is1));   
% * compute at bounding windspeeds and interpolate */ 
tf1=1.0+c(is1,1,:)*x+c(is1,2,:)*x2+c(is1,3,:)*x3+c(is1,4,:)*x4;
tf2=1.0+c(is2,1,:)*x+c(is1,2,:)*x2+c(is2,3,:)*x3+c(is2,4,:)*x4;                 
tf= tf1+slp*(tf2 - tf1); 
tf=squeeze(tf);

brdf_temp= tf0'./tf; 
brdf_sol =interp1(twave,brdf_temp,wl);
idx1=find(wl<412);
idx2=find(wl>670);

brdf_sol(idx1)=brdf_temp(1);
brdf_sol(idx2)=brdf_temp(end);


end


 