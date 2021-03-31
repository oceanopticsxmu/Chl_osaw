function [ap,aps,aph,aphs] = get_bricaud_aph(chl,wl,norm)

%
% Bricaud et al., JGR 103, 31033-31044, 1998
%
% ap = Ap Chl^Ep
% ap* = Ap Chl^(Ep-1)
% aph = Aph Chl^Eph
% aph* = Aph Chl^(Eph-1)
% bricaud_1998_aph.txt order of columns: wl,Ap,Ep,Aph,Eph
% empirical coefficients are slight different from Table 2 of Bricaud et al., 1995, as enlarged
% dataset is used here.
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013
%

dat = load('bricaud_1998_aph.txt');

ap0   = dat(:,2) .* chl.^dat(:,3);
% aps0  = dat(:,2) .* chl.^(dat(:,3) - 1);
% aph0  = dat(:,4) .* chl.^dat(:,5);
% aphs0 = dat(:,4) .* chl.^(dat(:,5) - 1);

% 
% if exist('norm') == 1
% 	idx = find(dat(:,1) == 442);
% 	aphs0 = aphs0 * 0.055 / aphs0(idx);
% end


ap   = interp1(dat(:,1),ap0,wl,'cubic');
% aps  = interp1(dat(:,1), aps0,wl,'cubic');
% aph  = interp1(dat(:,1), aph0,wl,'cubic');
% aphs = interp1(dat(:,1),aphs0,wl,'cubic');


