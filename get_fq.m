function fqint = get_fq(wl,solz,chl_in,senzp,phi,foq)
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013
% translated from the C code provided in SeaDAS/l2gen

fqint = -999;

wvl = [412.5,442.5,490,510,560,620,660];
sun = [0,15,30,45,60,75];
chl = [0.03,0.1,0.3,1,3,10];
nad = [1.078,3.411,6.289,9.278,12.3,15.33,18.37,21.41,24.45,27.5,30.54,33.59,36.64,39.69,42.73,45.78,48.83];	
azm = [0,15,30,45,60,75,90,105,120,135,150,165,180];

nw = 7;
ns = 6;
nc = 6;
nn = 17;
na = 13;

lchl = log(chl);
c = log(chl_in);

% locate nearest wavelength
%
iw = find(abs(wvl - wl) < 15);
if isempty(iw) 
	'unrecognized wavelength'
	return
else
	iw = iw(1);
end

% locate boundary zenith

if solz <= sun(1)
    solz = sun(1);
	js = 1;
elseif solz >= sun(ns)
    solz=sun(ns);
	js = ns-1;
else
	js = 1;
	while solz > sun(js + 1); js = js + 1; end
end


% locate boundary chl
%
if c <= lchl(1)
    c = lchl(1);
	kc = 1;
elseif c >= lchl(nc)
    c=lchl(nc);
	kc = nc-1;
else
	kc = 1;
	while c > lchl(kc + 1); kc = kc + 1; end
end

% locate boundary nadir 
if senzp <= nad(1) 
	senzp = nad(1); % force min nadir angle to min table value
	ln = 1;
elseif senzp >= nad(nn)
	ln = nn-1;
else
	ln = 1;
	while senzp > nad(ln + 1); ln = ln + 1; end
end

% locate boundary azimuth

if phi <= azm(1)
    phi = azm(1);
	ma = 1;
elseif phi >= azm(na)
    phi = azm(na);
	ma = na-1;
else
	ma = 1;
	while phi > azm(ma + 1); ma = ma + 1; end
end


% weight, interpolate, and iterate values to solve for f/Q

ds(1) = (sun(js+1) - solz) / (sun(js+1) - sun(js));
ds(2) = (solz - sun(js)) / (sun(js+1) - sun(js));

dc(1) = (lchl(kc+1) - c) / (lchl(kc+1) - lchl(kc));
dc(2) = (c - lchl(kc)) / (lchl(kc+1) - lchl(kc));

dn(1) = (nad(ln+1) - senzp) / (nad(ln+1) - nad(ln));
dn(2) = (senzp - nad(ln)) / (nad(ln+1) - nad(ln));

da(1) = (azm(ma+1) - phi) / (azm(ma+1) - azm(ma));
da(2) = (phi - azm(ma)) / (azm(ma+1) - azm(ma));



fqint = 0;

for j = 1:2
for k = 1:2
for l = 1:2
for m = 1:2
	fqint = fqint + ds(j) .* dc(k) .* dn(l) .* da(m) .* foq(iw, js+j-1, kc+k-1, ln+l-1, ma+m-1);
end
end
end
end
