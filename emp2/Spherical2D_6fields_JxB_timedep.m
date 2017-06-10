%% 2D spherical coordinates!

clear all; close all;

tic;

dopml = 1;
doplot = 1;
doJ = 0;
doioniz = 0;

fprintf('Initializing parameters and fields...\n');

e0 = 8.845e-12;
u0 = 4*pi*1e-7;
sig = 0; %1e-8;
sigm = 0;               % in this version we have updated the equations to include sigma_m
vp = 1/sqrt(e0*u0);

f0 = 30e3;              % 10 kHz
lam0 = vp/f0;           % about 30 km

r0 = 6370e3;            % a 1000 km radius "earth"
r1 = r0 + 110e3;        % a 200 km altitude "ionosphere"

dr1 = 1000;
dr2 = 500;
stepalt = 70e3;
stepi = stepalt/dr1;
r = [r0:dr1:(r0+stepalt) (r0+stepalt+dr2):dr2:r1];
rr = length(r);
dr = diff(r);

range = 100e3;          % 300 km along the ground
thmax = range / r0;     % gives max theta in radians
dth = dr1 / r0;          % steps in theta
th = 0:dth:thmax;
hh = length(th);

%dt = 1/vp/sqrt(1/dr1^2 + 1/(r0*dth/2)^2) / stabfac;      % r0*dth gives minimum dth in km
dt = 1e-7;

maxdist = max([(r1-r0) range]);

tsteps = round( 1.5 * maxdist/dt/vp );

c1e = (2*e0 - sig*dt)/(2*e0 + sig*dt);
c2e = 2*dt/(2*e0 + sig*dt);
c1h = (2*u0 - sigm*dt)/(2*u0 + sigm*dt);
c2h = 2*dt/(2*u0 + sigm*dt);

%% ----------------
%  initialize fields

Er = zeros(rr-1,hh);
Et = zeros(rr,hh-1);
Ep = zeros(rr,hh);
Hr = zeros(rr,hh-1);
Ht = zeros(rr-1,hh);
Hp = zeros(rr-1,hh-1);
Ep_old = Ep;
Er_old = Er;
Et_old = Et;

% time dependent quantities depend on E magnitude, which we will need.

Emag = zeros(rr,hh);
Eeff = zeros(rr,hh);
Epar = zeros(rr,hh);
Eperp = zeros(rr,hh);

% for Lee and Kalluri method, J located at integer grid locations.
% Collision frequency, Ne, B0, etc also at those locations.

Jr = zeros(rr,hh);
Jt = zeros(rr,hh);
Jp = zeros(rr,hh);

Jr0 = Jr;
Jt0 = Jt;
Jp0 = Jp;


%% initialize source

fprintf('Setting up current source...\n');

Jin = zeros(tsteps,1);

I0 = 200e3;

taur = 30e-6;
tauf = 70e-6;

for t = 1:tsteps,
    if (t*dt < taur),
        Jin(t) = I0 * t*dt/taur;
    else
        Jin(t) = I0 * exp(-((t*dt-taur)/tauf)^2);
    end
    
    % sinusoidal source
    Jin(t) = I0 * sin(2*pi*10e3*t*dt);

end


% bad ringing: let's filter it hard at 50 kHz. sample rate is 1/dt ~ 4 MHz
fcut = 25e3;
[b,a] = butter(2,fcut*dt/2);

Jin2 = filter(b,a,Jin);

Jin2 = I0 * Jin2 / max(Jin2);

sa = round(6e3/dr1);     % source altitude in dr steps, given 6 km
Js = zeros(2*sa,sa);


%% -------------------
% electron density and collisions (very preliminary!)

fprintf('Setting up Ionosphere parameters...\n');

beta = 0.5;         % per km
hp = 84.2;          % km
q = 1.6e-19;
me = 9.1e-31;

% Electron density from IRI

ne2 = zeros(size(r)); %IRIionosphere1((r-r0)/1000);      % per m^3

% below 60 km, use exponential
fi = find(ne2 > 0,1,'first');
logne = interp1([1 fi],[-5 log10(ne2(fi))],1:fi,'linear');
ne2(1:fi) = 10.^(logne);

ne2d = repmat(ne2',1,hh);
wp = sqrt(q^2 * ne2d / me / e0);

% get collision frequency from neutral density

ndcm = zeros(size(r)); %MSISatmosphere1((r-r0)/1000);      % per cm^3
nd = ndcm; %.total * 1e6;                   % per m^3
nd0 = nd(1);

% breakdown field
Ek = 3.2e6 * nd/nd0;
zeros(rr,hh);
% when electric field is zero,

mue = repmat((1.36 * nd0./nd)', 1, hh);
nu = (q/me) ./ mue;

%Bmag = 0;
Bmag = 50000e-9;        % 50,000 nT
Bvec = [1 0 0];         % B is in r-direction

B0 = Bmag*Bvec;
wc = -q*B0/me;

wcr = wc(1);
wct = wc(2);
wcp = wc(3);
wc0 = sqrt(wcr^2 + wct^2 + wcp^2);


%% ------------
% initialize ionization and attachment coefficients

vi = zeros(rr,hh);
va = zeros(rr,hh);

% optics

nN22P = zeros(rr,hh);
nN21P = zeros(rr,hh);

%% ------------
% set up optical integration to camera view

fprintf('Setting up elve parameters...\n');

camera.dist = 500e3;        % 500 km away
camera.th = camera.dist/r0;
camera.alt = 0;             % on the ground
camera.fovlr = 40*pi/180;
camera.fovud = 20*pi/180;
camera.upangle = 10*pi/180;
Nphi = 360;
elvemod = 20;       % fraction of total time steps to use in elve.

elvesteps = round(tsteps/elvemod);
et0 = (camera.dist - range) / vp;
etf = (camera.dist + range) / vp + dt*tsteps;
elvetimes = linspace(et0,etf,elvesteps);
elvedt = elvetimes(2) - elvetimes(1);

imsize = 128;

elve = zeros(imsize/2,imsize,elvesteps);

% elveparams are az, el, delay, and length, for each cell r, th, phi, where
% phi has been wrapped around 360 degrees.

%elveparams = opticalPropagation(r,th,Nphi,camera,elvetimes,elve,dr1,dr2);

n21pcube = zeros(rr,hh,floor(tsteps/elvemod));

%% ------------
% coordinate transformation

x = zeros(rr,hh);
y = zeros(rr,hh);
for m = 1:rr,
    for n = 1:hh,
        x(m,n) = r(m) * sin(th(n))/1000;
        y(m,n) = r(m) * cos(th(n))/1000;
    end
end

%% CPML stuff

fprintf('Setting up PML...\n');

if dopml,
    pmlm = 4;
    pmllen = 10;
    [A,B,P] = setupPML2d6f(pmlm,pmllen,dr1,dth,dt,rr,hh);
end

%% ----------------
%  setup figure

if doplot,
    
    c = colormap('jet');
    c3 = [0 0 0; 0 0 0.1; 0 0 0.2; 0 0 0.3; 0 0 0.4; c];
    c2 = [1 1 1; 0.8 0.8 0.9; 0.6 0.6 0.8; 0.4 0.4 0.7; 0.2 0.2 0.6; c];
    
    h1 = figure(1); 
    set(h1,'units','normalized','position',[0.1 0.2 0.8 0.7]);
    set(h1,'Renderer','painters');
    ax1 = subplot(231);
    ax2 = subplot(232);
    ax3 = subplot(233);
    ax4 = subplot(234);
    ax5 = subplot(235);
    ax6 = subplot(236);
    
    im1 = pcolor(x(1:end-1,:),y(1:end-1,:)-r0/1000,log10(abs(Er)+1e-6),'parent',ax1);
    im2 = pcolor(x,y-r0/1000,log10(abs(Jr)+1e-6),'parent',ax2);
    im3 = pcolor(x,y-r0/1000,log10(abs(Ep)+1e-8),'parent',ax3);
    im4 = pcolor(x,y-r0/1000,log10(abs(Jr)+1e-8),'parent',ax4);
    im5 = pcolor(x,y-r0/1000,log10(abs(Jt)+1e-8),'parent',ax5);
    im6 = pcolor(x,y-r0/1000,log10(abs(Jp)+1e-8),'parent',ax6);
    
    shading(ax1,'flat');
    shading(ax2,'flat');
    shading(ax3,'flat');
    shading(ax4,'flat');
    shading(ax5,'flat');
    shading(ax6,'flat');
    colormap(ax1,c2);
    colormap(ax2,c2);
    colormap(ax3,c2);
    colormap(ax4,c2);
    colormap(ax5,c2);
    colormap(ax6,c2);
    caxis(ax1,[-4 6]);
    caxis(ax2,[-8 -4]);
    caxis(ax3,[-4 1]);
    caxis(ax4,[4 11]);
    caxis(ax5,[-10 10]);
    caxis(ax6,[0 6]);
    xlabel(ax4,'Range along ground (km)');
    ylabel(ax1,'Altitude (km)');
    ylabel(ax4,'Altitude (km)');
    ti1 = title(ax1,'Er at time 0 us');
    ti2 = title(ax2,'Jr at time 0 us');
    ti3 = title(ax3,'Eeff/Ek at time 0 us');
    ti4 = title(ax4,'nu at time 0 us');
    ti5 = title(ax5,'Ne at time 0 us');
    ti6 = title(ax6,'N21P at time 0 us');
    
end;

if dopml,
    rshift = rr-pmllen-1;
    hshift = hh-pmllen-1;
end

%% ----------------
% update equations

sferic = zeros(1,tsteps);

fprintf('Beginning calculations...\n');

ii = 1;

for t = 1:tsteps,
    
    % Psi updates for H field --------------------
    
    if dopml,
        for i = 1:pmllen,
            for j = 1:hh-1,
                P.htr(i,j) = B.htr(i) * P.htr(i,j) + A.htr(i) * ( r(i+rshift+1) .* Ep(i+rshift+1,j) - r(i+rshift) * Ep(i+rshift,j) ) / dr(i+rshift);
            end
        end
        for i = 1:rr-1,
            for j = 1:pmllen,
                P.hrt(i,j) = B.hrt(j) * P.hrt(i,j) + ...
                    A.hrt(j) * ( sin(th(j+hshift+1)) * Ep(i,j+hshift+1) - sin(th(j+hshift)) * Ep(i,j+hshift) ) / dth;
            end
        end
        
        for i = 1:pmllen,
            for j = 1:hh-1,
                P.hpr(i,j) = B.hpr(i) * P.hpr(i,j) + A.hpr(i) * ( r(i+rshift+1) * Et(i+rshift+1,j) - r(i+rshift) * Et(i+rshift,j) ) / dr(i+rshift);
            end
        end
        for i = 1:rr-1,
            for j = 1:pmllen,
                P.hpt(i,j) = B.hpt(j) * P.hpt(i,j) + A.hpt(j) * (Er(i,j+hshift+1) - Er(i,j+hshift)) / dth;
            end
        end
    end
    
    % Hr update -------------------------------------
    for i = 1:rr,
        for j = 1:hh-1,
            Hr(i,j) = c1h * Hr(i,j) - c2h / r(i) / sin(th(j)+dth/2) * ...
                ( sin(th(j+1)) * Ep(i,j+1) - sin(th(j)) * Ep(i,j) ) / dth;
        end
    end
    
    if dopml,
        % PML correction
        for i = 1:rr-1,
            for j = 1:pmllen,
                Hr(i,j+hshift) = Hr(i,j+hshift) - c2h / r(i) / sin(th(j+hshift)+dth/2) * P.hrt(i,j);
            end
        end
    end
    
    % Ht update --------------------------------------
    for i = 1:rr-1,
        for j = 1:hh,
            Ht(i,j) = c1h * Ht(i,j) + c2h / (r(i)+dr(i)/2) * ...
                ( r(i+1) .* Ep(i+1,j) - r(i) * Ep(i,j) ) / dr(i);
        end
    end
    % PML correction
    if dopml,
        for i = 1:pmllen,
            for j = 1:hh-1,
                Ht(i+rshift,j) = Ht(i+rshift,j) + c2h / (r(i+rshift)+dr(i+rshift)/2) * P.htr(i,j);
            end
        end
    end
    
    % Hphi update ------------------------------
    for i = 1:rr-1,
        for j = 1:hh-1,
            Hp(i,j) = c1h * Hp(i,j) - c2h / (r(i)+dr(i)/2) * ...
                ( ( r(i+1) * Et(i+1,j) - r(i) * Et(i,j) ) / dr(i) - (Er(i,j+1) - Er(i,j)) / dth );
        end
    end
    % PML correction
    if dopml,
        for i = 1:pmllen,
            for j = 1:hh-1,
                Hp(i+rshift,j) = Hp(i+rshift,j) - c2h / (r(i+rshift)+dr(i+rshift)/2) * P.hpr(i,j);
            end
        end
        for i = 1:rr-1,
            for j = 1:pmllen,
                Hp(i,j+hshift) = Hp(i,j+hshift) + c2h / (r(i)+dr(i)/2) * P.hpt(i,j);
            end
        end
    end
    
    
    % Psi updates for E field --------------------
    
    if dopml,
        
        for i = 1:pmllen,
            for j = 1:hh-1,
                P.etr(i,j) = B.etr(i) * P.etr(i,j) + ...
                    A.etr(i) * ( (r(i+rshift)+dr(i+rshift)/2) * Hp(i+rshift,j) - (r(i+rshift)-dr(i+rshift-1)/2) * Hp(i+rshift-1,j) ) / ((dr(i+rshift)+dr(i+rshift-1))/2);
            end
        end
        for i = 1:rr-1,
            for j = 1:pmllen,
                P.ert(i,j) = B.ert(j) * P.ert(i,j) + ...
                    A.ert(j) * ( sin(th(j+hshift)+dth/2) * Hp(i,j+hshift) - sin(th(j+hshift)-dth/2) * Hp(i,j+hshift-1)) / dth;
            end
        end
        for i = 1:pmllen,
            for j = 1:hh-1,
                P.epr(i,j) = B.epr(i) * P.epr(i,j) + ...
                    A.epr(i) * ( (r(i+rshift)+dr(i+rshift)/2) * Ht(i+rshift,j) - (r(i+rshift)-dr(i+rshift-1)/2) * Ht(i+rshift-1,j) ) / ((dr(i+rshift)+dr(i+rshift-1))/2);
            end
        end
        for i = 1:rr-1,
            for j = 1:pmllen,
                P.ept(i,j) = B.ept(j) * P.ept(i,j) + ...
                    A.ept(j) * ( Hr(i,j+hshift) - Hr(i,j+hshift-1) ) / dth;
            end
        end
        
    end
    
    
    
    
    % Er update ---------------------------------
    for i = 1:rr-1,
        for j = 2:hh-1,
            Er(i,j) = c1e * Er(i,j) + c2e / (r(i)+dr(i)/2) / sin(th(j)) * ...
                ( sin(th(j)+dth/2) * Hp(i,j) - sin(th(j)-dth/2) * Hp(i,j-1)) / dth - c2e * (Jr(i+1,j) + Jr(i,j)) / 2;
        end
    end
    
    % on-axis correction
    
    for i = 1:rr-1,
        Er(i,1) = c1e * Er(i,1) + c2e * (4/dth) / (r(i)+dr(i)/2) * Hp(i,1) - c2e * (Jr(i+1,1) + Jr(i,1)) / 2;
        if thmax == pi,
            % need correction at south pole as well.
            Er(i,end) = c1e * Er(i,end) - c2e * sin(dth/2)/(r(i)+dr(i)/2)/(1 - cos(dth/2)) * Hp(i,1) - c2e * (Jr(i+1,end) + Jr(i,end)) / 2;
        end
    end
    
    % Er source
    for i = 1:2*sa,
        for j = 1:sa,
            if i < sa,
                Js(i,j) = Jin(t) * exp(-j^2/9) / (pi * 9 * dr1^2);
            else
                Js(i,j) = Jin(t) * exp(-j^2/9) * exp(-(i-sa)^2/9) / (pi * 9 * dr1^2);
            end
            Er(i,j) = Er(i,j) - c2e * Js(i,j);
        end
    end

    
    if dopml,
        % PML correction
        for i = 1:rr-1,
            for j = 1:pmllen,
                Er(i,j+hshift) = Er(i,j+hshift) + c2e / (r(i)+dr(i)/2) / sin(th(j+hshift)) * P.ert(i,j);
            end
        end
    end
    
    
    % Et update ---------------------------------------------
    for i = 2:rr-1,
        for j = 1:hh-1,
            Et(i,j) = c1e * Et(i,j) - c2e / r(i) * ...
                ( (r(i)+dr(i)/2) * Hp(i,j) - (r(i)-dr(i-1)/2) * Hp(i-1,j) ) / ((dr(i)+dr(i-1))/2) - c2e * (Jt(i,j+1) + Jt(i,j)) / 2;
        end
    end
    % PML correction
    if dopml,
        for i = 1:pmllen,
            for j = 1:hh-1,
                Et(i+rshift,j) = Et(i+rshift,j) - c2e / r(i+rshift) * P.etr(i,j);
            end
        end
    end
    
    % Ep update ------------------------------------
    for i = 2:rr-1,
        for j = 2:hh-1,
            Ep(i,j) = c1e * Ep(i,j) + c2e / r(i) * ...
                ( ( (r(i)+dr(i)/2) * Ht(i,j) - (r(i)-dr(i-1)/2) * Ht(i-1,j) ) / ((dr(i)+dr(i-1))/2) - ( Hr(i,j) - Hr(i,j-1) ) / dth ) - c2e * Jp(i,j);
        end
    end
    % PML correction
    if dopml,
        for i = 1:pmllen,
            for j = 1:hh-1,
                Ep(i+rshift,j) = Ep(i+rshift,j) + c2e / r(i+rshift) * P.epr(i,j);
            end
        end
        for i = 1:rr-1,
            for j = 1:pmllen,
                Ep(i,j+hshift) = Ep(i,j+hshift) - c2e / r(i) * P.ept(i,j);
            end
        end
    end
    
    % Ephi on the boundary using first-order Mur
    
%     if ~dopml,
%        
%         for i = 1:rr,
%            Ep(i,hh) = Ep_old(i,hh-1) + (vp*dt - r(i)*dth)/(vp*dt + r(i)*dth) * ...
%                ( Ep(i,hh-1) - Ep_old(i,hh) );
%         end
%         for i = 1:rr-1,
%            Er(i,hh) = Er_old(i,hh-1) + (vp*dt - r(i)*dth)/(vp*dt + r(i)*dth) * ...
%                ( Er(i,hh-1) - Er_old(i,hh) );
%         end
%         
%         for j = 1:hh,
%             Ep(rr,j) = Ep_old(rr-1,j) + (vp*dt - dr(rr-1))/(vp*dt + dr(rr-1)) * ...
%                 ( Ep(rr-1,j) - Ep_old(rr,j) );
%         end
%         for j = 1:hh-1,
%             Et(rr,j) = Et_old(rr-1,j) + (vp*dt - dr(rr-1))/(vp*dt + dr(rr-1)) * ...
%                 ( Et(rr-1,j) - Et_old(rr,j) );
%         end
%         
%     end
    
    
    % Update Electric field magnitude at each point --------------
    
    for j = 1:hh,
        for i = stepi:rr,
            if i == 1,
                % lower PEC boundary, Et and Ep must be zero, so only Er
                Emag(i,j) = Er(i,j);
                Epar(i,j) = Emag(i,j);
            elseif i == rr,
                % upper PEC boundary, same
                Emag(i,j) = Er(i-1,j);
                Epar(i,j) = Emag(i,j);
            else
                % everywhere else
                if j == 1,
                    % axis: due to symmetry, only Er
                    Emag(i,j) = (Er(i,j) + Er(i-1,j))/2;
                    Epar(i,j) = Emag(i,j);
                    % far PEC boundary, so only Et
                elseif j == hh,
                    Emag(i,j) = Et(i,j-1);
                    Eperp(i,j) = Emag(i,j);
                else
                    % the rest of the internal space.
                    Emag(i,j) = sqrt( ((Er(i,j) + Er(i-1,j))/2)^2 + ((Et(i,j) + Et(i,j-1))/2)^2 + Ep(i,j)^2 );
                    Epar(i,j) = (Er(i,j) + Er(i-1,j))/2;
                    Eperp(i,j) = sqrt( ((Et(i,j) + Et(i,j-1))/2)^2 + Ep(i,j)^2 );
                end
            end
            Eeff(i,j) = sqrt( Epar(i,j)^2 + Eperp(i,j)^2 * (nu(i,j)^2 / (nu(i,j)^2 + wc0^2)) );
        end
    end
    
    % -----------------------------------------------------------
    % now that I have Eeff, I can calculate vi, va, and change in electron
    % density
    if doioniz,
        vflag = 0;
        
        for i = stepi:rr,
            
            critmu = 1.36 * nd0/nd(i);
            for j = 1:hh,
                
                % mobility and collisions
                
                if (Eeff(i,j) < 1e-6),
                    mue(i,j) = critmu;
                else
                    xtemp = log10(Eeff(i,j)/nd(i));
                    mue(i,j) = 1/nd(i) * 10^(50.970 + 3.0260 * xtemp + 8.4733e-2 * xtemp^2);
                    if (mue(i,j) > critmu),
                        mue(i,j) = critmu;
                    end
                end
                nu(i,j) = (q/me) ./ mue(i,j);
                
                % ionization and optical rates
                
                xtemp = Eeff(i,j)/Ek(i);
                if xtemp < 0.2,         % Glukhov had this in here
                    vi(i,j) = 0;
                    va(i,j) = 0;
                    vN21P = 0;
                    vN22P = 0;
                else
                    va(i,j) = nd(i)/nd0 * 10^( (-0.226 * xtemp^2 + 9.332 * xtemp - 1.219)/(xtemp - 0.03009));
                    vi(i,j) = nd(i)/nd0 * 10^( (0.2018 * xtemp^2 + 10.72 * xtemp - 2.091)/(xtemp + 0.07843));
                    
                    vN21P = nd(i)/nd0 * 10^( (-0.05138 * xtemp^2 + 11.00*xtemp - 1.428)/(xtemp - 0.01424) );
                    vN22P = nd(i)/nd0 * 10^( (-0.01616 * xtemp^2 + 11.23*xtemp - 2.039)/(xtemp - 0.01504) );
                end
                
                ne2d(i,j) = ne2d(i,j) * exp((vi(i,j) - va(i,j))*dt);
                wp(i,j) = sqrt(q^2 * ne2d(i,j) / me / e0);
                
                % optics!
                
                tauN22P = 1 / (2e7 + 3e-16 * 0.21 * nd(i));
                tauN21P = 1 / (1.7e5 + 1e-17 * 0.78 * nd(i));
                
                nN22P(i,j) = tauN22P/(dt + tauN22P) * nN22P(i,j) ...
                    + (tauN22P*dt)/(dt + tauN22P) * (vN22P * ne2d(i,j));
                nN21P(i,j) = tauN21P/(dt + tauN21P) * nN21P(i,j) ...
                    + (tauN21P*dt)/(dt + tauN21P) * (vN21P * ne2d(i,j) + 2e7 * nN22P(i,j));
                
            end
        end
        
        % save the elve to process later
        
        if ~mod(t,elvemod),
            n21pcube(:,:,ii) = nN21P;
            ii = ii + 1;
        end
        
        if vflag,
            %fprintf('t = %d: WARNING! ionization coefficient is too large for time step!\n',t);
        end
        
    end
    
    % J updates ---------------------------------------------
    
    if doJ,
        
        for i = stepi:rr,
            
            if i == 1,
                Emid = Er(i,:);     % fix for bottom, perfect conductor
            elseif i == rr,
                Emid = Er(i-1,:);
            else
                Emid = ( dr(i-1) * Er(i,:) + dr(i) * Er(i-1,:) ) / (dr(i) + dr(i-1)) ;      % linear interp. across varying dr
            end
            
            for j = 1:hh-1,
                
                % update L&K matrices. I had this in a function, but it was
                % super slow! Much faster done directly here. The reason for
                % doing it at each time step, and at each i,j, is for when the
                % ionosphere parameters are time-dependent (i.e., nu and wp are
                % modified by the electric field).
                
                Ap = zeros(3,3);
                Kp = zeros(3,3);
                
                E1 = exp(-nu(i,j)*dt);
                E2 = 1/(wc0^2 + nu(i,j)^2);
                if wc0 == 0,
                    S1 = 1;
                    C1 = 0;
                else
                    S1 = sin(wc0*dt)/wc0;
                    C1 = (1 - cos(wc0*dt))/wc0^2;
                end
                C2 = (1 - E1)/nu(i,j) - E1*nu(i,j)*C1 - E1*S1;
                C3 = nu(i,j) * (1 - E1*cos(wc0*dt)) + E1*wc0*sin(wc0*dt);
                C4 = 1 - E1*cos(wc0*dt) - E1*nu(i,j)*S1;
                
                
                Ap(1,1) = E1 * ( C1*wcr^2 + cos(wc0*dt) );
                Ap(1,2) = E1 * ( C1*wcr*wct - S1*wcp );
                Ap(1,3) = E1 * ( C1*wcr*wcp + S1*wct );
                Ap(2,1) = E1 * ( C1*wct*wcr + S1*wcp );
                Ap(2,2) = E1 * ( C1*wct^2 + cos(wc0*dt) );
                Ap(2,3) = E1 * ( C1*wct*wcp - S1*wcr );
                Ap(3,1) = E1 * ( C1*wcp*wcr - S1*wct );
                Ap(3,2) = E1 * ( C1*wcp*wct + S1*wcr );
                Ap(3,3) = E1 * ( C1*wcp^2 + cos(wc0*dt) );
                
                Kp(1,1) = E2 * ( C2*wcr^2 + C3 );
                Kp(1,2) = E2 * ( C2*wcr*wct - C4*wcp );
                Kp(1,3) = E2 * ( C2*wcr*wcp + C4*wct );
                Kp(2,1) = E2 * ( C2*wct*wcr + C4*wcp );
                Kp(2,2) = E2 * ( C2*wct^2 + C3 );
                Kp(2,3) = E2 * ( C2*wct*wcp - C4*wcr );
                Kp(3,1) = E2 * ( C2*wcp*wcr - C4*wct );
                Kp(3,2) = E2 * ( C2*wcp*wct + C4*wcr );
                Kp(3,3) = E2 * ( C2*wcp^2 + C3 );
                
                % okay done with L&K matrices.
                
                
                if j == 1,      % on the axis
                    
                    if i == 1,
                        Jr(i,j) = Ap(1,1) * Jr0(i,j) + e0 * wp(i,j)^2 * ( Kp(1,1) * Er(i,j) );
                    elseif i == rr,
                        Jr(i,j) = Ap(1,1) * Jr0(i,j) + e0 * wp(i,j)^2 * ( Kp(1,1) * Er(i-1,j) );
                    else
                        Jr(i,j) = Ap(1,1) * Jr0(i,j) + ( e0 * wp(i,j)^2 * Kp(1,1) *  ( dr(i-1) * Er(i,j) + dr(i) * Er(i-1,j) ) / (dr(i) + dr(i-1)) );
                    end
                    
                else
                    
                    Jr(i,j) = Ap(1,1) * Jr0(i,j) + Ap(1,2) * Jt0(i,j) + Ap(1,3) * Jp0(i,j) ...
                        + e0 * wp(i,j)^2 * ( Kp(1,1) * Emid(j) + Kp(1,2) * 1/2 * (Et(i,j) + Et(i,j-1)) + Kp(1,3) * Ep(i,j) );
                    
                    Jt(i,j) = Ap(2,1) * Jr0(i,j) + Ap(2,2) * Jt0(i,j) + Ap(2,3) * Jp0(i,j) ...
                        + e0 * wp(i,j)^2 * ( Kp(2,1) * Emid(j) + Kp(2,2) * 1/2 * (Et(i,j) + Et(i,j-1)) + Kp(2,3) * Ep(i,j) );
                    
                    Jp(i,j) = Ap(3,1) * Jr0(i,j) + Ap(3,2) * Jt0(i,j) + Ap(3,3) * Jp0(i,j) ...
                        + e0 * wp(i,j)^2 * ( Kp(3,1) * Emid(j) + Kp(3,2) * 1/2 * (Et(i,j) + Et(i,j-1)) + Kp(3,3) * Ep(i,j) );
                    
                end
            end
        end
        
        Jr0 = Jr;
        Jt0 = Jt;
        Jp0 = Jp;
        
    end
    
    % store a value ----------------------------------
    
    sferic(t) = Er(1,round(0.5*hh));
    
    
    % plotting -------------------------------
    
    if doplot,
        
        ne2d0 = repmat(ne2',1,hh);
        deltadens = 100 * (ne2d - ne2d0) ./ ne2d0;        % relative change
        
        if ~mod(t,10),
            
            set(im1,'CData',log10(abs(Er)+1e-12));
            set(im2,'CData',log10(abs(Jr)+1e-12));
            set(im3,'CData',log10(abs(Ep)+1e-12));
            set(im4,'CData',log10(abs(nu)+1e-12));
            set(im5,'CData',deltadens);
            set(im6,'CData',log10(abs(nN21P)+1e-12));
            set(ti1,'String',sprintf('Er at time %d us',round(t*dt*1e6)));
            set(ti2,'String',sprintf('Jt at time %d us',round(t*dt*1e6)));
            set(ti3,'String',sprintf('Eeff/Ek at time %d us',round(t*dt*1e6)));
            set(ti4,'String',sprintf('nu at time %d us',round(t*dt*1e6)));
            set(ti5,'String',sprintf('Ne at time %d us',round(t*dt*1e6)));
            set(ti6,'String',sprintf('N21P at time %d us',round(t*dt*1e6)));
            
            cmax = max(max(abs(deltadens)));
            %caxis(ax5,[-cmax cmax]);
            
            drawnow;
        end
        
    end
    
    if ~mod(t,200),
        fprintf('Now %.1f%% done...\n',100*t/tsteps);
    end
   
    Ep_old = Ep;
    Er_old = Er;
    Et_old = Et;
    
end

toc;


%save n21pcube n21pcube

