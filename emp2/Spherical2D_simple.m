%% 2D spherical coordinates!

clear all; close all;

tic;

dopml = 1;
doplot = 1;

fprintf('Initializing parameters and fields...\n');

e0 = 8.845e-12;
u0 = 4*pi*1e-7;
sig = 0; %1e-8;
sigm = 0;               % in this version we have updated the equations to include sigma_m
vp = 1/sqrt(e0*u0);

f0 = 30e3;              % 10 kHz
lam0 = vp/f0;           % about 30 km

r0 = 6370e3;            % a 1000 km radius "earth"
r1 = r0 + 80e3;        % a 200 km altitude "ionosphere"

dr1 = 1000;
dr2 = 1000;
stepalt = 70e3;
stepi = stepalt/dr1;
r = [r0:dr1:(r0+stepalt) (r0+stepalt+dr2):dr2:r1];
rr = length(r);
dr = diff(r);

range = 150e3;          % 300 km along the ground
thmax = range / r0;     % gives max theta in radians
dth = dr1 / r0;          % steps in theta
th = 0:dth:thmax;
hh = length(th);

%dt = 1/vp/sqrt(1/dr1^2 + 1/(r0*dth/2)^2) / stabfac;      % r0*dth gives minimum dth in km
dt = 1e-6;

maxdist = max([(r1-r0) range]);

tsteps = round( 1.5 * maxdist/dt/vp );

c1e = (2*e0 - sig*dt)/(2*e0 + sig*dt);
c2e = 2*dt/(2*e0 + sig*dt);
c1h = (2*u0 - sigm*dt)/(2*u0 + sigm*dt);
c2h = 2*dt/(2*u0 + sigm*dt);

%% ----------------
%  initialize fields

Er = zeros(rr,hh);
Et = zeros(rr,hh);
Ep = zeros(rr,hh);
Hr = zeros(rr,hh);
Ht = zeros(rr,hh);
Hp = zeros(rr,hh);
Hpold = Hp;

%% initialize source

fprintf('Setting up current source...\n');

Jin = zeros(tsteps,1);
JinS = Jin;

I0 = 200e3;

taur = 10e-6;
tauf = 40e-6;

for t = 1:tsteps,
    if (t*dt < taur),
        Jin(t) = I0 * t*dt/taur;
    else
        Jin(t) = I0 * exp(-((t*dt-taur)/tauf)^2);
    end
    
    % sinusoidal source
    JinS(t) = I0 * sin(2*pi*10e3*t*dt);

end


% bad ringing: let's filter it hard at 50 kHz. sample rate is 1/dt ~ 4 MHz
fcut = 300e3;
[b,a] = butter(2,fcut*dt/2);

Jin2 = filter(b,a,Jin);

Jin2 = I0 * Jin2 / max(Jin2);

sa = round(10e3/dr1);     % source altitude in dr steps, given 6 km
Js = zeros(sa,sa);


%% surface impedance boundary condition

sig2 = 0.01;
eps2 = 4*e0;
a = sig2/eps2;
eta2 = sqrt(u0/eps2);

Ci = [1.22646e-8 2.56716e-6 1.51777e-4 4.42437e-3 6.98268e-2 0.42473];
omegai = [4.06981e-6 1.84651e-4 3.24245e-3 3.42849e-2 0.23606 0.83083];

K = a*omegai*dt;
pi3 = exp(-K);
pi1 = eta2 * (Ci./omegai) .* (1 + (exp(-K) - 1)./K);
pi2 = eta2 * (Ci./omegai) .* (1./K - exp(-K) .* (1 + 1./K));

Ai = zeros(size(Ci));
Aiold = Ai;

%% ------------------

% coordinate transformation

x = zeros(rr,hh);
y = zeros(rr,hh);
for m = 1:rr,
    for n = 1:hh,
        x(m,n) = r(m) * sin(th(n))/1000;
        y(m,n) = r(m) * cos(th(n))/1000;
    end
end

% on nansen pcolor seems to bonk
y = (r-r0)/1000;
x = th*r0/1000;


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
    
    im1 = imagesc(x,y,log10(abs(Er)+1e-6),'parent',ax1); 
    im2 = imagesc(x,y,log10(abs(Et)+1e-6),'parent',ax2);
    im3 = imagesc(x,y,log10(abs(Ep)+1e-6),'parent',ax3);
    im4 = imagesc(x,y,log10(abs(Hr)+1e-8),'parent',ax4);
    im5 = imagesc(x,y,log10(abs(Ht)+1e-8),'parent',ax5);
    im6 = imagesc(x,y,log10(abs(Hp)+1e-8),'parent',ax6);
    
    axis(ax1,'xy');
    axis(ax2,'xy');
    axis(ax3,'xy');
    axis(ax4,'xy');
    axis(ax5,'xy');
    axis(ax6,'xy');
    
    %shading(ax1,'flat');
    %shading(ax2,'flat');
    %shading(ax3,'flat');
    %shading(ax4,'flat');
    %shading(ax5,'flat');
    %shading(ax6,'flat');
    colormap(ax1,c2);
    colormap(ax2,c2);
    colormap(ax3,c2);
    colormap(ax4,c2);
    colormap(ax5,c2);
    colormap(ax6,c2);
    caxis(ax1,[-2 6]);
    caxis(ax2,[-2 6]);
    caxis(ax3,[-2 6]);
    caxis(ax4,[-4 4]);
    caxis(ax5,[-4 4]);
    caxis(ax6,[-4 4]);
    xlabel(ax4,'Range along ground (km)');
    ylabel(ax1,'Altitude (km)');
    ylabel(ax4,'Altitude (km)');
    ti1 = title(ax1,'Er at time 0 us');
    ti2 = title(ax2,'Et at time 0 us');
    ti3 = title(ax3,'Ep at time 0 us');
    ti4 = title(ax4,'Hr at time 0 us');
    ti5 = title(ax5,'Ht at time 0 us');
    ti6 = title(ax6,'Hp at time 0 us');
    
end

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
            Hp(i,j) = c1h * Hpold(i,j) - c2h / (r(i)+dr(i)/2) * ...
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
                ( sin(th(j)+dth/2) * Hp(i,j) - sin(th(j)-dth/2) * Hp(i,j-1)) / dth;
        end
    end
    
    % on-axis correction
    
    for i = 1:rr-1,
        Er(i,1) = c1e * Er(i,1) + c2e * (4/dth) / (r(i)+dr(i)/2) * Hp(i,1);
        if thmax == pi,
            % need correction at south pole as well.
            Er(i,end) = c1e * Er(i,end) - c2e * sin(dth/2)/(r(i)+dr(i)/2)/(1 - cos(dth/2)) * Hp(i,1);
        end
    end
    
    % Er source
    for i = 1:sa,
        for j = 1:sa,
            Js(i,j) = Jin(t) * exp(-j^2/9) / (pi * 9 * dr1^2) * (1 - i/sa);
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
                ( (r(i)+dr(i)/2) * Hp(i,j) - (r(i)-dr(i-1)/2) * Hp(i-1,j) ) / ((dr(i)+dr(i-1))/2);
        end
    end
    
    % do the Et boundary using SIBC
    for j = 1:hh-1,
        for m = 1:length(Ci),
            Ai(m) = - pi1(m) * Hp(1,j) - pi2(m) * Hpold(1,j) + pi3(m) * Aiold(m);
        end
        Et(1,j) = - eta2 * Hp(1,j) - sum(Ai);
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
                ( ( (r(i)+dr(i)/2) * Ht(i,j) - (r(i)-dr(i-1)/2) * Ht(i-1,j) ) / ((dr(i)+dr(i-1))/2) - ( Hr(i,j) - Hr(i,j-1) ) / dth );
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
    
    
    % update Hp and Ai
    Hpold = Hp;
    Aiold = Ai;
    

    
    % store a value ----------------------------------
    
    sferic(t) = Hp(2,151);
    
    
    % plotting -------------------------------
    
    if doplot,
        
        if ~mod(t,10),
            
            set(im1,'CData',log10(abs(Er)+1e-12));
            set(im2,'CData',log10(abs(Et)+1e-12));
            set(im3,'CData',log10(abs(Ep)+1e-12));
            set(im4,'CData',log10(abs(Hr)+1e-12));
            set(im5,'CData',log10(abs(Ht)+1e-12));
            set(im6,'CData',log10(abs(Hp)+1e-12));
            set(ti1,'String',sprintf('Er at time %d us',round(t*dt*1e6)));
            set(ti2,'String',sprintf('Et at time %d us',round(t*dt*1e6)));
            set(ti3,'String',sprintf('Ep at time %d us',round(t*dt*1e6)));
            set(ti4,'String',sprintf('Hr at time %d us',round(t*dt*1e6)));
            set(ti5,'String',sprintf('Ht at time %d us',round(t*dt*1e6)));
            set(ti6,'String',sprintf('Hp at time %d us',round(t*dt*1e6)));
            
            %cmax = max(max(abs(deltadens)));
            %caxis(ax5,[-cmax cmax]);
            
            drawnow;
        end
        
    end
    
    if ~mod(t,round(tsteps/10)),
        fprintf('Now %.1f%% done...\n',100*t/tsteps);
    end
    
end

toc;


%save n21pcube n21pcube

