% PML setup for the 2D codes

function [A,B,P] = setupPML2d6f(pmlm,pmllen,dr1,dr2,dt,rr,hh)

simaxr = 1.5 * (pmlm + 1) / (150 * pi * dr2);
simaxt = 1.5 * (pmlm + 1) / (150 * pi * dr1);
kamax = 1;% + 3 * simax;
almax = 0;% + 0.01 * simax;
e0 = 8.854e-12;

st = zeros(pmllen,1);
stm = zeros(pmllen,1);
sr = zeros(pmllen,1);
srm = zeros(pmllen,1);

kt = zeros(pmllen,1);
ktm = zeros(pmllen,1);
kr = zeros(pmllen,1);
krm = zeros(pmllen,1);

at = zeros(pmllen,1);
atm = zeros(pmllen,1);
ar = zeros(pmllen,1);
arm = zeros(pmllen,1);
    

% first set, for the TM-type wave
        
    A.etr = zeros(pmllen,1);
    B.etr = ones(pmllen,1);
    A.ert = zeros(pmllen,1);
    B.ert = ones(pmllen,1);
    A.hpt = zeros(pmllen,1);
    B.hpt = ones(pmllen,1);
    A.hpr = zeros(pmllen,1); 
    B.hpr = ones(pmllen,1);
    
    % second set, for the TE-type wave
    
    A.epr = zeros(pmllen,1);
    A.htr = zeros(pmllen,1);
    A.ept = zeros(pmllen,1);
    A.hrt = zeros(pmllen,1);
    B.epr = ones(pmllen,1);
    B.htr = ones(pmllen,1);  
    B.ept = ones(pmllen,1);
    B.hrt = ones(pmllen,1);
    
    for m = 1:pmllen,
        
        st(m) = simaxt*((m-0.5)/pmllen)^pmlm;      % theta boundary is at max theta
        stm(m) = simaxt*(m/pmllen)^pmlm;
        sr(m) = simaxr*((m-0.5)/pmllen)^pmlm;      % r boundary is at max r
        srm(m) = simaxr*(m/pmllen)^pmlm;
        
        kt(m) = 1 + (kamax-1)*((m-0.5)/pmllen)^pmlm;      % theta boundary is at max theta
        ktm(m) = 1 + (kamax-1)*(m/pmllen)^pmlm;
        kr(m) = 1 + (kamax-1)*((m-0.5)/pmllen)^pmlm;      % r boundary is at max r
        krm(m) = 1 + (kamax-1)*(m/pmllen)^pmlm;
        
        at(m) = almax*((pmllen-m+0.5)/pmllen)^pmlm;      % theta boundary is at max theta
        atm(m) = almax*((pmllen-m)/pmllen)^pmlm;
        ar(m) = almax*((pmllen-m+0.5)/pmllen)^pmlm;      % r boundary is at max r
        arm(m) = almax*((pmllen-m)/pmllen)^pmlm;
        
        B.etr(m) = exp(-((sr(m)/kr(m)) + ar(m))*dt/e0);
        B.hpr(m) = exp(-((srm(m)/krm(m)) + arm(m))*dt/e0);
        B.ert(m) = exp(-((st(m)/kt(m)) + at(m))*dt/e0);
        B.hpt(m) = exp(-((stm(m)/ktm(m)) + atm(m))*dt/e0);
        
        A.etr(m) = sr(m) / (sr(m)*kr(m) + kr(m)^2*ar(m)) * (B.etr(m) - 1);
        A.hpr(m) = srm(m) / (srm(m)*krm(m) + krm(m)^2*arm(m)) * (B.hpr(m) - 1);
        A.ert(m) = st(m) / (st(m)*kt(m) + kt(m)^2*at(m)) * (B.ert(m) - 1);
        A.hpt(m) = stm(m) / (stm(m)*ktm(m) + ktm(m)^2*atm(m)) * (B.hpt(m) - 1);     
        
        B.epr(m) = exp(-((sr(m)/kr(m)) + ar(m))*dt/e0);
        B.htr(m) = exp(-((srm(m)/krm(m)) + arm(m))*dt/e0);
        B.ept(m) = exp(-((st(m)/kt(m)) + at(m))*dt/e0);
        B.hrt(m) = exp(-((stm(m)/ktm(m)) + atm(m))*dt/e0);
        
        A.epr(m) = sr(m) / (sr(m)*kr(m) + kr(m)^2*ar(m)) * (B.epr(m) - 1);
        A.htr(m) = srm(m) / (srm(m)*krm(m) + krm(m)^2*arm(m)) * (B.htr(m) - 1);
        A.ept(m) = st(m) / (st(m)*kt(m) + kt(m)^2*at(m)) * (B.ept(m) - 1);
        A.hrt(m) = stm(m) / (stm(m)*ktm(m) + ktm(m)^2*atm(m)) * (B.hrt(m) - 1);
        
    end
    
    % initialize Psi functions
    
    P.ert = zeros(rr,pmllen);    % Psi_Er_theta - modifies Er equation, using dHphi/dtheta - used to be rr-1 in length!
    P.etr = zeros(pmllen,hh);    % Psi_Etheta_r - modified Etheta equation, using dHphi/dr
    
    P.hpt = zeros(rr,pmllen);    % Psi_Hphi_theta - modifies Hphi equation, using dEr/dtheta  - used to be rr-1 in length!
    P.hpr = zeros(pmllen,hh);    % Psi_Hphi_r - modified Hphi equation, using dEtheta/dr
    
    P.epr = zeros(pmllen,hh);    % Psi_Er_theta - modifies Er equation, using dHphi/dtheta
    P.ept = zeros(rr,pmllen);    % Psi_Etheta_r - modified Etheta equation, using dHphi/dr
    
    P.hrt = zeros(rr,pmllen);    % Psi_Hphi_theta - modifies Hphi equation, using dEr/dtheta
    P.htr = zeros(pmllen,hh);    % Psi_Hphi_r - modified Hphi equation, using dEtheta/dr
    
