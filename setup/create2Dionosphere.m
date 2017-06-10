
function ne = create2Dionosphere(in)

loadconstants;

% read UF h', beta values and create ne array in 2D.

% need to create distance array
thmax = in.range / in.Re;
dth = in.drange / in.Re;
hh = round(thmax / dth) + 1;
th = 0:dth:(hh-1)*dth;
dvec = th*in.Re;

% test for comparison
dvec2 = 0:in.drange:in.range;

% I have inputs called in.hbrange, in.hparray, in.betaarray.

% initialize
ne0 = YukiIonosphere((in.r-in.Re)/1e3,in.betaarray(1),in.hparray(1));
ne02d = repmat(ne0,1,hh);

for m = 2:length(in.hbrange),
    ne1 = YukiIonosphere((in.r-in.Re)/1e3,in.betaarray(m),in.hparray(m));
    ind = find(dvec > in.hbrange(m),1,'first');
    if ~isempty(ind),
        ne02d(:,ind:end) = repmat(ne1,1,(hh-ind+1));
    end
end

ne = ne02d;
