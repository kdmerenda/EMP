function txinfo = GetTxInfo(txname)

% take in transmitter moniker (3-letter), and return info including
% transmitter latitude, longitude, frequency, and power
%
% if no inputs, return the entire matrix of transmitters. This is useful if
% we want to plot them on a map with labels.

txtag = [{'RA1'} {'RA2'} {'RA3'} {'JXN'} {'VTX'} ...
    {'HW3'} {'NST'} {'GQD'} {'NWC'} {'ICV'} ...
    {'3SB'} {'HW2'} {'HWV'} {'3SA'} {'NPM'} {'HW1'} ...
    {'GBZ'} {'JJI'} {'HWU'} {'DHO'} {'NAA'} ...
    {'NLK'} {'NLM'} {'TBB'} {'NRK'} {'NAU'} ...
    {'NSC'} {'VLF'} {'TBA'}];

% latitude north
txlat = [55.760 45.403 50.070 59.910 8.387 ...
    46.713 -38.481 52.911 -21.816 40.923 ...
    39.600 48.544 48.544 25.030 21.420 46.713 ...
    52.911 32.040 46.713 53.079 44.646 ...
    48.203 46.366 37.430 63.851 18.399 ...
    38.000 0 0];

% longitude east
txlon = [84.450 38.158 135.600 10.520 77.753 ...
    1.245 146.935 -3.280 114.166 9.731 ...
    103.330 2.576 2.576 111.670 -158.154 1.245 ...
    -3.280 130.810 1.245 7.614 -67.281 ...
    -121.917 -98.335 27.550 -22.459 -67.178 ...
    13.500 0 0];

% frequency in Hz
txfreq = [12000 14000 15000 16400 18200 ...
    18300 18600 19600 19800 20270 ...
    20600 20900 20900 21100 21400 21750 ...
    22100 22200 22600 23400 24000 ...
    24800 25200 26700 37500 40750 ...
    45900 10000 0];

% power in Watts
txpower = 1e3 * [500 500 500 45 NaN ...
    400 NaN 100 1000 20 ...
    NaN 400 NaN NaN 424 400 ...
    200 200 400 800 1000 ...
    192 NaN NaN NaN 100 ...
    NaN NaN NaN];

if nargin == 0,
    for m = 1:length(txtag),
        txinfo(m).tag = txtag{m};
        txinfo(m).lat = txlat(m);
        txinfo(m).lon = txlon(m);
        txinfo(m).freq = txfreq(m);
        txinfo(m).power = txpower(m);
    end
else
    
    for m = 1:length(txname),
        ind = find(strcmp(txname{m},txtag),1,'first');
        
        txinfo(m).tag = txtag{ind};
        txinfo(m).lat = txlat(ind);
        txinfo(m).lon = txlon(ind);
        txinfo(m).freq = txfreq(ind);
        txinfo(m).power = txpower(ind);
    end
    
end