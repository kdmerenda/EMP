%% launch batch of EMP 3D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
emp3Ddefaults;

% master directory for set of runs
toprundir = '/home/kmerenda/simelve/runs/';

% variable for batch of runs. name must match an input!
var1.name = 'Auger3D_5kmSource_FINE_500m_100ns_SIBC';
var1.values = [100]*1e3

% submit jobs

for m = 1:length(var1.values),
    
    % change variables as requested
    evalstr = ['inputs.' var1.name ' = ' num2str(var1.values(m)) ';'];
    eval(evalstr);
    
    inputs.runname = [var1.name '_' sprintf('%03.3g',var1.values(m))];
    inputs.runname(strfind(inputs.runname,'+')) = '';
    inputs.rundir = [toprundir inputs.runname];
    
    % launch job
    [in,jobid] = emp3Drun(inputs);
    
end
