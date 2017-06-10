%% launch batch of EMP 2D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
emp2Ddefaults;

% anything you want to change?
inputs.submitjob = 1;
inputs.I0 = 100e3;

% master directory for set of runs
toprundir = '/home/kmerenda/simelve/runs/junk/';

% variable for batch of runs. name must match an input!
var1.name = 'I0';
var1.values = 1e3 * [50 75 100 150];

% submit jobs

for m = 1:length(var1.values),
    
    % change variables as requested
    evalstr = ['inputs.' var1.name ' = ' num2str(var1.values(m)) ';'];
    eval(evalstr);
    
    inputs.runname = [var1.name '_' sprintf('%.2g',var1.values(m))];
    inputs.runname(strfind(inputs.runname,'+')) = '';
    inputs.rundir = [toprundir inputs.runname];
    
    % launch job
    [in,jobid] = emp2Drun(inputs);
    
    drawnow;
    
end
