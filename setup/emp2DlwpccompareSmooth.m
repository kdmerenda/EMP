%% launch batch of EMP 2D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
emp2Ddefaults_lwpc;

% anything you want to change?
inputs.submitjob = 1;
inputs.I0 = 10e3;
inputs.beta = 0.5;
inputs.Re = 6370e3;

% master directory for set of runs
%toprundir = '/Users/Bob/Documents/temp/';
toprundir = '/projects/roma8490/emp/runs/lwpcFullDaySmooth/';

% read files from UF and generate h', beta vectors

%ufdata = load('/Users/Bob/Google Drive/STOIC/UF_data/lwpc_configurations.mat');
ufdata = load('/projects/roma8490/emp/setup/lwpc_configurations.mat');

% submit jobs

for m = 1:length(ufdata.smooth.time),
    
    % need to create input vectors of h', beta, hbrange to go into 2D runs.
    
  seg = ufdata.smooth.ionosphere(m);
inputs.hbrange = seg.distance * 1e3;
inputs.hparray = seg.hprime;
inputs.betaarray = seg.beta;
    
inputs.runname = sprintf('lwpcFD%02d',m);
inputs.rundir = [toprundir inputs.runname];
    
    % launch job
    [in,jobid] = emp2Drun(inputs);
    
%drawnow;
    
end
