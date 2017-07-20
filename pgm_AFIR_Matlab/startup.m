%Dock figures by default
    set(0,'DefaultFigureWindowStyle','docked');

%set debugger on (I always forget)
    %dbstop if error

%Set path back to default and add working directory and subdirectories
    path(pathdef)
    addpath(genpath(pwd));
    addpath(genpath('../../MatlabUtilities'));

%{        
%Add SBPOP to path
    if strcmp(computer,'PCWIN')
        addpath(genpath('C:\svn\SBPOP\SBPOP PACKAGE'));    
    elseif strcmp(computer,'GLNXA64')
        addpath(genpath('~/svn/SBPOP_Package/'));
    end

%Add Monolix to path
    if strcmp(computer,'GLNXA64')
        addpath(genpath('/CHBS/apps/Monolix/latest/Matlab/matlab'))
    end    
%}
