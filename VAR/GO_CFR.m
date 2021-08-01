%--------------------------------------------------------------------------%
%                         Replication codes for                            %
%                   "International Credit Supply Shocks"                   %
%            by A. Cesa-Bianchi, A. Ferrero, and A. Rebucci                %
%                             December, 2017                               %
%                                                                          %
%--------------------------------------------------------------------------%
% This code replicates the results from the panel VAR model in A. 
% Cesa-Bianchi, A. Ferrero, and A. Rebucci (2018) "International Credit 
% Supply Shocks," published in the the Journal of International Economics. 
%--------------------------------------------------------------------------%
% BEFORE RUNNING THE CODE: 
% Please make sure to add the folder "Codes" and all its subfolders to your 
% Matlab path. The code has been tested with Matlab R2020a on a Mac.
%--------------------------------------------------------------------------%


%% Prelims
clear all; clear session; close all; clc
warning off all
addpath(genpath('Codes'))
savexls = 1; % if =1 saves results to excel, =0 otherwise
saveplt = 1; % if =1 saves charts to pdf, =0 otherwise

%% Codes
GO0_Data
GO1_VAR


