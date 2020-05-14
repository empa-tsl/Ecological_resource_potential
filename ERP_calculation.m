% calculting ERP (ecological resource potentials) for materials
tic
clc
clear all

%%%%%%%%%%%%%%%%%%
% Inupt parameter:
%%%%%%%%%%%%%%%%%%
P_v = 0.01;                  % probability of violation
n_runs = 100000;             % number of simulation runs for Monte Carlo
%%%%%%%%%%%%%%%%%%

addpath('functions/');

%%%%%%%%%%%%%%%%%%
% reading input files

filename = ['ERP.xlsx'];      
%figurename = ['display/ERA_' sector_name];

P_cut = P_v/2;              % cut off probability for quantile calculation

%%%%%%%%%%%%%%%%%
% Earth system boundaries
% load and generate random numbers for MC
ESB = MCrand2(xlsread(filename,'ESB','D6:G17'),n_runs);

%%%%%%%%%%%%%%%%%
% Unit impacts, imput data
UI_raw = xlsread(filename,'UI','G5:BB104'); 
UI_raw(isnan(UI_raw))=0;

n_ESB = size(UI_raw,2)/4;
n_resources = size(UI_raw,1);


LCI_uncertainty = MCrand2(xlsread(filename,'LCI_uncertainty','E5:H30'),n_runs);

%%%%%%%%%%%%%%%%%%
% Unit impacts (UI)
UI=zeros(n_resources,n_ESB,n_runs);
for i=1:n_ESB
    UI(:,i,:) = MCrand2(UI_raw(:,(1+4*(i-1)):(4+4*(i-1))),n_runs);
end
for i=1:n_runs
    UI(:,:,i) = UI(:,:,i) * diag(LCI_uncertainty(:,:,i));
end


ESB_limit = quantile(ESB,P_cut,3);
UI_limit = quantile(UI,1-P_cut,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERP calculation for every material individually
% How much of each material can be produced within all material's allocated
% ESB, if only this material would be produced?
ERP = zeros(n_resources,1);

for i=1:n_resources
    ERP_init(i,:) = ESB_limit' ./ UI_limit(i,:);
    [ERP_limit_init(i,:), limiting_boundary_init(i,:)] = min(ERP_init(i,:));
    [ERP(i,:), P_v_check(i,:), count(i,:)] = ERA_adjustment_to_P_v(ERP_limit_init(i,:), UI(i,:,:), ESB, P_v, n_ESB, n_runs);
    [~,limiting_boundary_rev(i,:)] = max(P_v_check(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save data


[~,boundary,~] = xlsread(filename,'dimension','G5:G30');
for i=1:n_resources
    limiting_boundary_name(i,1) = boundary(limiting_boundary_rev(i,1));
    limiting_boundary_name_max(i,1) = boundary(limiting_boundary_rev(i,1));
end

xlswrite(filename,ERP,'ERP','I5')
xlswrite(filename,limiting_boundary_name_max,'ERP','J5')

toc   