%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tracking continious flash suppression (tCFS) simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Christopher Whyte 02/10/25

set(0,'defaulttextinterpreter','latex')
rng('default')
close all
clear

%% parameters 

% parameter struct
p = {};

% --- simulation params
p.T = 80000; % time in ms
p.DT = 0.1;
p.sim_length = size(0:p.DT:p.T,2)-1;

% --- model parameters

% recurrent excitation (epsilon in paper)
r = .05; 
% competative inhibition
a = 3.4;

p.W = [r, -a; 
       -a, r];

% gain of population transfer function
p.M = 1;

% time constants
p.tau = 15; p.tau_H = 1000;

% adaptation strength 
gL = 1.7; gR = 3; 
p.g = [gL,gR]';

% noise in adaptation current
p.sigma = 0;

% error bound on "perceptual switches"
p.percept_bound = 0;

% rate of contrast change 
n_cont = 30; % number of intervals between max and min contrast rate
contrast_range = linspace(0.000021,0.000063,n_cont)*p.DT;

%% Run simulation

% initialise storage arrays
X_store = zeros(2,p.sim_length,n_cont);
H_store = zeros(2,p.sim_length,n_cont);
input_store = zeros(2,p.sim_length,n_cont);
percept = zeros(2,p.sim_length,n_cont);

for c_idx = 1:n_cont
    disp(['simulation: ', num2str(c_idx)]);
    p.contrast_rate = contrast_range(c_idx);
    [X_store(:,:,c_idx), H_store(:,:,c_idx), input_store(:,:,c_idx), percept(:,:,c_idx)] = tCFS_Simulator(p);
end 

%% find break through and re-suppression times

for c_idx = 1:n_cont
    BT_idx = find(diff(percept(2,:,c_idx))==1);
    ST_idx = find(diff(percept(1,:,c_idx))==1);
    % remove values whilst model is converging to equilibrium
    R_BT(c_idx) = mean(input_store(2,BT_idx(4:end-1),c_idx),2);
    R_ST(c_idx) = mean(input_store(2,ST_idx(4:end-1),c_idx),2);
end 

%% Figures

time = linspace(0,p.T/1000,p.sim_length);

% contrast rate index for figure 1
c_idx = 15;

% exzample synaptic drive for figure 1
synaptic_drive = p.W*X_store(:,:,c_idx) + input_store(:,:,c_idx) - p.g.*H_store(:,:,c_idx);

% example breakthrough times for figure 1
BT_idx = find(diff(percept(2,:,c_idx))==1);
ST_idx = find(diff(percept(1,:,c_idx))==1);
BT = time(BT_idx);
ST = time(ST_idx);

f=figure(1);
f.Position = [0,0,550,700];

subplot(4,1,1); hold on
plot(time,input_store(2,:,c_idx),'b','LineWidth',2)
xline(BT,'b--','LineWidth',1.5)
xline(ST,'m--','LineWidth',1.5)
ylabel('Contrast (a.u.)')
xlabel('Time (ms)')
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 10;
ax.LineWidth = 1.5;
% ylim([.74 .9])

subplot(4,1,2); hold on
plot(time,X_store(1,:,c_idx),'m','LineWidth',2)
plot(time,X_store(2,:,c_idx),'b','LineWidth',2)
ylabel('Firing rate (a.u.)')
xlabel('Time (s)')
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 10;
ax.LineWidth = 1.5;
xlim([0,p.T/1000])

subplot(4,1,3); hold on
plot(time,gL*H_store(1,:,c_idx),'m','LineWidth',2)
plot(time,gR*H_store(2,:,c_idx),'b','LineWidth',2)
ylabel('Adaptation current (a.u.)')
xlabel('Time (s)')
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 10;
ax.LineWidth = 1.5;
xlim([0,p.T/1000])

subplot(4,1,4); hold on
plot(time,synaptic_drive(1,:),'m','LineWidth',2)
plot(time,synaptic_drive(2,:),'b','LineWidth',2)
ylabel('synaptic drive (a.u.)')
xlabel('Time (ms)')
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 10;
ax.LineWidth = 1.5;
xlim([0,p.T/1000])

% suppression depth line plot
SD = R_BT - R_ST;

figure(2); hold on
plot(contrast_range,SD,'k-','LineWidth',3)
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax = gca;
ax.FontSize = 15;
ax.LineWidth = 1.5;
xlabel("Rate of contrast change (1/ms)")
ylabel("Suppression depth (a.u.)")
xline([contrast_range(10),contrast_range(20),contrast_range(30)],'k--','LineWidth',2)
axis padded
ylim([0 0.2])

% suppression depth bar plot for comparison with empirical data
suppression_bar = [SD(1), SD(15), SD(30)];
contrast_x = [contrast_range(10),contrast_range(20),contrast_range(30)];

figure(3);
x = ["Slow","Mid","Fast"];
bar(x,suppression_bar,'FaceColor',[.7 .7 .7],'EdgeColor',[0 0 0])
ylim([0 0.2])
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 15;
ax.LineWidth = 1.5;
xlabel("Rate of contrast change (a.u.)")
ylabel("Suppression depth (a.u.)")
box off

