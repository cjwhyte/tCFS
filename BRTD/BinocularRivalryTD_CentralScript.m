%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Binocular rivalry threshold detection simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Christopher Whyte 22/05/25

set(0,'defaulttextinterpreter','latex')
rng(0) 
close all
clear

%% Parameters 

% parameter struct
p = {};

% number of simulations 
n_sims = 40;

% --- simulation parameters
p.DT = .1;
p.burn_out = 5000/p.DT; % additional 5s of simulation time to prevent probes being presented at end simulation
p.T = 180000/p.DT + p.burn_out; % 3 min in ms + burn out
p.sim_length = length(0:p.T)-1;

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
p.tau = 15; p.tau_H = 1950;

% rivalry stimulus
L = .85; R = .85;
p.input = ones(2,p.sim_length).*[L;R];

% adaptation strength 
p.g = 2.75; 

% noise in adaptation current
p.sigma = 0.0025;

% --- probe stimulus params
t_interval = 3000/p.DT;
N_probes = (p.sim_length-p.burn_out)/t_interval;

% length of probe in ms
pulse_length = 10/p.DT;

% probe structure
probes = zeros(1,t_interval,N_probes,n_sims);
% probe onset timings structure
probe_onset = zeros(1,t_interval,N_probes,n_sims);

% randomise probe times for each simulation
for n = 1:n_sims
    probe_times = randi([1,t_interval-pulse_length],[1,N_probes]);
    for ii = 1:N_probes
        probe_time_idx = probe_times(ii);
        probes(1,probe_time_idx:probe_time_idx+pulse_length-1,ii,n) = 1;
        probe_onset(1,probe_time_idx,ii,n) = 1;
    end
end 

% reshape for use in simulations
probe_onset = [reshape(probe_onset,[1,t_interval*N_probes,n_sims]),zeros(1,p.burn_out,n_sims)];
probes = reshape(probes,[1,t_interval*N_probes,n_sims]);
probes = [probes,zeros(1,p.burn_out,n)];

% strength of probes
probe_strength_range = 0.1:0.1:1;

%% Run simulation 

% initialise storage arrays
X = zeros(2,p.T,size(probe_strength_range,2),n_sims);
H = zeros(2,p.T,size(probe_strength_range,2),n_sims);
 
sim_counter = 0;
for p_idx = 1:size(probe_strength_range,2)
    for n = 1:n_sims
        sim_counter = sim_counter + 1;
        disp(['simulation ',num2str(sim_counter), ' of ', num2str(n_sims*length(probe_strength_range))]);
        p.input = ones(2,p.sim_length).*[L;R];
        p.input(2,:) = p.input(2,:) + probe_strength_range(p_idx).*probes(:,:,n);
        [X(:,:,p_idx,n),H(:,:,p_idx,n)] = MinimalRivalry_Simulator(p);
    end 
end 

%% Calculate dominance durations 

% initialise storage 
m_durationL = zeros(size(probe_strength_range,2),n_sims);
m_durationR = zeros(size(probe_strength_range,2),n_sims);
percept = zeros(1,p.T,size(probe_strength_range,2),n_sims);
durationL = {}; durationR = {};

for p_idx = 1:size(probe_strength_range,2)
    for n = 1:n_sims
        [m_durationL(p_idx,n), m_durationR(p_idx,n),...
         durationL{p_idx,n}, durationR{p_idx,n},...
         percept(:,:,p_idx,n)] = MinimalRivalry_DominanceDurations(X(:,:,p_idx,n));
    end 
end 

%% extract and normalise dominant and suppressed time series

X_dom_FR_struct = {}; X_sup_FR_struct = {};
X_dom_probe_FR_struct = {}; X_sup_probe_FR_struct = {};

% concatenate activity, dominance, and probe info
for p_idx = 1:size(probe_strength_range,2)
    X_dom_FR = []; X_sup_FR = [];
    X_dom_probe_FR = []; X_sup_probe_FR = [];
    for n = 1:n_sims
        X_store = [X(:,:,p_idx,n); percept(:,:,p_idx,n); probe_onset(:,:,n)];
        [X_dom_FR_int, X_sup_FR_int,...
         X_dom_probe_FR_int, X_sup_probe_FR_int] = MinimalRivalry_IntervalSorting(X_store,p);
        X_dom_FR = [X_dom_FR; X_dom_FR_int];
        X_sup_FR = [X_sup_FR; X_sup_FR_int];
        X_dom_probe_FR = [X_dom_probe_FR; X_dom_probe_FR_int];
        X_sup_probe_FR = [X_sup_probe_FR; X_sup_probe_FR_int];
    end
    X_dom_FR_struct{p_idx} = X_dom_FR;
    X_sup_FR_struct{p_idx} = X_sup_FR;
    X_dom_probe_FR_struct{p_idx} = X_dom_probe_FR;
    X_sup_probe_FR_struct{p_idx} = X_sup_probe_FR;
end 

%% exctract suppression depth

p_idx = 1; 
SD_dom = []; SD_sup = [];
for n = 1:n_sims
    SD_store = r.*X(1,:,p_idx,n) + p.input(1,:) - a.*(X(2,:,p_idx,n)) - p.g*H(1,:,p_idx,n);
    X_store = [X(:,:,p_idx,n); percept(:,:,p_idx,n); probe_onset(:,:,n)];
    [SD_dom_int, SD_sup_int] = MinimalRivalry_SynIntervalSorting(X_store,SD_store,p);
    SD_dom = [SD_dom; SD_dom_int];
    SD_sup = [SD_sup; SD_sup_int];
end
SD_dom_mean = mean(SD_dom,1);
SD_sup_mean = mean(SD_sup,1);
SD_dom_ste = std(SD_dom,1);
SD_sup_ste = std(SD_sup,1);

%% compute optimal threshold across intervals for dominant and suppressed. 

% compute model response across criterions:
% - noise/stimulus absent trials are calculated from the unperturbed monocular population
% - nan entries correspond to trials where probes caused switches and
%   are therefore excluded

% dominant 
criterions_dominant = 0:.1:15;
for p_idx = 1:size(probe_strength_range,2) 
    for interval = 1:6
        for criterion = 1:length(criterions_dominant)
            signal_dominant(:,criterion,interval,p_idx) = mean(X_dom_probe_FR_struct{p_idx}(X_dom_probe_FR_struct{p_idx}(:,2)==interval,1) >  criterions_dominant(criterion),1);
            noise_dominant(:,criterion,interval,p_idx) = mean(X_dom_FR_struct{p_idx}(:,interval) >  criterions_dominant(criterion),1);
        end 
    end 
end 

% response
criterions_suppressed = 0:.1:5;
for p_idx = 1:size(probe_strength_range,2)
    for interval = 1:6
        for criterion = 1:length(criterions_suppressed)
            signal_suppressed(:,criterion,interval,p_idx) = mean(X_sup_probe_FR_struct{p_idx}(X_sup_probe_FR_struct{p_idx}(:,2)==interval,1) > criterions_suppressed(criterion),1,'omitnan');
            noise_suppressed(:,criterion,interval,p_idx) = mean(X_sup_FR_struct{p_idx}(:,interval) > criterions_suppressed(criterion),1,'omitnan');
        end 
    end 
end 

% compute misses (1-hits) and false alarms across criterions
signal_dominant = mean(mean(signal_dominant,3),4);
signal_suppressed = mean(mean(signal_suppressed,3),4);
noise_dominant = mean(mean(noise_dominant,3),4);
noise_suppressed = mean(mean(noise_suppressed,3),4);

% choose criterion that best minimises misses + false alarms
misses_plus_FA_dominant = 1-signal_dominant + noise_dominant;
misses_plus_FA_suppressed = 1-signal_suppressed + noise_suppressed;
[~,optimal_dominant_criterion] = min(misses_plus_FA_dominant,[],"All");
[~,optimal_suppressed_criterion] = min(misses_plus_FA_suppressed,[],"All");

% optimal criterion across intervals
criterion_sup = criterions_suppressed(optimal_suppressed_criterion);
criterion_dom = criterions_dominant(optimal_dominant_criterion);

% compute mean hits for each interval with optimal criterion
for p_idx = 1:size(probe_strength_range,2)
    for interval = 1:6

        hits_dominant(interval,p_idx) = mean(X_dom_probe_FR_struct{p_idx}(X_dom_probe_FR_struct{p_idx}(:,2)==interval,1) > criterion_dom,1,'omitnan');
        hits_suppressed(interval,p_idx) = mean(X_sup_probe_FR_struct{p_idx}(X_sup_probe_FR_struct{p_idx}(:,2)==interval,1) > criterion_sup,1,'omitnan');

    end 

end 

% normalise so that chance is 50%;
hits_dominant = 50 + (100-50).*hits_dominant;
hits_suppressed = 50 + (100-50).*hits_suppressed;

% compute hit standard error for each interval with optimal criterion
for p_idx = 1:size(probe_strength_range,2)
    for interval = 1:6

        hits_dominant_ste(interval,p_idx) = std(X_dom_probe_FR_struct{p_idx}(X_dom_probe_FR_struct{p_idx}(:,2)==interval,1) > criterion_dom,1,'omitnan')./sqrt(n_sims);
        hits_suppressed_ste(interval,p_idx) = std(X_sup_probe_FR_struct{p_idx}(X_sup_probe_FR_struct{p_idx}(:,2)==interval,1) > criterion_sup,1,'omitnan')./sqrt(n_sims);
  
    end 
end 

% convert to percentage
hits_dominant_ste = hits_dominant_ste * 100; 
hits_suppressed_ste = hits_suppressed_ste * 100;

%% Figures

time = linspace(0,p.T/1000*p.DT,p.sim_length);

% indices for example trial
n = 8; p_idx = 10;
p.input = ones(2,p.sim_length).*[L;R];

% synaptic_drive
synaptic_drive(1,:) = p.input(1,:) + r.*X(1,:,p_idx,n) - a.*X(2,:,p_idx,n) - p.g*H(1,:,p_idx,n);
synaptic_drive(2,:) = p.input(2,:) + r.*X(2,:,p_idx,n) + probe_strength_range(p_idx).*probes(:,:,n) - a.*X(1,:,p_idx,n) - p.g*H(2,:,p_idx,n);

f=figure(1);
subplot(2,1,1); hold on
plot(time,X(1,:,p_idx,n),'m','LineWidth',2)
plot(time,X(2,:,p_idx,n),'b','LineWidth',2)
ylabel('Firing rate')
xlabel('Time (s)')
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 20;
xlim([0,60])
ylim([0,1.25])

subplot(2,1,2); hold on
plot(time,synaptic_drive(1,:),'m','LineWidth',2)
plot(time,synaptic_drive(2,:),'b','LineWidth',2)
ylabel('Synaptic drive')
xlabel('Time (s)')
yline(0,'k--','LineWidth',2)
ax = gca;
ax.LineWidth = 1.5;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 20;
xlim([1,60])

f = figure(2); hold on 
errorbar(1:6,hits_dominant(:,p_idx),hits_dominant_ste(:,p_idx),'k-o', 'linewidth',2)
errorbar(1:6,hits_suppressed(:,p_idx),hits_suppressed_ste(:,p_idx),'k--o','linewidth',2)
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
axis padded
ylim([45,102])
yticks([60 80 100])

f.Position = [0,0,550,200];

f = figure(3); hold on
color = linspace(.4,.8,length(probe_strength_range));
errorbar(1:6,SD_dom_mean,SD_dom_ste,'k-o', 'linewidth',2)
errorbar(1:6,SD_sup_mean,SD_sup_ste,'k--o','linewidth',2)
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
axis padded
f.Position = [0,0,550,200];

f=figure(4); hold on
for ii = 1:length(hits_dominant)
    errorbar(1:6,hits_dominant(:,ii),hits_dominant_ste(:,ii),'-o', 'linewidth',2, 'Color',[1-ii/10 1-ii/10 1-ii/10])
end 
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
axis padded
f.Position = [0,0,550,200];

f = figure(5); hold on
for ii = 1:length(hits_dominant)
    errorbar(1:6,hits_suppressed(:,ii),hits_suppressed_ste(:,ii),'-o', 'linewidth',2, 'Color',[1-ii/10 1-ii/10 1-ii/10])
end
ax = gca;
set(gca, 'FontName', 'Times')
set(gcf,'Color','w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
axis padded
f.Position = [0,0,550,200];


