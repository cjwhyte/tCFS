%% tCFS Behavioural Data Analysis

% This script loads and processes behavioural data from three related tCFS
% (tracking Continuous Flash Suppression) experiments to analyze
% percept durations.

% data from Alais et al., 2024 can be downloaded from the OSF repo: 
% https://osf.io/nzp9v/

% Christopher Whyte 21/10/25 

%% Setup and Initialization
% -------------------------------------------------------------------------
close all;  
clear;      

% --- Analysis Parameters
start = 3; % Index to start analysis (e.g., exclude first two percepts)

% --- Experiment Parameters
% Contrast rates, converted (e.g., from %/s to %/minute or similar)
contrast_slow = 0.035 * 60;
contrast_mid = 0.07 * 60;
contrast_fast = 0.105 * 60;

% Threshold for trial exclusion in decibels
dB_threshold = 33.98;

%% Process Experiment 1 (E1)
% -------------------------------------------------------------------------
fprintf('Processing Experiment 1...\n');

% Pre-allocate matrix for all participant data
% Dimensions: [trial, percept_timestamp, subject]
times_all_1 = []; 

% Load data for each participant
for i = 1:20
    % Load the .mat file for participant 'i'
    load(sprintf("p_%d_tracking.mat", i));
    times_all_1(:,:,i) = allProfileCnts;
end 

% --- Calculate Percept Durations for E1
perceptDurs1 = []; % Initialize array to store valid durations

% Loop through each subject
for sub = 1:size(times_all_1, 3)
    trl_counter = 0; % Counter for valid trials for this subject
    
    % Loop through each trial for the current subject
    for trial = 1:size(times_all_1, 1)
        trialPs = diff([1, squeeze(times_all_1(trial, :, sub))]) ./ 60;
        
        % --- Trial Exclusion
        if max(trialPs .* contrast_mid) < dB_threshold
            trl_counter = trl_counter + 1;
            perceptDurs1 = [perceptDurs1; trialPs];
        end 
    end
end 

% --- Separate Dominant and Suppressed Percepts for E1
% Define indices for alternating percepts
TargVisible = 2:2:size(times_all_1, 2);   % Even indices: Target visible (bCFS) -> Suppressed
TargInvisible = 1:2:size(times_all_1, 2); % Odd indices: Target invisible (reCFS) -> Dominant

% Extract durations based on percept type
dominant = perceptDurs1(:, TargInvisible, :);
suppressed = perceptDurs1(:, TargVisible, :);

% Exclude the first percept (using 'start' variable)
dominant_E1 = dominant(:, start:end, :);     
suppressed_E1 = suppressed(:, start:end, :); 

%% Process Experiment 2 (E2)
% -------------------------------------------------------------------------
% This section is a repeat of E1 for a different dataset (18 subjects)
fprintf('Processing Experiment 2...\n');

times_all_2 = []; % Pre-allocate

% Load data for each participant
for i = 1:18
    load(sprintf("p_%d_tracking_Exp2.mat", i));
    times_all_2(:,:,i) = allProfileCnts;
end 

% --- Calculate Percept Durations for E2
perceptDurs2 = [];
for sub = 1:size(times_all_2, 3)
    trl_counter = 0;
    for trial = 1:size(times_all_2, 1)
        trialPs = diff([1, squeeze(times_all_2(trial, :, sub))]) ./ 60;
        
        % Trial exclusion criterion 
        if max(trialPs .* contrast_mid) < dB_threshold
            trl_counter = trl_counter + 1;
            perceptDurs2 = [perceptDurs2; trialPs];
        end 
    end
end 

% --- Separate Dominant and Suppressed Percepts for E2
TargVisible = 2:2:size(times_all_2, 2);   % bCFS (Suppressed)
TargInvisible = 1:2:size(times_all_2, 2); % reCFS (Dominant)

dominant = perceptDurs2(:, TargInvisible, :);
suppressed = perceptDurs2(:, TargVisible, :);

% Exclude the first percept
dominant_E2 = dominant(:, start:end, :);     
suppressed_E2 = suppressed(:, start:end, :); 


%% Process Experiment 3 (E3)
% -------------------------------------------------------------------------
% This section processes E3 (17 subjects), which has variable contrast rates
% per trial.
fprintf('Processing Experiment 3...\n');

times_all_3 = []; % Pre-allocate
contrast_rate = []; % Stores the contrast rate for each trial

% Load data for each participant
for i = 1:17
    load(sprintf("p_%d_tracking_Exp3.mat", i));
    % Store the contrast rate for each trial (1=slow, 2=mid, 3=fast)
    contrast_rate(:, i) = randDecrStepsIms(1, :);
    % Store the percept timestamps
    times_all_3(:,:,i) = allProfileCnts;
end 

% --- Calculate Percept Durations for E3
perceptDurs3 = [];     
percept_dur_slow = []; 
percept_dur_mid = [];  
percept_dur_fast = []; 

% Loop through each subject
for sub = 1:size(times_all_3, 3)
    % Counters for each contrast rate
    counter_slow = 0;
    counter_fast = 0;
    counter_mid = 0;
    trl_counter = 0; % Counter for all valid trials
    
    % Loop through each trial
    for trial = 1:size(times_all_3, 1)
        % Calculate percept durations (frames to seconds)
        trialPs = diff([1, squeeze(times_all_3(trial, :, sub))]) ./ 60;
        
        % --- Sort by Contrast Rate and Apply Rate-Specific Exclusion ---
        if contrast_rate(trial, sub) == 1 % Slow rate
            counter_slow = counter_slow + 1;
            percept_dur_slow(counter_slow, :, sub) = trialPs;
            if max(trialPs .* contrast_slow) < dB_threshold
                trl_counter = trl_counter + 1;
                perceptDurs3 = [perceptDurs3; trialPs];
            end 
        elseif contrast_rate(trial, sub) == 2 % Mid rate
            counter_mid = counter_mid + 1;
            percept_dur_mid(counter_mid, :, sub) = trialPs;
            if max(trialPs .* contrast_mid) < dB_threshold
                trl_counter = trl_counter + 1;
                perceptDurs3 = [perceptDurs3; trialPs];
            end 
        elseif contrast_rate(trial, sub) == 3 % Fast rate
            counter_fast = counter_fast + 1;
            percept_dur_fast(counter_fast, :, sub) = trialPs;
            if max(trialPs .* contrast_fast) < dB_threshold
                trl_counter = trl_counter + 1;
                perceptDurs3 = [perceptDurs3; trialPs];
            end 
        end 
    end
end 

% --- Separate Dominant and Suppressed Percepts for E3
TargVisible = 2:2:size(times_all_3, 2);   % bCFS (Suppressed)
TargInvisible = 1:2:size(times_all_3, 2); % reCFS (Dominant)

% Extract durations from the combined 'perceptDurs3' array
dominant = perceptDurs3(:, TargInvisible, :);
suppressed = perceptDurs3(:, TargVisible, :);

% Exclude the first percept
dominant_E3 = dominant(:, start:end, :);     
suppressed_E3 = suppressed(:, start:end, :); 


% --- Sanity Check (E3 only) ---
% Reshape rate-specific data and calculate mean median durations
percept_dur_slow = reshape(percept_dur_slow, [4*20, 17]);
percept_dur_mid = reshape(percept_dur_mid, [4*20, 17]);
percept_dur_fast = reshape(percept_dur_fast, [4*20, 17]);
slow_mean =  mean(median(percept_dur_slow, 1), 2);
mid_mean =  mean(median(percept_dur_mid, 1), 2);
fast_mean =  mean(median(percept_dur_fast, 1), 2);


%% Aggregate data from all experiments
% -------------------------------------------------------------------------
fprintf('Aggregating all experiments...\n');

% Combine all suppressed durations into a single vector
all_suppressed = [reshape(suppressed_E1, 1, []), ...
                  reshape(suppressed_E2, 1, []), ...
                  reshape(suppressed_E3, 1, [])];

% Combine all dominant durations into a single vector
all_dominant = [reshape(dominant_E1, 1, []), ...
                reshape(dominant_E2, 1, []), ... 
                reshape(dominant_E3, 1, [])];

% --- Data Filtering ---
% Exclude noise-driven reversals / very short button presses (< 0.5s)
all_suppressed = all_suppressed(all_suppressed > 0.5);
all_dominant = all_dominant(all_dominant > 0.5);


%% Statistical analysis and distribution plotting
% -------------------------------------------------------------------------
fprintf('Running statistical tests and plotting distributions...\n');

% --- KS Test ---
% Run a Kolmogorov-Smirnov test to see if the two distributions
% (dominant and suppressed) are significantly different.
[h_ks, p_ks] = kstest2(all_dominant, all_suppressed);

% --- Suppressed Durations Plot (Figure 1) ---
Nbins = round(sqrt(length(all_suppressed))); % Calculate optimal bin number
figure(1); 
hold on;

% Plot histogram with a fitted Gamma distribution
h = histfit(all_suppressed, Nbins, 'gamma');
set(h(1), 'facecolor', [76, 109, 172]/255); % Set histogram color
set(h(2), 'color', 'k');                     % Set gamma fit line color

% Plot a fitted Lognormal distribution on the same axes
h1 = histfit(all_suppressed, Nbins, 'lognormal');
set(h1(1), 'facealpha', 0); % Make lognormal histogram invisible
set(h1(2), 'color', 'k');    % Set lognormal fit line color
set(h1(2), 'linestyle', '--'); % Make lognormal fit dashed

% Normalize Y-axis to show probability density
y = get(gca, 'YTick'); 
set(gca, 'YTick', y, 'YTickLabel', round(y/numel(all_suppressed), 3));

% --- Plot Formatting ---
xlim([0.5, 12]);
ax = gca;
ax.FontSize = 20; 
set(gca, 'FontName', 'Times');

% --- Dominant Durations Plot (Figure 2) ---
Nbins = round(sqrt(length(all_dominant))); % Recalculate bins
figure(2); 
hold on;

% Plot histogram with a fitted Gamma distribution
h = histfit(all_dominant, Nbins, 'gamma');
set(h(1), 'facecolor', [76, 109, 172]/255); 
set(h(2), 'color', 'k');

% Plot a fitted Lognormal distribution
h1 = histfit(all_dominant, Nbins, 'lognormal');
set(h1(1), 'facealpha', 0); 
set(h1(2), 'color', 'k');
set(h1(2), 'linestyle', '--'); 

% Normalize Y-axis
y = get(gca, 'YTick');
set(gca, 'YTick', y, 'YTickLabel', round(y/numel(all_dominant), 3));

% --- Plot Formatting ---
xlim([0.5, 12]);
ax = gca;
ax.FontSize = 20; 
set(gca, 'FontName', 'Times');

hold off;


%% Distribution fitting and model comparison
% -------------------------------------------------------------------------
fprintf('Fitting distributions and calculating likelihoods...\n');

% --- Fit Distributions ---
% Fit Gamma distributions
p_suppresed_gamma = fitdist(all_suppressed', 'Gamma');
p_dominant_gamma = fitdist(all_dominant', 'Gamma');
% Get Gamma confidence intervals
suppressed_gamma_ci =  paramci(p_suppresed_gamma);
dominant_gamma_ci =  paramci(p_dominant_gamma);
% Get Gamma negative log-likelihoods
LL_suppressed_gamma = negloglik(p_suppresed_gamma);
LL_dominant_gamma = negloglik(p_dominant_gamma);

% Fit Lognormal distributions
p_suppresed_lnormal = fitdist(all_suppressed', 'lognormal');
p_dominant_lnormal = fitdist(all_dominant', 'lognormal');
% Get Lognormal confidence intervals
suppressed_lnormal_ci =  paramci(p_suppresed_lnormal);
dominant_lnormal_ci =  paramci(p_dominant_lnormal );
% Get Lognormal negative log-likelihoods
LL_suppressed_lognormal = negloglik(p_suppresed_lnormal);
LL_dominant_lognormal = negloglik(p_dominant_lnormal);

% Fit Normal distributions (for comparison)
p_suppresed_normal = fitdist(all_suppressed', 'normal');
p_dominant_normal = fitdist(all_dominant', 'normal');
% Get Normal negative log-likelihoods
LL_suppressed_normal = negloglik(p_suppresed_normal);
LL_dominant_normal = negloglik(p_dominant_normal);

% --- Log-Likelihood Ratio Test ---
% (To compare which fit, Gamma or Lognormal, is better)
LL_ratio_suppressed = LL_suppressed_gamma / LL_suppressed_lognormal;
LL_ratio_dominant = LL_dominant_gamma / LL_dominant_lognormal;

% Perform the formal test
[h_sup, p_sup] = lratiotest(-LL_suppressed_gamma, -LL_suppressed_lognormal, 2, .05);
[h_dom, p_dom] = lratiotest(-LL_dominant_gamma, -LL_dominant_lognormal, 2, .05);


%% Plot fitted parameters
% -------------------------------------------------------------------------
fprintf('Plotting fitted parameters...\n');

% --- Plot Gamma Parameters (Figure 3) ---
% Extract shape (a) and scale (b) parameters
params_a = [p_dominant_gamma.a, p_suppresed_gamma.a];
param_a_ci = abs([dominant_gamma_ci(1,1) - dominant_gamma_ci(2,1), ...
                suppressed_gamma_ci(1,1) - suppressed_gamma_ci(2,1)]);
params_b = [p_dominant_gamma.b, p_suppresed_gamma.b];
param_b_ci = abs([dominant_gamma_ci(1,2) - dominant_gamma_ci(2,2), ...
                suppressed_gamma_ci(1,2) - suppressed_gamma_ci(2,2)]);
            
x1 = [0.25, 0.5]; % X-positions for param 'a'
x2 = [0.75, 1];   % X-positions for param 'b'

f = figure(3); 
hold on;
% Plot parameter 'a' on the left Y-axis
errorbar(x1, params_a, param_a_ci, 'ko', 'linewidth', 1, 'LineStyle', 'none');
% Use the right Y-axis for parameter 'b'
yyaxis right;
errorbar(x2, params_b, param_b_ci, 'ko', 'linewidth', 1, 'LineStyle', 'none');
ylim([0 20]); % Set limit for right Y-axis

% --- Plot Formatting ---
ax = gca;
set(gca, 'FontName', 'Times');
set(gcf, 'Color', 'w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
axis padded;
xticks([]); % Remove x-ticks
box off;
f.Position = [0, 0, 300, 300];

% --- Plot Lognormal Parameters (Figure 4) ---
% Extract mu and sigma parameters
params_mu = [p_dominant_lnormal.mu, p_suppresed_lnormal.mu];
param_mu_ci = abs([dominant_lnormal_ci(1,1) - dominant_lnormal_ci(2,1), ...
                suppressed_lnormal_ci(1,1) - suppressed_lnormal_ci(2,1)]);
params_sigma = [p_dominant_lnormal.sigma, p_suppresed_lnormal.sigma];
param_sigma_ci = abs([dominant_lnormal_ci(1,2) - dominant_lnormal_ci(2,2), ...
                suppressed_lnormal_ci(1,2) - suppressed_lnormal_ci(2,2)]);

x1 = [0.25, 0.5]; % X-positions for param 'mu'
x2 = [0.75, 1];   % X-positions for param 'sigma'

f = figure(4); 
hold on;
% Plot parameter 'mu' on the left Y-axis
errorbar(x1, params_mu, param_mu_ci, 'ko', 'linewidth', 1, 'LineStyle', 'none');
% Use the right Y-axis for parameter 'sigma'
yyaxis right;
errorbar(x2, params_sigma, param_sigma_ci, 'ko', 'linewidth', 1, 'LineStyle', 'none');
ylim([0 20]); % Set limit for right Y-axis

% --- Plot Formatting ---
ax = gca;
set(gca, 'FontName', 'Times');
set(gcf, 'Color', 'w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
axis padded;
xticks([]); % Remove x-ticks
box off;
f.Position = [0, 0, 300, 300];

hold off;


%% Plot example trials
% -------------------------------------------------------------------------
% This section plots the raw contrast profile for two example trials
% from a single participant in Experiment 3.

fprintf('Plotting example trials...\n');

participant = 1; % Select participant 1
load(sprintf("p_%d_tracking_Exp3.mat", participant));

% --- Plot Example Fast Trial ---
fasttrial = 5; % Select trial 5
contrData = contrProfile{1, fasttrial}; % Get contrast data
time = linspace(1, times_all_3(fasttrial, end, participant) ./ 60, length(contrData));

f = figure(5); 
plot(time, 20*log10(contrData)); % Plot contrast in dB
xlim([0, 120]);

% --- Plot Formatting ---
ax = gca;
set(gca, 'FontName', 'Times');
set(gcf, 'Color', 'w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
f.Position = [0, 0, 500, 150];
box off;

% --- Plot Example Slow Trial ---
slowtrial = 4; % Select trial 4
contrData = contrProfile{1, slowtrial}; % Get contrast data
time = linspace(1, times_all_3(slowtrial, end, participant) ./ 60, length(contrData));

f = figure(6); 
plot(time, 20*log10(contrData)); % Plot contrast in dB
xlim([0, 120]);

% --- Plot Formatting ---
ax = gca;
set(gca, 'FontName', 'Times');
set(gcf, 'Color', 'w');
ax.FontSize = 20;
ax.LineWidth = 1.5;
f.Position = [0, 0, 500, 150];
box off



