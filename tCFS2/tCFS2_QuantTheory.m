%% tCFS Minimal Quant Model of Visual Rivalry Theory Script

% Christopher Whyte 21/10/25 

%% ================== SETUP ==================

close all;
clear;
clc;
rng('shuffle'); 
set(0,'defaulttextinterpreter','latex');

%% ================== CONFIGURATION ==================

% Choose the analysis mode: 'sweep' or 'theory'
config.analysis_mode = 'sweep'; 

%% ================== MODEL & SIMULATION PARAMETERS ==================

% --- Shared Simulation Parameters ---
params.T = 200000;      % Total simulation time in ms
params.DT = 0.1;          % Time step in ms
params.sim_length = length(0:params.DT:params.T) - 1;

% --- Shared Model Parameters ---
params.W_r = .05; % Recurrent excitation
params.W_i = 3.4; % Mutual inhibition
params.W = [params.W_r, -params.W_i; -params.W_i, params.W_r];

% Neuronal gain
params.gain = 1;

% Time constants
params.tau = 15;      % Firing rate time constant (ms)
params.tau_H = 1000;  % Adaptation time constant (ms)

% Error bound for defining a perceptual switch
params.percept_bound = 0;

% --- Mode-Specific Parameters ---
switch config.analysis_mode
    case 'sweep'
        disp('Running in PARAMETER SWEEP mode...');
        % Adaptation strength for the right eye (fixed)
        params.gR = 3;
        % Range of adaptation strengths for the left eye to test
        params.gL_range = 1:.1:3;
        % Range of mask strengths (input to the left eye) to test
        params.m_range = .1:.1:1;
        % Discrete contrast rates (slow, mid, fast)
        params.contrast_range = [0.000021, 0.000042, 0.000063];
        params.n_cont = length(params.contrast_range);
        
    case 'theory'
        disp('Running in THEORY COMPARISON mode...');
        % Fixed adaptation strengths
        params.gR = 3;
        params.gL = 1.7;
        % Fixed mask strength (input to the left eye)
        params.L = 0.8;
        % Continuous range of contrast rates
        params.n_cont = 30;
        params.contrast_range = linspace(0.000021, 0.000063, params.n_cont);
end

%% ================== INITIALISE STORAGE ARRAYS ==================

if strcmp(config.analysis_mode, 'sweep')
    X_store = zeros(2, params.sim_length, params.n_cont, length(params.gL_range), length(params.m_range));
    H_store = zeros(2, params.sim_length, params.n_cont, length(params.gL_range), length(params.m_range));
    input_store = zeros(2, params.sim_length, params.n_cont, length(params.gL_range), length(params.m_range));
    percept = zeros(2, params.sim_length, params.n_cont, length(params.gL_range), length(params.m_range));
else % 'theory' mode
    X_store = zeros(2, params.sim_length, params.n_cont);
    H_store = zeros(2, params.sim_length, params.n_cont);
    input_store = zeros(2, params.sim_length, params.n_cont);
    percept = zeros(2, params.sim_length, params.n_cont);
end

%% ================== SIMULATION LOOP ==================

% Loop through the specified parameter ranges and run simulation.
if strcmp(config.analysis_mode, 'sweep')
    total_sims = length(params.m_range) * length(params.gL_range) * params.n_cont;
    counter = 0;
    for m_idx = 1:length(params.m_range)
        sim_params = params;
        sim_params.L = params.m_range(m_idx);
        for g_idx = 1:length(params.gL_range)
            sim_params.gL = params.gL_range(g_idx);
            for c_idx = 1:params.n_cont
                counter = counter + 1;
                disp(['Running simulation: ', num2str(counter), ' of ', num2str(total_sims)]);
                
                % Set current contrast rate
                sim_params.contrast_rate = params.contrast_range(c_idx);
                
                % Run the simulation using the helper function
                [X, H, input_ts, percept_ts] = runCFSSimulation(sim_params);
                
                % Store results
                X_store(:,:,c_idx,g_idx,m_idx) = X;
                H_store(:,:,c_idx,g_idx,m_idx) = H;
                input_store(:,:,c_idx,g_idx,m_idx) = input_ts;
                percept(:,:,c_idx,g_idx,m_idx) = percept_ts;
            end
        end
    end
else % 'theory' mode
    for c_idx = 1:params.n_cont
        disp(['Running simulation: ', num2str(c_idx), ' of ', num2str(params.n_cont)]);
        
        sim_params = params;
        sim_params.contrast_rate = params.contrast_range(c_idx);
        
        % Run the simulation using the helper function
        [X, H, input_ts, percept_ts] = runCFSSimulation(sim_params);
        
        % Store results
        X_store(:,:,c_idx) = X;
        H_store(:,:,c_idx) = H;
        input_store(:,:,c_idx) = input_ts;
        percept(:,:,c_idx) = percept_ts;
    end
end

disp('All simulations complete.');

%% ================== ANALYSIS ==================

disp('Analyzing results...');

% This struct will hold all analysis results for plotting
analysis_results = struct(); 

if strcmp(config.analysis_mode, 'sweep')
    % --- Parameter Sweep Analysis ---
    mean_dur = zeros(params.n_cont, length(params.gL_range), length(params.m_range));
    for m_idx = 1:length(params.m_range)
        for g_idx = 1:length(params.gL_range)
            for c_idx = 1:params.n_cont
                % Calculate dominance durations for this parameter set
                durations = calculateDominanceDurations(squeeze(percept(:,:,c_idx,g_idx,m_idx)), params.DT);
                mean_dur(c_idx,g_idx,m_idx) = durations.mean_total_dur;
            end
        end
    end
    
    % Compare model results to empirical data to calculate a "loss" value
    slow_model = squeeze(mean_dur(1,:,:)) / 1000;
    mid_model = squeeze(mean_dur(2,:,:)) / 1000;
    fast_model = squeeze(mean_dur(3,:,:)) / 1000;
    
    % hard coded values from Alais et al., 2024
    slow_empirical = 4.61;
    mid_empirical = 3.45;
    fast_empirical = 3.08;
    
    % Total absolute difference between model and data
    total_delta = abs(fast_empirical-fast_model) + abs(mid_empirical-mid_model) + abs(slow_empirical-slow_model);
    total_delta(isnan(total_delta)) = max(total_delta,[],'all'); % Replace NaNs
    
    analysis_results.total_delta = total_delta;
    
else % 'theory' mode
    % --- Theory Comparison Analysis ---
    mean_durL = zeros(1, params.n_cont);
    mean_durR = zeros(1, params.n_cont);
    R0_BT = zeros(1, params.n_cont); % Contrast at breakthrough
    R0_ST = zeros(1, params.n_cont); % Contrast at suppression
    E_ST = zeros(2, params.n_cont);  % firing rate at suppression
    
    for c_idx = 1:params.n_cont
        % Calculate dominance durations
        durations = calculateDominanceDurations(squeeze(percept(:,:,c_idx)), params.DT);
        mean_durL(c_idx) = durations.mean_L_dur;
        mean_durR(c_idx) = durations.mean_R_dur;
        
        % Find breakthrough (BT) and suppression (ST) event indices
        correction = 40/params.DT; % Time correction in ms steps (approximately when the synaptic drive passes through zero)
        start_event = 5; % Ignore initial transient events
        BT_idx = find(diff(percept(2,:,c_idx))==1);
        ST_idx = find(diff(percept(1,:,c_idx))==1);

        
        % Calculate mean contrast at BT and ST events
        R0_BT(c_idx) = mean(input_store(2, BT_idx(start_event:end-1)-correction, c_idx), 2);
        R0_ST(c_idx) = mean(input_store(2, ST_idx(start_event:end-1)-correction, c_idx), 2);
        E_ST(:,c_idx) = mean(X_store(:,ST_idx(start_event:end-1)-correction,c_idx),2);
    end
    
    % Store simulation results for plotting
    analysis_results.sim_dur.L = mean_durL;
    analysis_results.sim_dur.R = mean_durR;
    analysis_results.suppression_depth_empirical = R0_BT - R0_ST;

    % Calculate correction value 
    params = findCorrectionFactor(params, R0_ST, E_ST);
    
    % Calculate theoretical durations 
    theory_durations = calculateTheoreticalDurations_Iterative(params);
    
    analysis_results.theory_dur = theory_durations;
end

%% ================== PLOTs ==================

disp('Generating plots...');
generatePlots(config, params, analysis_results, X_store, H_store, input_store);
disp('Done.');


%% ================== HELPER FUNCTIONS ==================


function [X_store, H_store, input_store, percept] = runCFSSimulation(params)
%runCFSSimulation Runs a single simulation of the tCFS model.
%
%   Inputs:
%   - params: A struct containing all model and simulation parameters.
%
%   Outputs:
%   - X_store: Time series of firing rates for left (1) and right (2) eyes.
%   - H_store: Time series of adaptation currents.
%   - input_store: Time series of contrast inputs.
%   - percept: Time series of the dominant percept (1 for left, 2 for right).

% Unpack necessary parameters for clarity
DT = params.DT;
sim_length = params.sim_length;
W = params.W;
tau = params.tau;
tau_H = params.tau_H;
gain = params.gain;
gL = params.gL;
gR = params.gR;
L = params.L;
percept_bound = params.percept_bound;
contrast_rate_per_ms = params.contrast_rate;

% Convert contrast rate to be per time step DT
contrast_rate = contrast_rate_per_ms * DT;

% --- Initial Conditions ---
X = zeros(2,1); % Initial firing rates
H = zeros(2,1); % Initial adaptation currents

% Set initial contrast for the right eye (target)
R = 1.2;

% --- Pre-allocate Storage ---
X_store = zeros(2, sim_length);
H_store = zeros(2, sim_length);
input_store = zeros(2, sim_length);
percept = zeros(2, sim_length);

% Firing rate function (rectified linear unit)
F = @(x, gain) gain * max(x, 0);
% F = @(x, gain) (x+sqrt(x.^2+eps))/2;

% sim loop
for t = 1:sim_length
    
    % Define the input vector for the current time step
    input = [L, R]';
    
    % dX/dt = (-X + F(W*X + input - g.*H)) / tau
    aggregate_synaptic_drive = W*X + input - [gL; gR].*H;
    DX = (-X + F(aggregate_synaptic_drive, gain)) / tau;
    
    % dH/dt = (-H + X) / tau_H
    DH = (-H + X) / tau_H;
    
    % Update state variables
    X = X + DT*DX;
    H = H + DT*DH;
    
    % Determine which eye is "dominant" and updates the
    % contrast of the right-eye (target) stimulus accordingly.
    if X(1) > X(2) + percept_bound
       % Left eye is dominant: target is suppressed, increase its contrast
       R = R + contrast_rate;
       percept(1,t) = 1;
    elseif X(2) > X(1) + percept_bound
       % Right eye is dominant: target is perceived, decrease its contrast
       R = R - contrast_rate;
       percept(2,t) = 1;
    end
    
    % Store current state
    X_store(:,t) = X;
    H_store(:,t) = H;
    input_store(:,t) = input;
end

end

% -------------------------------------------------------------------------

function results = calculateDominanceDurations(percept, DT)
%calculateDominanceDurations Computes the duration of perceptual dominance.
%   This function takes the percept time series and calculates the average
%   time the left eye was dominant, the right eye was dominant, and the
%   overall average.
%
%   Inputs:
%   - percept: A 2xT matrix indicating the dominant percept at each time step.
%   - DT: The simulation time step in ms.
%
%   Outputs:
%   - results: A struct containing mean_L_dur, mean_R_dur, and mean_total_dur.

results = struct('mean_L_dur', NaN, 'mean_R_dur', NaN, 'mean_total_dur', NaN);
n_ignore = 5; % Number of initial/final switches to ignore to avoid edge effects

% --- Calculate duration for Left Eye Dominance ---
vals = find(percept(1,:) == 1); % Find time indices where left eye is dominant
if length(vals) > 1
    a = diff(vals);
    b = find([a inf] > 1); % Find breaks in consecutive dominance
    durationL = diff([0 b]); % Calculate lengths of each dominance period
    if length(durationL) > (2 * n_ignore)
        results.mean_L_dur = mean(durationL(n_ignore:end-n_ignore)) * DT;
    end
end

% --- Calculate duration for Right Eye Dominance ---
vals = find(percept(2,:) == 1); % Find time indices where right eye is dominant
if length(vals) > 1
    a = diff(vals);
    b = find([a inf] > 1);
    durationR = diff([0 b]);
    if length(durationR) > (2 * n_ignore)
        results.mean_R_dur = mean(durationR(n_ignore:end-n_ignore)) * DT;
    end
end

% --- Calculate overall average dominance duration ---
if ~isnan(results.mean_L_dur) && ~isnan(results.mean_R_dur)
    results.mean_total_dur = 0.5 * (results.mean_L_dur + results.mean_R_dur);
end

end

% -------------------------------------------------------------------------

function params = findCorrectionFactor(params, R0_ST, E_ST)
%findCorrectionFactor Finds the optimal correction factor for the theoretical model.
%
%   Inputs:
%   - params:      The main parameters struct. Must contain fields:
%                  gR, W_r, contrast_range.
%   - R0_ST:       A vector of simulated contrast values at the point of
%                  suppression for each contrast rate.
%   - E_ST:        A 2xN matrix of simulated firing rates at suppression,
%                  where N is the number of contrast rates. The second row
%                  (E_ST(2,:)) is used.
%   - show_plot:   (Optional) A boolean (true/false) to control whether to
%                  display a plot of the loss function. Defaults to false.
%
%   Outputs:
%   - params:      The input parameters struct, updated with a new field:
%                  'correction_factor'.


% Unpack Necessary Parameters
gR = params.gR;
W_r = params.W_r;
contrast_range = params.contrast_range;

% Define Search Range and Initialize
correction_range = 400:1:1100; % The range of values to test for the factor
loss = zeros(size(correction_range)); % Pre-allocate the loss vector

% Calculate Loss for Each Candidate Factor
% Loop through each possible correction factor and calculate the resulting loss.
for i = 1:length(correction_range)
    current_factor = correction_range(i);
    
    % This equation represents the theoretically predicted firing rate of
    % the suppressed population (eye 2), which depends on the correction factor.
    predicted_E_ST2 = (R0_ST - current_factor .* contrast_range) / (1 + gR - W_r);
    
    % The loss is the sum of the absolute errors between the prediction and
    % the actual simulated firing rates (from E_ST).
    loss(i) = sum(abs(predicted_E_ST2 - E_ST(2,:)));
end

% Find the Optimal Correction Factor
% The best factor is the one that corresponds to the minimum calculated loss.
[~, min_loss_idx] = min(loss);
correction_factor = correction_range(min_loss_idx);

fprintf('Optimal correction factor found: %d\n', correction_factor);

% Add the Factor to the Params Struct
% This saves the result back into the main parameters struct for later use.
params.correction_factor = correction_factor;

end

% -------------------------------------------------------------------------

function results = calculateTheoreticalDurations_Iterative(params)
%   This function uses a self-consistent, iterative approach to solve for
%   the durations of perceptual suppression and dominance. It starts with
%   the initial conditions of the target contrast and iterates until the 
%   values converge.
%
%   It also analyzes the contribution of two components to the duration:
%   1.  A 'stationary' component
%   2.  A 'time-dependent' component
%
%   Inputs:
%   - params: A struct containing all model parameters, including a 'correction_factor'.
%
%   Outputs:
%   - results: A struct with the final theoretical durations and contribution analysis.

% Unpack Parameters from Struct
W_r = params.W_r;
W_i = params.W_i;
gain = params.gain;
tau_H = params.tau_H;
gL = params.gL;
gR = params.gR;
L = params.L;
contrast_range = params.contrast_range;
correction_factor = params.correction_factor; % This must be pre-calculated and added to params

n_cont = length(contrast_range);
n_iterations = 20; % Number of iterations to find a stable solution

% Initial Conditions for Iteration
% These initial conditions do not depend on the contrast rate (gamma).
E_L_inf = L / (1 + gL - W_r); % Steady-state firing rate of the left population
S_B0 = 1.2; % stimulus initial conditions

% Pre-allocate arrays
S_b = zeros(n_cont,1);
S_s = zeros(n_cont,1);
T_dom = zeros(n_cont,1);
T_sup = zeros(n_cont,1);
T_dom_const = zeros(n_cont, 1);
T_dom_time = zeros(n_cont, 1);
T_sup_const = zeros(n_cont, 1);
T_sup_time = zeros(n_cont, 1);

S_b(:, 1) = S_B0;
L_const_init = L - W_i * (S_B0 / (1 + gain*gR - gain*W_r));
alpha_init = W_i * gain * contrast_range' ./ (1 + gain*gR - gain*W_r);
T_dom(:, 1) = -L_const_init ./ alpha_init;

% Iterated Map Loop
% This loop iteratively solves for the breakthrough (S_b), suppression (S_s),
% dominant (T_dom), and suppressed (T_sup) values until they converge.

for n = 1:n_iterations-1
    for c_idx = 1:n_cont
        gamma = contrast_range(c_idx);
        
        % 1. Calculate the suppression contrast value based on the previous breakthrough value
        S_s(c_idx) = S_b(c_idx) - gamma * T_dom(c_idx);
        
        % 2. Calculate the time spent in the suppressed state (T_sup)
        R_const = S_s(c_idx) - W_i * E_L_inf;
        beta_sup = gR * (S_s(c_idx) / (1 + gR*gain - W_r*gain));
        arg_sup = (beta_sup / (gamma*tau_H)) * exp(R_const / (gamma*tau_H));
        
        T_sup_const(c_idx) = -R_const / gamma;
        T_sup_time(c_idx) = tau_H * lambertw(0, arg_sup);
        T_sup(c_idx) = T_sup_time(c_idx) + T_sup_const(c_idx);
        
        % 3. Calculate the new breakthrough contrast value
        S_b(c_idx) = S_s(c_idx) + gamma * T_sup(c_idx);
        
        % 4. Calculate the time spent in the dominant state (T_dom)
        correction = correction_factor * gamma;
        L_const = L - W_i * ((S_b(c_idx) - correction) / (1 + gain*gR - gain*W_r));
        alpha = W_i*gain*gamma / (1 + gain*gR - gain*W_r);
        beta_dom = gL * (L / (1 + gain*gL - gain*W_r));
        arg_dom = (beta_dom / (alpha*tau_H)) * exp(L_const / (alpha*tau_H));
        
        T_dom_const(c_idx) = -L_const / alpha;
        T_dom_time(c_idx) = tau_H * lambertw(0, arg_dom);
        T_dom(c_idx) = T_dom_time(c_idx) + T_dom_const(c_idx);
    end
end

% Final Calculations and Results

% --- Calculate relative contribution of stationary vs. time-dependent terms 

% For dominant durations
total_dom = abs(T_dom_const) + abs(T_dom_time);
avg_dom_stationary_contrib = mean(abs(T_dom_const) ./ total_dom);
avg_dom_time_dependent_contrib = mean(abs(T_dom_time) ./ total_dom);

% For suppressed durations
total_sup = abs(T_sup_const) + abs(T_sup_time);
avg_sup_stationary_contrib = mean(abs(T_sup_const) ./ total_sup);
avg_sup_time_dependent_contrib = mean(abs(T_sup_time) ./ total_sup);

% Overall average contributions
avg_stationary = 0.5 * (avg_dom_stationary_contrib + avg_sup_stationary_contrib);
avg_time_dependent = 0.5 * (avg_dom_time_dependent_contrib + avg_sup_time_dependent_contrib);


% Store results in the output struct
results.L = T_sup; % Suppressed durations
results.R = T_dom; % Dominant durations
results.suppression_depth_theory = contrast_range' .* T_sup;
results.avg_stationary_contribution = avg_stationary;
results.avg_time_dependent_contribution = avg_time_dependent;

end

% -------------------------------------------------------------------------

% Generate plots 
function generatePlots(config, params, results, X, H, I)
% generatePlots Creates all figures for the tCFS model analysis.
%   Inputs:
%   - config: The main configuration struct.
%   - params: The parameters struct.
%   - results: The analysis results struct.
%   - X, H, I: The raw simulation data (X_store, etc.).

if strcmp(config.analysis_mode, 'sweep')
    % --- Plot for Parameter Sweep Mode ---
    figure('Name', 'Parameter Sweep Loss Function');

    imagesc(params.m_range, params.gL_range, squeeze(results.total_delta));

    colormap(flipud(winter));
    cbar = colorbar;
    ylabel(cbar, 'Total Error (a.u.)', 'FontName', 'Times', 'FontSize', 14);
    ylabel('Adaptation Strength ($g_L$)');
    xlabel('Mask Strength ($L$)');
    title('Model Error vs. Parameters');

    set(gca, 'FontName', 'Times', 'FontSize', 14);
    set(gcf, 'Color', 'w');
    
else % 'theory' mode
    time_axis = linspace(0, params.T/1000, params.sim_length);
    c_idx = 15; % Index for example time series plot (mid contrast rate)

    % --- Figure 1: Example Time Series ---
    figure('Name', 'Example Simulation', 'Position', [0,0,600,800]);
    subplot(4,1,1);
    plot(time_axis, squeeze(I(2,:,c_idx)),'b','LineWidth',2);
    title('Slow Contrast Rate Example');
    ylabel('Contrast (a.u.)');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);

    subplot(4,1,2);
    plot(time_axis, squeeze(X(1,:,c_idx)),'m','LineWidth',2); hold on;
    plot(time_axis, squeeze(X(2,:,c_idx)),'b','LineWidth',2);
    ylabel('Firing rate (a.u.)');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);
    
    subplot(4,1,3);
    plot(time_axis, params.gL*squeeze(H(1,:,c_idx)),'m','LineWidth',2); hold on;
    plot(time_axis, params.gR*squeeze(H(2,:,c_idx)),'b','LineWidth',2);
    ylabel('Adaptation (a.u.)');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);
    
    subplot(4,1,4);
    net_input = params.W * squeeze(X(:,:,c_idx)) + squeeze(I(:,:,c_idx)) - [params.gL; params.gR] .* squeeze(H(:,:,c_idx));
    plot(time_axis, net_input(1,:), 'm', 'LineWidth', 2); hold on;
    plot(time_axis, net_input(2,:), 'b', 'LineWidth', 2);
    ylabel('Net Input (a.u.)');
    xlabel('Time (s)');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);
    set(gcf,'Color','w');

    % --- Figure 2 & 3: Theory vs. Simulation for Dominance Durations ---
    figure('Name', 'Dominance Duration Comparison', 'Position', [0,0,600,250]);
    plot(params.contrast_range, results.theory_dur.R/1000, 'k', 'LineWidth', 2); hold on;
    plot(params.contrast_range, results.sim_dur.R/1000, 'ko', 'LineWidth', 1, 'MarkerFaceColor','k');
    legend('Theory', 'Simulation', 'Location', 'northeast');
    ylabel('Dominant Duration (s)');
    xlabel('Contrast Rate');
    ylim([2, 4.5]);
    title('Target dominance duration');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);
    set(gcf,'Color','w');
    axis padded;
    
    figure('Name', 'Suppression Duration Comparison', 'Position', [0,0,600,250]);
    plot(params.contrast_range, results.theory_dur.L/1000, 'r', 'LineWidth', 2); hold on;
    plot(params.contrast_range, results.sim_dur.L/1000, 'ro', 'LineWidth', 1, 'MarkerFaceColor','r');
    legend('Theory', 'Simulation', 'Location', 'northeast');
    ylabel('Suppressed Duration (s)');
    xlabel('Contrast Rate');
    title('Target suppression duration');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);
    set(gcf,'Color','w');
    axis padded;
    
    % --- Figure 4: Suppression Depth ---
    figure('Name', 'Suppression Depth Comparison', 'Position', [0,600,600,250]);
    plot(params.contrast_range, results.theory_dur.suppression_depth_theory, 'k', 'LineWidth', 2); hold on;
    plot(params.contrast_range, results.suppression_depth_empirical, 'ko', 'LineWidth', 1, 'MarkerFaceColor', 'k');
    legend('Theory', 'Simulation', 'Location', 'northwest');
    ylabel('Suppression Depth (a.u.)');
    xlabel('Contrast Rate');
    title('Suppression Depth vs. Contrast Rate');
    set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1.5);
    set(gcf,'Color','w');
    axis padded;
    ylim([0, 0.2]);
end

end