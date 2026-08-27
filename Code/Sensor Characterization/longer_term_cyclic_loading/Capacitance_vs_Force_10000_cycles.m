% =========================================================================
% MATLAB Script to Synchronize and Plot Sensor Response vs. Force (Channel 7)
% Linear two-point time calibration + capacitance conversion
% =========================================================================
clear; clc; close all;

%% 1. Define File Names
force_file = "Cyclic measurement of sensor_test_8.txt"; % Zwick force file
cap_file   = "Sensor_response_test_8_10000_cycles.txt";  % Sensor response file
T3_idx = 7; % Taxel / Channel 7

%% 2. Load Zwick Force Data Safely
if ~isfile(force_file)
    error("Force file '%s' not found.", force_file);
end

fid = fopen(force_file, 'r');
if fid == -1
    error('Could not open file: %s', force_file);
end

lineIdx = 0;
headerLineNum = -1;
while ~feof(fid)
    thisLine = fgetl(fid);
    lineIdx = lineIdx + 1;
    if contains(thisLine, 'Test time')
        headerLineNum = lineIdx;
    end
end
fclose(fid);

if headerLineNum == -1
    error('Could not find header row containing "Test time" in %s.', force_file);
end
dataStartLine = headerLineNum + 1;

opts = delimitedTextImportOptions('Delimiter', ';');
opts.DataLines = [dataStartLine, Inf];
opts.VariableNames = {'TestTime', 'StandardForce', 'Time', ...
    'StandardTravel', 'AbsCrossheadTravel', 'CycleNumber'};
opts.VariableTypes = repmat({'double'}, 1, 6);
opts.ConsecutiveDelimitersRule = 'join';

T = readtable(force_file, opts);
T = rmmissing(T, 'DataVariables', {'TestTime', 'StandardForce'});

if isempty(T)
    error("Force data is empty after stripping headers. Check delimiter or file formatting.");
end

force_time = T.TestTime;
force_val  = T.StandardForce;

%% 3. Load Sensor Response Data Safely (raw ADC counts, no conversion yet)
if ~isfile(cap_file)
    error("Sensor response file '%s' not found.", cap_file);
end
cap_raw = readmatrix(cap_file);

if size(cap_raw, 2) == 17
    cap_time_raw = cap_raw(:, 17) / 1000; % ms to seconds
    cap_data = cap_raw(:, 1:16);
else
    t_max = max(force_time);
    if t_max <= 0
        error("Maximum force time is <= 0. Check Zwick file import.");
    end
    % Initial (uncorrected) assumed-constant-rate time axis
    cap_time_raw = linspace(0, t_max, size(cap_raw, 1))';
    cap_data = cap_raw(:, 1:16);
end

%% 3b. Linear Two-Point Time Calibration
% Known corresponding peak pairs (raw sensor time -> true force time)
cap_peak1_raw = 54.88;      % sensor time of 1st peak
cap_peak1_true = 1.86;      % force time of 1st peak

cap_peakN_raw = 12732.7;    % sensor time of last peak
cap_peakN_true = 12934.3;   % force time of last peak

% Solve true_time = a_time * cap_time_raw + b_time
a_time = (cap_peakN_true - cap_peak1_true) / (cap_peakN_raw - cap_peak1_raw);
b_time = cap_peak1_true - a_time * cap_peak1_raw;

fprintf('Time calibration: true_time = %.6f * raw_time + %.6f\n', a_time, b_time);
fprintf('  (slope a_time corrects sensor sample-rate error, intercept b_time corrects offset)\n');

% Apply correction to get calibrated sensor time axis
cap_time = a_time * cap_time_raw + b_time;

%% 4. Convert ADC Counts to Capacitance (pF)
% Linear calibration coefficients (ADC counts -> pF)
m = -0.0017790184866274645;
b = 16.904751517474693;

% Full affine conversion -> absolute capacitance (this dataset has no
% pre-subtracted baseline, unlike the earlier experiment, so the
% intercept b must be included here)
cap_pF_absolute = m * cap_data + b;

% Baseline-relative capacitance (equivalent to old "C - C0" convention)
% Estimated from the first N samples, assumed to be at rest / no load
baseline_window = 1:100; % adjust if the test doesn't start at rest
cap_baseline = mean(cap_pF_absolute(baseline_window, :), 1);
cap_pF_relative = cap_pF_absolute - cap_baseline;

%% 5. Synchronize onto Common Time Grid (driven by sensor response range)
time_step = 0.02; % 50 Hz grid
t_start = min(cap_time);
t_end   = max(cap_time);
time_sync = (t_start : time_step : t_end)';

% Median filter channel 7 to remove spikes (choose absolute or relative)
cap_filtered = medfilt1(cap_pF_relative(:, T3_idx), 3);

% Force is extrapolated linearly outside its own real range —
% treat those extrapolated regions with caution, they're not measured data
force_interp = interp1(force_time, force_val, time_sync, 'linear', NaN);
cap_interp   = interp1(cap_time, cap_filtered, time_sync, 'linear', 'extrap');

%% 6. Time Domain Synchronization Plot
blue   = [0/255, 114/255, 178/255];
orange = [230/255, 159/255, 0/255];

figure('Name', 'Synchronization Check', 'Color', 'w', 'Position', [200, 200, 900, 400]);
yyaxis left
plot(time_sync, cap_interp * 1000, '-', 'Color', blue, 'LineWidth', 1.2);
ylabel('Capacitance, C - C_0 [fF]');
set(gca, 'YColor', blue);

yyaxis right
plot(time_sync, force_interp, '-', 'Color', orange, 'LineWidth', 1.2);
ylabel('Force [N]');
set(gca, 'YColor', orange);

xlabel('Time [s]');
title(sprintf('Channel %d: Synchronized Capacitance and Force (Time-Calibrated)', T3_idx));
grid on;

%% 7. Hysteresis Loop: Capacitance vs. Force, Colored by Cycle
% Interpolate cycle number onto the same time_sync grid
% Use 'nearest' since cycle number is a discrete label, not continuous
cycle_interp = interp1(force_time, T.CycleNumber, time_sync, 'nearest', NaN);

unique_cycles = unique(cycle_interp(~isnan(cycle_interp)));
num_cycles = numel(unique_cycles);

figure('Name', 'Hysteresis Loop', 'Color', 'w', 'Position', [200, 200, 700, 550]);
hold on;

cmap = turbo(num_cycles); % color gradient across cycles

for k = 1:num_cycles
    this_cycle = unique_cycles(k);
    idx = (cycle_interp == this_cycle);

    plot(force_interp(idx), cap_interp(idx) * 1000, '-', ...
        'Color', cmap(k, :), 'LineWidth', 0.8);
end

hold off;
grid on;
xlabel('Applied Force [N]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Capacitance, C - C_0 [fF]', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Channel %d: Capacitance-Force Hysteresis (%d cycles)', T3_idx, num_cycles));

% Colorbar to show cycle progression
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Cycle Number';
clim([min(unique_cycles), max(unique_cycles)]);

%% 8. Peak-to-Peak Amplitude vs. Time (Cycle-Based Method)
% For each labeled cycle, take the max and min of force and capacitance
% within that cycle's window, and compute peak-to-peak (max - min).
% This uses the Zwick's own ground-truth cycle boundaries (via
% cycle_interp, already computed in Section 7) rather than inferring
% peaks from signal shape.

num_valid_cycles = numel(unique_cycles);

p2p_force = nan(num_valid_cycles, 1);
p2p_cap   = nan(num_valid_cycles, 1);
cycle_time = nan(num_valid_cycles, 1); % representative time for each cycle (midpoint)

for k = 1:num_valid_cycles
    this_cycle = unique_cycles(k);
    idx = (cycle_interp == this_cycle);

    if ~any(idx)
        continue
    end

    % Force peak-to-peak within this cycle
    f_vals = force_interp(idx);
    f_vals = f_vals(~isnan(f_vals));
    if ~isempty(f_vals)
        p2p_force(k) = max(f_vals) - min(f_vals);
    end

    % Capacitance peak-to-peak within this cycle
    c_vals = cap_interp(idx) * 1000; % convert to fF, matching your plots
    c_vals = c_vals(~isnan(c_vals));
    if ~isempty(c_vals)
        p2p_cap(k) = max(c_vals) - min(c_vals);
    end

    % Use the midpoint time of this cycle's samples as its x-axis position
    t_vals = time_sync(idx);
    cycle_time(k) = mean(t_vals);
end

%% 8b. Plot Peak-to-Peak vs. Cycle Number
figure('Name', 'Peak-to-Peak vs Cycle Number', 'Color', 'w', 'Position', [200, 200, 900, 450]);

yyaxis left
plot(unique_cycles, p2p_cap, '-', 'Color', blue, 'LineWidth', 1.0);
ylabel('Capacitance Peak-to-Peak [fF]');
set(gca, 'YColor', blue);

yyaxis right
plot(unique_cycles, p2p_force, '-', 'Color', orange, 'LineWidth', 1.0);
ylabel('Force Peak-to-Peak [N]');
set(gca, 'YColor', orange);

xlabel('Cycle Number');
title(sprintf('Channel %d: Peak-to-Peak Amplitude vs. Cycle Number (%d cycles)', T3_idx, num_valid_cycles));
grid on;

%% 9. Synchronized Time Series for First 10 Cycles, Split into Loading/Unloading
first_n_cycles = unique_cycles(unique_cycles <= 10);
num_first_cycles = numel(first_n_cycles);

cmap10 = lines(num_first_cycles); % distinguishable color per cycle

figure('Name', 'First 10 Cycles: Loading vs Unloading', 'Color', 'w', 'Position', [200, 200, 1000, 500]);

% --- Capacitance (left axis) ---
yyaxis left
hold on;
for k = 1:num_first_cycles
    this_cycle = first_n_cycles(k);
    idx = find(cycle_interp == this_cycle);
    if isempty(idx)
        continue
    end

    t_vals = time_sync(idx);
    c_vals = cap_interp(idx) * 1000;
    f_vals = force_interp(idx);

    [~, peak_idx] = max(f_vals);
    loading   = 1:peak_idx;
    unloading = peak_idx:numel(f_vals);

    % Loading: solid line
    plot(t_vals(loading), c_vals(loading), '-', ...
        'Color', cmap10(k, :), 'LineWidth', 1.2);
    % Unloading: dotted line, same color
    plot(t_vals(unloading), c_vals(unloading), ':', ...
        'Color', cmap10(k, :), 'LineWidth', 1.6);
end
ylabel('Capacitance, C - C_0 [fF]');
set(gca, 'YColor', 'k');

% --- Force (right axis) ---
yyaxis right
hold on;
for k = 1:num_first_cycles
    this_cycle = first_n_cycles(k);
    idx = find(cycle_interp == this_cycle);
    if isempty(idx)
        continue
    end

    t_vals = time_sync(idx);
    f_vals = force_interp(idx);

    [~, peak_idx] = max(f_vals);
    loading   = 1:peak_idx;
    unloading = peak_idx:numel(f_vals);

    % Loading: dashed line
    plot(t_vals(loading), f_vals(loading), '--', ...
        'Color', cmap10(k, :), 'LineWidth', 1.2);
    % Unloading: dash-dot line, same color
    plot(t_vals(unloading), f_vals(unloading), '-.', ...
        'Color', cmap10(k, :), 'LineWidth', 1.6);
end
ylabel('Force [N]');
set(gca, 'YColor', 'k');

% --- Mark cycle boundaries with vertical lines ---
for k = 1:num_first_cycles
    this_cycle = first_n_cycles(k);
    idx = (cycle_interp == this_cycle);
    t_vals = time_sync(idx);
    if ~isempty(t_vals)
        xline(t_vals(1), ':', sprintf('C%d', this_cycle), ...
            'Color', [0.5 0.5 0.5], 'LabelOrientation', 'horizontal', 'FontSize', 7);
    end
end

hold off;
grid on;
xlabel('Time [s]');
title(sprintf('Channel %d: First %d Cycles - Capacitance and Force, Split by Loading (solid/dashed) vs Unloading (dotted/dash-dot)', ...
    T3_idx, num_first_cycles));

%% 10. Peak-to-Peak Using Closest Available Value to Target Force (Loading vs Unloading)
target_high = 1.0;   % N
target_low  = 0.01;  % N
max_gap_allowed = 0.05; % N - flag if closest match is too far from target (adjust as needed)

num_valid_cycles = numel(unique_cycles);

% Preallocate: columns = [value_at_high, value_at_low], rows = cycles
cap_load   = nan(num_valid_cycles, 2);
force_load = nan(num_valid_cycles, 2);
cap_unload   = nan(num_valid_cycles, 2);
force_unload = nan(num_valid_cycles, 2);

n_flagged_load_high = 0; n_flagged_load_low = 0;
n_flagged_unload_high = 0; n_flagged_unload_low = 0;

for k = 1:num_valid_cycles
    this_cycle = unique_cycles(k);
    idx = find(cycle_interp == this_cycle);
    if isempty(idx)
        continue
    end

    t_vals = time_sync(idx);
    f_vals = force_interp(idx);
    c_vals = cap_interp(idx) * 1000; % fF

    valid = ~isnan(f_vals) & ~isnan(c_vals);
    f_vals = f_vals(valid);
    c_vals = c_vals(valid);
    if numel(f_vals) < 2
        continue
    end

    [~, peak_idx] = max(f_vals);
    loading_range   = 1:peak_idx;
    unloading_range = peak_idx:numel(f_vals);

    % ---- Loading phase ----
    fL = f_vals(loading_range);
    cL = c_vals(loading_range);

    [gap_high, i_high] = min(abs(fL - target_high));
    [gap_low, i_low]    = min(abs(fL - target_low));

    force_load(k, 1) = fL(i_high);
    cap_load(k, 1)   = cL(i_high);
    force_load(k, 2) = fL(i_low);
    cap_load(k, 2)   = cL(i_low);

    if gap_high > max_gap_allowed, n_flagged_load_high = n_flagged_load_high + 1; end
    if gap_low  > max_gap_allowed, n_flagged_load_low  = n_flagged_load_low  + 1; end

    % ---- Unloading phase ----
    fU = f_vals(unloading_range);
    cU = c_vals(unloading_range);

    [gap_high, i_high] = min(abs(fU - target_high));
    [gap_low, i_low]    = min(abs(fU - target_low));

    force_unload(k, 1) = fU(i_high);
    cap_unload(k, 1)   = cU(i_high);
    force_unload(k, 2) = fU(i_low);
    cap_unload(k, 2)   = cU(i_low);

    if gap_high > max_gap_allowed, n_flagged_unload_high = n_flagged_unload_high + 1; end
    if gap_low  > max_gap_allowed, n_flagged_unload_low  = n_flagged_unload_low  + 1; end
end

% Peak-to-peak using closest-available-value method
p2p_cap_load   = cap_load(:,1)   - cap_load(:,2);
p2p_force_load = force_load(:,1) - force_load(:,2);
p2p_cap_unload   = cap_unload(:,1)   - cap_unload(:,2);
p2p_force_unload = force_unload(:,1) - force_unload(:,2);

fprintf('Loading phase: %d/%d cycles had closest match to %.2fN off by > %.3fN\n', ...
    n_flagged_load_high, num_valid_cycles, target_high, max_gap_allowed);
fprintf('Loading phase: %d/%d cycles had closest match to %.2fN off by > %.3fN\n', ...
    n_flagged_load_low, num_valid_cycles, target_low, max_gap_allowed);
fprintf('Unloading phase: %d/%d cycles had closest match to %.2fN off by > %.3fN\n', ...
    n_flagged_unload_high, num_valid_cycles, target_high, max_gap_allowed);
fprintf('Unloading phase: %d/%d cycles had closest match to %.2fN off by > %.3fN\n', ...
    n_flagged_unload_low, num_valid_cycles, target_low, max_gap_allowed);

%% 10a. Diagnostic: Cycles That Never Reach Target Force Bounds
target_high = 1.0;    % N
target_low  = 0.02;  % N

all_cycles = unique(T.CycleNumber);
num_all_cycles = numel(all_cycles);

never_reached_high = false(num_all_cycles, 1);
never_reached_low  = false(num_all_cycles, 1);

for k = 1:num_all_cycles
    this_cycle = all_cycles(k);
    idx = (T.CycleNumber == this_cycle);
    f_this_cycle = T.StandardForce(idx);

    if isempty(f_this_cycle)
        continue
    end

    never_reached_high(k) = max(f_this_cycle) < target_high;
    never_reached_low(k)  = min(f_this_cycle) > target_low;
end

n_never_high = sum(never_reached_high);
n_never_low  = sum(never_reached_low);

pct_never_high = 100 * n_never_high / num_all_cycles;
pct_never_low  = 100 * n_never_low  / num_all_cycles;

fprintf('--- Per-Cycle Force Target Diagnostic (raw data) ---\n');
fprintf('Cycles that NEVER reach >= %.3f N: %d / %d (%.2f%%)\n', ...
    target_high, n_never_high, num_all_cycles, pct_never_high);
fprintf('Cycles that NEVER reach <= %.3f N: %d / %d (%.2f%%)\n', ...
    target_low, n_never_low, num_all_cycles, pct_never_low);

% Optional: list which specific cycle numbers were affected
fprintf('Cycle numbers never reaching %.3f N: %s\n', target_high, ...
    mat2str(all_cycles(never_reached_high)'));
fprintf('Cycle numbers never reaching %.3f N: %s\n', target_low, ...
    mat2str(all_cycles(never_reached_low)'));
%% 10b. Plot: Loading Phase Peak-to-Peak vs Cycle Number
figure('Name', 'Peak-to-Peak (Loading) vs Cycle Number', 'Color', 'w', 'Position', [200, 200, 900, 450]);

yyaxis left
plot(unique_cycles, p2p_cap_load, '-', 'Color', blue, 'LineWidth', 1.0);
ylabel('Capacitance Peak-to-Peak [fF]');
set(gca, 'YColor', blue);

yyaxis right
plot(unique_cycles, p2p_force_load, '-', 'Color', orange, 'LineWidth', 1.0);
ylabel('Force Peak-to-Peak [N]');
set(gca, 'YColor', orange);

xlabel('Cycle Number');
title(sprintf('Channel %d: Loading Phase - Peak-to-Peak at Closest Match to %.2fN/%.2fN', ...
    T3_idx, target_high, target_low));
grid on;

%% 10c. Plot: Unloading Phase Peak-to-Peak vs Cycle Number
figure('Name', 'Peak-to-Peak (Unloading) vs Cycle Number', 'Color', 'w', 'Position', [200, 200, 900, 450]);

yyaxis left
plot(unique_cycles, p2p_cap_unload, '-', 'Color', blue, 'LineWidth', 1.0);
ylabel('Capacitance Peak-to-Peak [fF]');
set(gca, 'YColor', blue);

yyaxis right
plot(unique_cycles, p2p_force_unload, '-', 'Color', orange, 'LineWidth', 1.0);
ylabel('Force Peak-to-Peak [N]');
set(gca, 'YColor', orange);

xlabel('Cycle Number');
title(sprintf('Channel %d: Unloading Phase - Peak-to-Peak at Closest Match to %.2fN/%.2fN', ...
    T3_idx, target_high, target_low));
grid on;