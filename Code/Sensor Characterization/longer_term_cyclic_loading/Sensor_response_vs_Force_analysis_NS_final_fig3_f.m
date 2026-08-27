% =========================================================================
% Sensor Response vs. Force Analysis (Channel 7)
%
%   Figure 1 - Capacitance and force plotted together over the whole
%              test, so you can see the raw signals side by side.
%
%   Figure 2 - Zooms into just the first 10 cycles. Loading (force going
%              up) and unloading (force coming back down) are drawn
%              separately, with markers showing exactly where each cycle
%              crossed the 1N and 0.02N targets.
%
%   Figure 3 - The capacitance value measured at 1N and at 0.02N, for
%              every cycle in the test, using only the LOADING (rising
%              force) half of each cycle.
%
%   Figure 4 - Same as Figure 3, but using only the UNLOADING (falling
%              force) half of each cycle.
%
%   Figure 5 - How much the sensor's response range (high minus low,
%              loading side only) has grown or shrunk compared to the
%              start of the test, cycle by cycle.
%
%   Figure 6 - How consistent the peak-to-peak (P2P) response has been, checked
%              every 1000 cycles as more data is included (cycles 1-1000,
%              then 1-2000, then 1-3000...). Three versions shown: the
%              combined average P2P, loading only P2P, and unloading only P2P.
%
%   Figure 7 - Same idea as Figure 6, but for the capacitance reading at
%              1N during loading by itself (not the P2P loading/unloading difference) --
%              how consistent has the sensor response been at 1N as more cycles
%              get included.
%
% =========================================================================

clear; clc; close all;

%% 1. File Names & Config
force_file = "Cyclic measurement of sensor_test_8.txt";
cap_file   = "Sensor_response_test_8_10000_cycles.txt";
T3_idx = 7; % Taxel / Channel 7

target_high = 1.0;   % N
target_low  = 0.02;  % N

blue   = [0/255, 114/255, 178/255];
orange = [230/255, 159/255, 0/255];

%% 2. Load Zwick Force Data
if ~isfile(force_file)
    error("Force file '%s' not found.", force_file);
end

fid = fopen(force_file, 'r');
if fid == -1
    error('Could not open file: %s', force_file);
end
lineIdx = 0; headerLineNum = -1;
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

opts = delimitedTextImportOptions('Delimiter', ';');
opts.DataLines = [headerLineNum + 1, Inf];
opts.VariableNames = {'TestTime', 'StandardForce', 'Time', ...
    'StandardTravel', 'AbsCrossheadTravel', 'CycleNumber'};
opts.VariableTypes = repmat({'double'}, 1, 6);
opts.ConsecutiveDelimitersRule = 'join';

T = readtable(force_file, opts);
T = rmmissing(T, 'DataVariables', {'TestTime', 'StandardForce'});
if isempty(T)
    error("Force data is empty after stripping headers. Check delimiter or file formatting.");
end

force_time = T.TestTime;   % raw, native Zwick clock -- used for ALL cycle detection & crossing math
force_val  = T.StandardForce;

%% 3. Load Raw Capacitance Data
if ~isfile(cap_file)
    error("Sensor response file '%s' not found.", cap_file);
end
cap_raw = readmatrix(cap_file);

if size(cap_raw, 2) == 17
    cap_time_raw = cap_raw(:, 17) / 1000; % ms to seconds
    cap_data = cap_raw(:, 1:16);
else
    t_max = max(force_time);
    cap_time_raw = linspace(0, t_max, size(cap_raw, 1))';
    cap_data = cap_raw(:, 1:16);
end

%% 3b. Two-Point Linear Time Calibration (sensor clock -> approx. force clock)
% Only used for: (a) background curves in Figures 1 & 2, and (b) the
% single capacitance-at-crossing-time lookup in Section 6. NEVER used for
% cycle detection or crossing-time math.
cap_peak1_raw = 54.88;      cap_peak1_true = 1.86;
cap_peakN_raw = 12732.7;    cap_peakN_true = 12934.3;

a_time = (cap_peakN_true - cap_peak1_true) / (cap_peakN_raw - cap_peak1_raw);
b_time = cap_peak1_true - a_time * cap_peak1_raw;
cap_time = a_time * cap_time_raw + b_time;

%% 4. ADC Counts -> Capacitance (pF)
m = -0.0017790184866274645;
b = 16.904751517474693;

cap_pF_absolute = m * cap_data + b;
baseline_window = 1:100;
cap_baseline = mean(cap_pF_absolute(baseline_window, :), 1);
cap_pF_relative = (cap_pF_absolute - cap_baseline)./cap_baseline(:,:); % divided over baseline T3

%% 5. Resample onto a Common Grid (for background plotting only)
time_step = 0.02; % 50 Hz
t_start = min(cap_time);
t_end   = max(cap_time);
time_sync = (t_start : time_step : t_end)';

cap_filtered = medfilt1(cap_pF_relative(:, T3_idx), 3);
force_interp = interp1(force_time, force_val, time_sync, 'linear', NaN);
cap_interp   = interp1(cap_time, cap_filtered, time_sync, 'linear', 'extrap');

%% 6. Cycle Detection on RAW Force Clock (single source of truth)
min_cycle_period_s  = 1.0;   % somewhat less than the true cycle period
smooth_window_s     = 0.1;   % smoothing window in seconds (raw clock)
min_peak_prominence = 0.5 * (target_high - target_low); % rejects tiny preload wiggles

force_smooth = smoothdata(force_val, 'movmean', smooth_window_s, 'SamplePoints', force_time);

[~, peak_time_raw]   = findpeaks(force_smooth,  force_time, ...
    'MinPeakDistance', min_cycle_period_s, 'MinPeakProminence', min_peak_prominence);
[~, trough_time_raw] = findpeaks(-force_smooth, force_time, ...
    'MinPeakDistance', min_cycle_period_s, 'MinPeakProminence', min_peak_prominence);

num_cycles = max(0, numel(trough_time_raw) - 1);
cycle_num  = (1:num_cycles)';

fprintf('Detected %d cycles from raw force clock.\n', num_cycles);

% Map trough times onto the nearest time_sync index, for background-curve
% slicing in Figure 2 only (visualization, not used for any measurement).
trough_idx_sync = nan(numel(trough_time_raw), 1);
for i = 1:numel(trough_time_raw)
    [~, trough_idx_sync(i)] = min(abs(time_sync - trough_time_raw(i)));
end

%% 7. Exact Target-Crossing Per Cycle (linear interpolation, trough-centered pairing)
interp_method = 'linear';

high_force_exact = nan(num_cycles, 1);
high_cap_exact    = nan(num_cycles, 1);
low_force_exact   = nan(num_cycles, 1);
low_cap_exact     = nan(num_cycles, 1);

t_high_load_exact   = nan(num_cycles, 1);
t_high_unload_exact = nan(num_cycles, 1);
t_low_load_exact    = nan(num_cycles, 1);
t_low_unload_exact  = nan(num_cycles, 1);

cap_high_load_exact   = nan(num_cycles, 1);
cap_high_unload_exact = nan(num_cycles, 1);
cap_low_load_exact    = nan(num_cycles, 1);
cap_low_unload_exact  = nan(num_cycles, 1);

prev_tU = []; prev_fU = []; % carries the preceding cycle's unloading segment forward,
                             % so each trough's "falling in" crossing pairs with the
                             % segment that actually leads into THAT trough

for k = 1:num_cycles
    t_start_k = trough_time_raw(k);
    t_end_k   = trough_time_raw(k+1);

    this_peak = peak_time_raw(peak_time_raw > t_start_k & peak_time_raw < t_end_k);
    if isempty(this_peak)
        continue
    end
    t_peak = this_peak(1);

    mask_load   = force_time >= t_start_k & force_time <= t_peak;
    mask_unload = force_time >= t_peak    & force_time <= t_end_k;

    tL = force_time(mask_load);   fL = force_val(mask_load);
    tU = force_time(mask_unload); fU = force_val(mask_unload);

    % --- HIGH target (~1N): both crossings bracket the SAME peak ---
    [t_h_load, f_h_load]   = find_crossing(tL, fL, target_high, interp_method);
    [t_h_unload, f_h_unload] = find_crossing(tU, fU, target_high, interp_method);

    high_force_exact(k) = mean([f_h_load, f_h_unload]);
    t_high_load_exact(k)   = t_h_load;
    t_high_unload_exact(k) = t_h_unload;

    c_h_load   = interp1(cap_time, cap_pF_relative(:, T3_idx), t_h_load,   interp_method, 'extrap');
    c_h_unload = interp1(cap_time, cap_pF_relative(:, T3_idx), t_h_unload, interp_method, 'extrap');
    high_cap_exact(k) = mean([c_h_load, c_h_unload]) * 1000; % fF
    cap_high_load_exact(k)   = c_h_load * 1000;
    cap_high_unload_exact(k) = c_h_unload * 1000;

    % --- LOW target (~0.02N): TROUGH-CENTERED pairing ---
    % Rising-out-of-trough(k): this cycle's own loading segment.
    [t_l_load, f_l_load] = find_crossing(tL, fL, target_low, interp_method);

    % Falling-into-trough(k): the PRECEDING cycle's unloading segment
    % (the data that actually leads into THIS trough).
    if ~isempty(prev_fU)
        [t_l_unload, f_l_unload] = find_crossing(prev_tU, prev_fU, target_low, interp_method);
    else
        t_l_unload = t_l_load; f_l_unload = f_l_load; % first cycle: no preceding data
    end

    low_force_exact(k) = mean([f_l_load, f_l_unload]);
    t_low_load_exact(k)   = t_l_load;
    t_low_unload_exact(k) = t_l_unload;

    c_l_load   = interp1(cap_time, cap_pF_relative(:, T3_idx), t_l_load,   interp_method, 'extrap');
    c_l_unload = interp1(cap_time, cap_pF_relative(:, T3_idx), t_l_unload, interp_method, 'extrap');
    low_cap_exact(k) = mean([c_l_load, c_l_unload]) * 1000;
    cap_low_load_exact(k)   = c_l_load * 1000;
    cap_low_unload_exact(k) = c_l_unload * 1000;

    prev_tU = tU; prev_fU = fU;
end

p2p_cap_exact = high_cap_exact - low_cap_exact; % kept for reference

%% ---- Continuous "cycle number" x-axis for the whole recording ----
% Troughs delimit cycles: trough_time_raw(k) marks the start of cycle k.
% Linearly interpolate a fractional cycle coordinate for every synced sample,
% so force/capacitance can be plotted against CYCLES instead of time.
if numel(trough_time_raw) >= 2
    cycle_edge_coord = (1:numel(trough_time_raw))';               % trough k -> cycle k
    cyc_axis_sync = interp1(trough_time_raw, cycle_edge_coord, time_sync, 'linear', 'extrap');
    cyc_axis_lim  = [1, numel(trough_time_raw)];
else
    cyc_axis_sync = time_sync;                                    % fallback if too few troughs
    cyc_axis_lim  = [min(time_sync), max(time_sync)];
end

%% ===================== FIGURE 1 =====================
% Force (top) and capacitance (bottom) vs. cycle number, shared x-axis.
figure('Name', 'Figure 1: Force and Capacitance vs Cycle', 'Color', 'w', 'Position', [200, 200, 900, 600]);

ax1_f = subplot(2, 1, 1);
plot(cyc_axis_sync, force_interp, '-', 'Color', orange, 'LineWidth', 1.2);
ylabel('Force [N]');
grid on; xlim(cyc_axis_lim);
title(sprintf('Channel %d: Force and Capacitance vs. Cycle', T3_idx));
set(ax1_f, 'XTickLabel', []);   % hide top x labels (shared with bottom panel)

ax1_c = subplot(2, 1, 2);
plot(cyc_axis_sync, cap_interp * 1000, '-', 'Color', blue, 'LineWidth', 1.2);
ylabel('Capacitance, C - C_0 [fF]');
xlabel('Cycle Number');
grid on; xlim(cyc_axis_lim);

linkaxes([ax1_f, ax1_c], 'x');

%% ============ FIGURE: Capacitance vs Cycle ============
% -------------------- EDITABLE panel settings --------------------
panel_width_mm  = 55;                 % panel width  (mm)
panel_height_mm = 50;                 % panel height (mm)
panel_font      = 'Arial';
panel_fontsize  = 7;                  % pt at final print size (Nature: 5-7 pt)
data_linewidth  = 0.75;               % data line width (pt)
axis_linewidth  = 0.5;                % axis/tick line width (pt)
export_dir      = 'figures';   % output subfolder ('.' = current folder)
% -----------------------------------------------------------------

if ~isfolder(export_dir); mkdir(export_dir); end

fig_capA = figure('Name', 'Capacitance vs Cycle (Nature panel)', 'Color', 'w');
ax_capA  = axes(fig_capA);
plot(ax_capA, cyc_axis_sync, cap_interp, '-', 'Color', "#ea7878", 'LineWidth', data_linewidth);
xlim(ax_capA, cyc_axis_lim);
xlabel(ax_capA, 'Cycle number');
ylabel(ax_capA, '$\Delta C/C_0$', 'Interpreter', 'latex');
apply_nature_style(fig_capA, ax_capA, panel_width_mm, panel_height_mm, panel_font, panel_fontsize, axis_linewidth);
export_vector(fig_capA, fullfile(export_dir, sprintf('capacitance_vs_cycle_ch%d', T3_idx)));

%% ============ FIGURE(S): Close-up insets of capacitance ============
% Two small zoomed panels of the same capacitance signal, exported as vector.
%
% -------------------- EDITABLE close-up settings --------------------
closeup_centers    = [100, 10000];   % cycle numbers to zoom into
closeup_halfwidth  = 2;             % +/- this many cycles around each center
closeup_width_mm   = 55;            % close-up panel width  (mm)
closeup_height_mm  = 55;            % close-up panel height (mm)
% --------------------------------------------------------------------

for ci = 1:numel(closeup_centers)
    ctr  = closeup_centers(ci);
    x_lo = ctr - closeup_halfwidth;
    x_hi = ctr + closeup_halfwidth;

    if ctr < cyc_axis_lim(1) || ctr > cyc_axis_lim(2)
        warning('Close-up center %g is outside the detected cycle range [%g %g]; skipping.', ...
            ctr, cyc_axis_lim(1), cyc_axis_lim(2));
        continue
    end

    fig_zoom = figure('Name', sprintf('Capacitance close-up ~cycle %g', ctr), 'Color', 'w');
    ax_zoom  = axes(fig_zoom);
    plot(ax_zoom, cyc_axis_sync, cap_interp, '-',  'Color', "#ea7878", 'LineWidth', data_linewidth);
    xlim(ax_zoom, [x_lo, x_hi]);

    in_win = cyc_axis_sync >= x_lo & cyc_axis_sync <= x_hi;   % autoscale y to the window
    if any(in_win)
        yv  = cap_interp(in_win);
        pad = 0.05 * (max(yv) - min(yv) + eps);
        ylim(ax_zoom, [min(yv) - pad, max(yv) + pad]);
    end
    xlabel(ax_zoom, 'Cycle number');
    ylabel(ax_capA, '$\Delta C/C_0$', 'Interpreter', 'latex');
    apply_nature_style(fig_zoom, ax_zoom, closeup_width_mm, closeup_height_mm, panel_font, panel_fontsize, axis_linewidth);
    export_vector(fig_zoom, fullfile(export_dir, sprintf('capacitance_closeup_cycle%g_ch%d', ctr, T3_idx)));
end

% %% ===================== FIGURE 2 =====================
% % First 10 cycles, loading vs unloading, with exact crossing points marked
% first_n = min(10, num_cycles);
% cmap10 = lines(max(first_n, 1));
% 
% figure('Name', 'Figure 8: First 10 Cycles', 'Color', 'w', 'Position', [200, 200, 1000, 500]);
% 
% yyaxis left
% hold on;
% for k = 1:first_n
%     idx_range = trough_idx_sync(k):trough_idx_sync(k+1);
%     t_vals = time_sync(idx_range);
%     c_vals = cap_interp(idx_range) * 1000;
%     f_vals = force_interp(idx_range);
% 
%     [~, rel_peak] = max(f_vals);
%     loading   = 1:rel_peak;
%     unloading = rel_peak:numel(f_vals);
% 
%     plot(t_vals(loading), c_vals(loading), '-', 'Color', cmap10(k, :), 'LineWidth', 1.2);
%     plot(t_vals(unloading), c_vals(unloading), ':', 'Color', cmap10(k, :), 'LineWidth', 1.6);
% end
% plot(t_high_load_exact(1:first_n),   cap_high_load_exact(1:first_n),   '^', ...
%     'MarkerSize', 7, 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
% plot(t_high_unload_exact(1:first_n), cap_high_unload_exact(1:first_n), 'v', ...
%     'MarkerSize', 7, 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
% plot(t_low_load_exact(1:first_n),    cap_low_load_exact(1:first_n),    '^', ...
%     'MarkerSize', 7, 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
% plot(t_low_unload_exact(1:first_n),  cap_low_unload_exact(1:first_n),  'v', ...
%     'MarkerSize', 7, 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
% ylabel('Capacitance, C - C_0 [fF]');
% set(gca, 'YColor', 'k');
% 
% yyaxis right
% hold on;
% for k = 1:first_n
%     idx_range = trough_idx_sync(k):trough_idx_sync(k+1);
%     t_vals = time_sync(idx_range);
%     f_vals = force_interp(idx_range);
%     [~, rel_peak] = max(f_vals);
%     loading   = 1:rel_peak;
%     unloading = rel_peak:numel(f_vals);
%     plot(t_vals(loading), f_vals(loading), '--', 'Color', cmap10(k, :), 'LineWidth', 1.2);
%     plot(t_vals(unloading), f_vals(unloading), '-.', 'Color', cmap10(k, :), 'LineWidth', 1.6);
% end
% ylabel('Force [N]');
% set(gca, 'YColor', 'k');
% yline(target_high, ':', sprintf('%.2fN', target_high), 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
% yline(target_low, ':', sprintf('%.2fN', target_low), 'Color', [0.15 0.15 0.15], 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
% 
% for k = 1:first_n
%     idx_range = trough_idx_sync(k):trough_idx_sync(k+1);
%     t_vals = time_sync(idx_range);
%     if ~isempty(t_vals)
%         xline(t_vals(1), ':', sprintf('C%d', k), 'Color', [0.5 0.5 0.5], ...
%             'LabelOrientation', 'horizontal', 'FontSize', 7);
%     end
% end
% hold off;
% grid on;
% xlabel('Time [s]');
% title(sprintf('Channel %d: First %d Cycles - Loading (solid/dashed) vs Unloading (dotted/dash-dot)', T3_idx, first_n));
% 


% %% ===================== FIGURES 3 & 4 =====================
% % Cycle number vs. capacitance, loading-only and unloading-only
% cap_low_unload_shifted = nan(num_cycles, 1);
% cap_low_unload_shifted(1:end-1) = cap_low_unload_exact(2:end);
% 
% cap_ylim = [min([high_cap_exact; low_cap_exact], [], 'omitnan'), ...
%             max([high_cap_exact; low_cap_exact], [], 'omitnan')];
% 
% figure('Name', 'Figure 13: Cycle vs Capacitance (Loading)', 'Color', 'w', 'Position', [200, 200, 900, 450]);
% plot(cycle_num, cap_high_load_exact, '-', 'Color', orange, 'LineWidth', 0.75, 'DisplayName', 'High (~1N)');
% hold on;
% plot(cycle_num, cap_low_load_exact, '-', 'Color', blue, 'LineWidth', 0.75, 'DisplayName', 'Low (~0.02N)');
% hold off;
% grid on; ylim(cap_ylim);
% xlabel('Cycle Number'); ylabel('Capacitance, C - C_0 [fF]');
% legend('Location', 'best');
% title(sprintf('Channel %d: High/Low Capacitance vs. Cycle Number (Loading Only)', T3_idx));
% 
% figure('Name', 'Figure 14: Cycle vs Capacitance (Unloading)', 'Color', 'w', 'Position', [200, 200, 900, 450]);
% plot(cycle_num, cap_high_unload_exact, '-', 'Color', orange, 'LineWidth', 0.75, 'DisplayName', 'High (~1N)');
% hold on;
% plot(cycle_num, cap_low_unload_shifted, '-', 'Color', blue, 'LineWidth', 0.75, 'DisplayName', 'Low (~0.02N)');
% hold off;
% grid on; ylim(cap_ylim);
% xlabel('Cycle Number'); ylabel('Capacitance, C - C_0 [fF]');
% legend('Location', 'best');
% title(sprintf('Channel %d: High/Low Capacitance vs. Cycle Number (Unloading Only)', T3_idx));

%% ===================== FIGURE 5 =====================
% Peak-to-peak drift, loading-only, vs. baseline.
% Each plotted point is the AVERAGE of the p2p values over a SLIDING window
% of p2p_window cycles (moving average, stepping one cycle at a time),
% instead of a single cycle. Edit p2p_window to change the smoothing window.
p2p_cap_loading = cap_high_load_exact - cap_low_load_exact;

p2p_window = 10;   % <-- EDITABLE: sliding-window length (in cycles) for the moving average
baseline_n = 10;   % baseline = mean of the first baseline_n cycles (raw per-cycle p2p)

% --- Robustness guards so the values are always calculable ---
p2p_window = max(1, round(p2p_window));   % keep the window a positive whole number
bn = min(baseline_n, num_cycles);         % never index past the available cycles

% Sliding moving average of the per-cycle p2p. 'omitnan' skips missing
% cycles inside each window; movmean's default endpoint rule shrinks the
% window at the very start/end, so the first and last cycles still get a
% value even when a full window of p2p_window cycles isn't available yet.
p2p_loading_window = movmean(p2p_cap_loading, p2p_window, 'omitnan');

baseline_loading = mean(p2p_cap_loading(1:bn), 'omitnan');
drift_loading_window_pct = (p2p_loading_window - baseline_loading) / baseline_loading * 100;

% Consistency metric reported on the RAW per-cycle p2p (smoothing would
% artificially deflate the spread), so it reflects true cycle-to-cycle variation.
cv_loading_pct = std(p2p_cap_loading, 'omitnan') / mean(p2p_cap_loading, 'omitnan') * 100;

figure('Name', 'Figure 18: P2P Drift, Sliding-Window (Loading-Only)', 'Color', 'w', 'Position', [200, 200, 900, 450]);
plot(cycle_num, drift_loading_window_pct, '-', 'Color', blue, 'LineWidth', 1.5);
hold on;
yline(0, '--', 'Color', [0.4 0.4 0.4]);
hold off;
grid on;
xlabel('Cycle Number');
ylabel(sprintf('Capacitance P2P Drift [%% vs. mean of first %d cycles]', baseline_n));
legend({sprintf('%d-cycle moving average', p2p_window)}, 'Location', 'best');
title(sprintf('Channel %d: Loading-Only P2P Drift, %d-cycle window (per-cycle CV = %.2f%%)', T3_idx, p2p_window, cv_loading_pct));

% %% ===================== FIGURE 6 =====================
% % CV of peak-to-peak, computed over the growing window from cycle 1 up
% % through each comparison cycle (1000, 2000, 3000, ...) -- three rows:
% % average (loading+unloading combined), loading-only, unloading-only.
% 
% p2p_cap_average = high_cap_exact - low_cap_exact; % combined average, as in Figure 1's underlying data
% 
% p2p_cap_unloading = nan(num_cycles, 1);
% p2p_cap_unloading(1:end-1) = cap_high_unload_exact(1:end-1) - cap_low_unload_exact(2:end);
% 
% step_size_p2p = 1000;
% compare_cycles_p2p = step_size_p2p:step_size_p2p:num_cycles;
% if compare_cycles_p2p(end) < num_cycles
%     compare_cycles_p2p = [compare_cycles_p2p, num_cycles];
% end
% 
% cv_average_growing   = nan(size(compare_cycles_p2p));
% cv_loading_growing   = nan(size(compare_cycles_p2p));
% cv_unloading_growing = nan(size(compare_cycles_p2p));
% 
% for i = 1:numel(compare_cycles_p2p)
%     n_end = compare_cycles_p2p(i);
%     w_avg = p2p_cap_average(1:n_end);
%     w_load = p2p_cap_loading(1:n_end);
%     w_unload = p2p_cap_unloading(1:n_end);
%     cv_average_growing(i)   = std(w_avg, 'omitnan')    / mean(w_avg, 'omitnan')    * 100;
%     cv_loading_growing(i)   = std(w_load, 'omitnan')   / mean(w_load, 'omitnan')   * 100;
%     cv_unloading_growing(i) = std(w_unload, 'omitnan') / mean(w_unload, 'omitnan') * 100;
% end
% 
% figure('Name', 'Figure 7: P2P CV Growing Window (Average/Loading/Unloading)', 'Color', 'w', 'Position', [200, 200, 800, 750]);
% 
% subplot(3, 1, 1);
% plot(compare_cycles_p2p, cv_average_growing, '-o', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3]);
% ylim([1, 7]);
% grid on; ylabel('CV% (Average)');
% title(sprintf('Channel %d: P2P CV, Growing Window From Cycle 1', T3_idx));
% 
% subplot(3, 1, 2);
% plot(compare_cycles_p2p, cv_loading_growing, '-o', 'Color', blue, 'MarkerFaceColor', blue);
% ylim([1, 7]);
% grid on; ylabel('CV% (Loading)');
% 
% subplot(3, 1, 3);
% plot(compare_cycles_p2p, cv_unloading_growing, '-o', 'Color', orange, 'MarkerFaceColor', orange);
% ylim([1, 7]);
% grid on; ylabel('CV% (Unloading)');
% xlabel('Cycle Number (cumulative through this cycle)');
% 
% fprintf('P2P CV%% (average): starts %.2f%%, ends %.2f%%.\n', cv_average_growing(1), cv_average_growing(end));
% fprintf('P2P CV%% (loading): starts %.2f%%, ends %.2f%%.\n', cv_loading_growing(1), cv_loading_growing(end));
% fprintf('P2P CV%% (unloading): starts %.2f%%, ends %.2f%%.\n', cv_unloading_growing(1), cv_unloading_growing(end));
 

% %% ===================== FIGURE 7 =====================
% % CV of 1N capacitance, computed over the window from cycle 1 up through
% % each comparison cycle (1000, 2000, 3000, ...) -- i.e. "how variable has
% % the 1N reading been, from the very start up to this point."
% step_size_1N = 1000;
% 
% compare_cycles_1N = step_size_1N:step_size_1N:num_cycles;
% if compare_cycles_1N(end) < num_cycles
%     compare_cycles_1N = [compare_cycles_1N, num_cycles];
% end
% 
% cap_1N_cv_pct = nan(size(compare_cycles_1N));
% for i = 1:numel(compare_cycles_1N)
%     window_vals = high_cap_exact(1:compare_cycles_1N(i));
%     cap_1N_cv_pct(i) = std(window_vals, 'omitnan') / mean(window_vals, 'omitnan') * 100;
% end
% 
% figure('Name', 'Figure 21: 1N Capacitance Coefficient of Variation vs Cycle', 'Color', 'w', 'Position', [200, 200, 850, 450]);
% plot(compare_cycles_1N, cap_1N_cv_pct, '-o', 'Color', orange, 'MarkerFaceColor', orange, 'LineWidth', 1.2);
% grid on;
% ylim([0, 5]);
% xlabel('Cycle Number');
% ylabel('Coefficient of Variation % of 1N Capacitance (cycle 1 through this cycle)');
% title(sprintf('Channel %d: 1N Capacitance CV, Growing Window From Cycle 1', T3_idx));
% 
% fprintf('1N Loading Capacitance Coefficient of Variation%%: starts at %.2f%% (cycle 1-%d), ends at %.2f%% (cycle 1-%d).\n', ...
%     cap_1N_cv_pct(1), compare_cycles_1N(1), cap_1N_cv_pct(end), compare_cycles_1N(end));


%% Local Functions
function [t_cross, f_cross] = find_crossing(t_seg, f_seg, target, method)
% Finds the time t_cross at which a smooth interpolant of f_seg vs t_seg
% crosses "target". If no crossing exists (target never reached in this
% segment), falls back to the segment's own closest achievable sample.

    t_cross = NaN; f_cross = NaN;
    if numel(t_seg) < 2
        return
    end

    diffs = f_seg - target;
    sign_change_idx = find(diffs(1:end-1) .* diffs(2:end) < 0, 1, 'first');

    if isempty(sign_change_idx)
        [~, i_closest] = min(abs(diffs));
        t_cross = t_seg(i_closest);
        f_cross = f_seg(i_closest);
        return
    end

    interp_fun = @(tq) interp1(t_seg, f_seg, tq, method, 'extrap') - target;
    t_lo = t_seg(sign_change_idx);
    t_hi = t_seg(sign_change_idx + 1);
    t_cross = fzero(interp_fun, [t_lo, t_hi]);
    f_cross = target;
end

function apply_nature_style(figH, axH, w_mm, h_mm, fname, fsize, alw)
% Style an axes as a single Nature-style panel and set the figure to the
% requested physical size in mm. Sans-serif (Arial), 5-7 pt, thin axes.
    set(figH, 'Color', 'w', 'Units', 'centimeters');
    p = get(figH, 'Position');
    set(figH, 'Position', [p(1), p(2), w_mm/10, h_mm/10]);      % mm -> cm
    set(figH, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'auto');

    set(axH, 'FontName', fname, 'FontSize', fsize, 'LineWidth', alw, ...
        'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02 0.02], ...
        'Layer', 'top', 'Units', 'normalized', 'Position', [0.20 0.21 0.75 0.72]);

    set([get(axH,'XLabel'), get(axH,'YLabel'), get(axH,'Title')], ...
        'FontName', fname, 'FontSize', fsize);
end

function export_vector(figH, out_base)
% Export a figure to vector PDF (+ SVG) with editable text (no outlining).
% out_base is a full path WITHOUT extension.
    out_dir = fileparts(out_base);
    if ~isempty(out_dir) && ~isfolder(out_dir); mkdir(out_dir); end

    % Primary: tight vector PDF (Nature-preferred), text kept editable.
    try
        exportgraphics(figH, [out_base '.pdf'], 'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    catch
        print(figH, [out_base '.pdf'], '-dpdf', '-painters');   % older-MATLAB fallback
    end

    % Also an editable-layer SVG for figure assembly (Illustrator/Inkscape).
    try
        print(figH, [out_base '.svg'], '-dsvg', '-painters');
    catch
        % SVG optional; ignore if unsupported on this install.
    end

    fprintf('Exported: %s.pdf (+ .svg)\n', out_base);
end