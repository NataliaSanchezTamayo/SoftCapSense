function results = calculate_T3_area_hysteresis(force_T3,cap_T3,T3_idx)
%CALCULATE_T3_AREA_HYSTERESIS Area-based hysteresis for the pushed taxel.
%
% Usage:
%   results = calculate_T3_area_hysteresis(force_T3,cap_T3,7);
%
% The loading and unloading branches are averaged independently in force
% bins, then interpolated onto the same force vector. Consequently,
% C_load_q(k) and C_unload_q(k) always correspond to exactly the same force.

    if nargin < 3
        T3_idx = 7;
    end

    if T3_idx ~= 7
        warning('This dataset was validated for the pressed taxel');
    end

    % Force maxima verified from the uploaded dataset.
    peaks_idx = [5218,11835,18425,24978,31525, ...
                 38056,44575,51083,57556,64030];

    half_window = 3200;
    force_bin_width = 0.02;      % N
    integration_step = 0.005;    % N
    min_samples_per_bin = 3;

    n_cycles = numel(peaks_idx);
    Harea_percent = nan(n_cycles,1);
    Aloop_pF_N = nan(n_cycles,1);
    Aload_pF_N = nan(n_cycles,1);
    Aunload_pF_N = nan(n_cycles,1);
    Fmin_N = nan(n_cycles,1);
    Fmax_N = nan(n_cycles,1);
    N_matched_force_points = zeros(n_cycles,1);

    for cycle_id = 1:n_cycles

        peak_idx = peaks_idx(cycle_id);

        idx_load = (peak_idx-half_window):peak_idx;
        idx_unload = peak_idx:(peak_idx+half_window);

        F_load = force_T3(idx_load);
        C_load = cap_T3(idx_load,T3_idx);
        F_unload = force_T3(idx_unload);
        C_unload = cap_T3(idx_unload,T3_idx);

        F_load = F_load(:);
        C_load = C_load(:);
        F_unload = F_unload(:);
        C_unload = C_unload(:);

        valid_load = isfinite_sample(F_load,C_load);
        valid_unload = isfinite_sample(F_unload,C_unload);

        F_load = F_load(valid_load);
        C_load = C_load(valid_load);
        F_unload = F_unload(valid_unload);
        C_unload = C_unload(valid_unload);

        % Construct identical force bins over the shared force range.
        common_min = max(min(F_load),min(F_unload));
        common_max = min(max(F_load),max(F_unload));
        edges = common_min:force_bin_width:(common_max+force_bin_width);

        [F_load_bin,C_load_bin] = bin_branch( ...
            F_load,C_load,edges,min_samples_per_bin);

        [F_unload_bin,C_unload_bin] = bin_branch( ...
            F_unload,C_unload,edges,min_samples_per_bin);

        % This shared vector enforces exact force alignment.
        Fq_min = max(min(F_load_bin),min(F_unload_bin));
        Fq_max = min(max(F_load_bin),max(F_unload_bin));
        Fq = (Fq_min:integration_step:Fq_max)';

        C_load_q = interp1(F_load_bin,C_load_bin,Fq,'linear');
        C_unload_q = interp1(F_unload_bin,C_unload_bin,Fq,'linear');

        % Cycle-specific zero-force baseline.
        baseline_values = C_load(F_load>=0 & F_load<=0.05);
        if isempty(baseline_values)
            C0 = C_load(1);
        else
            C0 = mean(baseline_values,'omitnan');
        end

        deltaC_load = C_load_q-C0;
        deltaC_unload = C_unload_q-C0;

        A_load = trapz(Fq,deltaC_load);
        A_unload = trapz(Fq,deltaC_unload);

        % Absolute pointwise separation is robust if the curves cross.
        A_loop = trapz(Fq,abs(C_unload_q-C_load_q));
        H_area = 100*A_loop/abs(A_load);

        Harea_percent(cycle_id) = H_area;
        Aloop_pF_N(cycle_id) = A_loop;
        Aload_pF_N(cycle_id) = A_load;
        Aunload_pF_N(cycle_id) = A_unload;
        Fmin_N(cycle_id) = Fq(1);
        Fmax_N(cycle_id) = Fq(end);
        N_matched_force_points(cycle_id) = numel(Fq);

        % Cycle 1 diagnostic: shaded area is integrated at matched forces.
        if cycle_id == 1
            figure('Name','Taxel 7 hysteresis alignment');
            fill([Fq;flipud(Fq)], ...
                 [C_load_q;flipud(C_unload_q)], ...
                 [0.75 0.82 0.90], ...
                 'FaceAlpha',0.35,'EdgeColor','none');
            hold on
            plot(Fq,C_load_q,'b-','LineWidth',1.5);
            plot(Fq,C_unload_q,'r--','LineWidth',1.5);
            xlabel('Force [N]');
            ylabel('C-C_0 [pF]');
            legend('Integrated area','Loading','Unloading', ...
                   'Location','best');
            title('Cycle 1: loading and unloading on identical force grid');
        end
    end

    Cycle = (1:n_cycles)';
    results = table(Cycle,Harea_percent,Aloop_pF_N,Aload_pF_N, ...
                    Aunload_pF_N,Fmin_N,Fmax_N, ...
                    N_matched_force_points);

    fprintf(['Taxel 7 area hysteresis = %.2f +/- %.2f %% ', ...
             '(mean +/- SD across %d cycles)\n'], ...
            mean(Harea_percent,'omitnan'), ...
            std(Harea_percent,'omitnan'),n_cycles);
end


function valid = isfinite_sample(F,C)
    valid = ~isnan(F) & ~isinf(F) & ~isnan(C) & ~isinf(C);
end


function [F_bin,C_bin] = bin_branch(F,C,edges,min_samples)

    bin_id = discretize(F,edges);
    n_bins = numel(edges)-1;

    valid = ~isnan(bin_id) & ~isinf(bin_id) & ...
            ~isnan(F) & ~isinf(F) & ~isnan(C) & ~isinf(C);

    samples_per_bin = accumarray( ...
        bin_id(valid),1,[n_bins,1],@sum,0);

    F_bin = accumarray( ...
        bin_id(valid),F(valid),[n_bins,1],@mean,NaN);

    C_bin = accumarray( ...
        bin_id(valid),C(valid),[n_bins,1],@mean,NaN);

    keep = samples_per_bin>=min_samples & ...
           ~isnan(F_bin) & ~isnan(C_bin);

    F_bin = F_bin(keep);
    C_bin = C_bin(keep);
end
