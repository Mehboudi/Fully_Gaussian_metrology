% ---------------------------------------------------------
% Parameter Settings
% ---------------------------------------------------------
alpha0 = sqrt(2);      % Tunable parameter alpha (related to input coherent state's amplitude)
nu = 100;              % Total number of samples per simulation
N_sims = 100000;       % Number of independent simulations
true_r = 0.1;          % The true hidden parameter r

% ---------------------------------------------------------
% Run Simulations (Left Column & Bottom Right)
% ---------------------------------------------------------
disp('Running Simulation for m=2...');
[est_m2, sq_err_m2, avg_err_m2, bias_m2, var_m2] = run_simulation(2, nu, N_sims, alpha0, true_r);

disp('Running Simulation for m=5...');
[est_m5, sq_err_m5, avg_err_m5, bias_m5, var_m5] = run_simulation(5, nu, N_sims, alpha0, true_r);

% ---------------------------------------------------------
% Numerical uncertainty of the squared bias (bootstrap)
% ---------------------------------------------------------
B_boot = 1000;   % number of bootstrap resamples
[sb_lo_m2, sb_hi_m2] = bootstrap_sqbias(est_m2, true_r, B_boot);
[sb_lo_m5, sb_hi_m5] = bootstrap_sqbias(est_m5, true_r, B_boot);

% ---------------------------------------------------------
% Prepare Theoretical Bounds
% ---------------------------------------------------------
sample_indices = (1:nu)';

% -- Bounds for m=2 --
m_val = 2;
CRB_m2       = 1 ./ (sample_indices .* (2*(m_val-1) + 4 * m_val * alpha0^2));
QCRB_m2      = 1 ./ (sample_indices .* (m_val * (4 * alpha0^2 + 2)));
Local_CRB_m2 = 1 ./ (sample_indices .* (4 * m_val * alpha0^2));

% -- Bounds for m=5 --
m_val = 5;
CRB_m5       = 1 ./ (sample_indices .* (2*(m_val-1) + 4 * m_val * alpha0^2));
QCRB_m5      = 1 ./ (sample_indices .* (m_val * (4 * alpha0^2 + 2)));
Local_CRB_m5 = 1 ./ (sample_indices .* (4 * m_val * alpha0^2));

% -- Fisher Information vs m (for Panel 2) --
m_range = 1:100;
FI_Global  = (2.*(m_range-1) + 4 .* m_range .* alpha0^2); 
FI_Local   = (4 .* m_range .* alpha0^2);
FI_Quantum = (m_range .* (4 * alpha0^2 + 2));

% ---------------------------------------------------------
% Plotting Figure 1: MLE Analysis: Convergence and Scaling
% ---------------------------------------------------------
figure('Name', 'MLE Analysis: Convergence and Scaling', 'Color', 'w', 'Position', [100, 100, 1200, 900]);

% --- Panel 1: MLE Trajectory (m=2) ---
subplot(2,2,1);
plot(1:nu, est_m2(:,1), 'b-', 'LineWidth', 3); hold on;
yline(true_r, 'r--', 'LineWidth', 3); 
title('MLE Estimate (Single Traj, m=2)');
xlabel('$\nu$','Interpreter','latex');
xlim([1,nu]);
legend('Estimate', 'True r', 'Location', 'best');
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% --- Panel 2: Fisher Information vs m ---
subplot(2,2,2);
loglog(m_range, FI_Global./m_range, 'r--o', 'LineWidth', 3, 'MarkerSize', 6); hold on;
loglog(m_range, FI_Local./m_range, 'g-.s', 'LineWidth', 3, 'MarkerSize', 6);
loglog(m_range, FI_Quantum./m_range, 'b:^', 'LineWidth', 3, 'MarkerSize', 6);
xlabel('m');
legend('${\cal F}_{\rm Global}^C/m$', '${\cal F}_{\rm Local}^C/m$', '${\cal F}^Q/m$', 'Interpreter','Latex');
xlim([1, m_range(end)]);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% --- Panel 3: Convergence (m=2) ---
subplot(2,2,3);
% Plot gray trajectories first so they are in background
for i = 1:min(50, N_sims) 
    loglog(1:nu, sq_err_m2(:, i), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); hold on;
end
% Plot major curves with LineWidth 3
h_avg   = loglog(1:nu, avg_err_m2, 'k-', 'LineWidth', 3); 
h_crb   = loglog(1:nu, CRB_m2, 'r--', 'LineWidth', 3);
h_qcrb  = loglog(1:nu, QCRB_m2, 'b:', 'LineWidth', 3);
h_local = loglog(1:nu, Local_CRB_m2, 'g-.', 'LineWidth', 3);

xlabel('$\nu$','Interpreter','latex');
title('CRB (m=2)');
legend([h_avg, h_crb, h_qcrb, h_local], ...
       'MSE', 'CRB Global', 'QCRB', 'CRB Local', ...
       'Location', 'best');
ylim([1e-3, .5]); 
xlim([1, nu]);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% --- Panel 4: Convergence (m=5) ---
subplot(2,2,4);
for i = 1:min(50, N_sims)
    loglog(1:nu, sq_err_m5(:, i), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); hold on;
end
h_avg   = loglog(1:nu, avg_err_m5, 'k-', 'LineWidth', 3);
h_crb   = loglog(1:nu, CRB_m5, 'r--', 'LineWidth', 3);
h_qcrb  = loglog(1:nu, QCRB_m5, 'b:', 'LineWidth', 3);
h_local = loglog(1:nu, Local_CRB_m5, 'g-.', 'LineWidth', 3);

xlabel('$\nu$','Interpreter','latex');
title('CRB (m=5)');
legend([h_avg, h_crb, h_qcrb, h_local], ...
       'MSE', 'CRB Global', 'QCRB', 'CRB Local', ...
       'Location', 'best');
ylim([5e-4, .1]); 
xlim([1, nu]);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% ---------------------------------------------------------
% Plotting Figure 2: MSE Components (Squared Bias vs Variance)
% ---------------------------------------------------------
figure('Name', 'MLE Analysis: MSE Breakdown', 'Color', 'w', 'Position', [150, 150, 800, 600]);

nu_ax   = (1:nu)';
y_floor = 1e-12;   % log axes cannot show 0: guards against a lower edge of exactly 0

% Log axes must be set before calling fill, otherwise fill makes them linear
ax = gca; set(ax, 'XScale', 'log', 'YScale', 'log'); hold(ax, 'on');

% Shaded bands first, so the curves are drawn on top
shade_band(nu_ax, max(sb_lo_m2, y_floor), sb_hi_m2, [0 0 1]);
shade_band(nu_ax, max(sb_lo_m5, y_floor), sb_hi_m5, [1 0 0]);

% Using loglog since squared errors span multiple orders of magnitude
h_b2 = loglog(1:nu, bias_m2.^2, 'b-', 'LineWidth', 3);
h_v2 = loglog(1:nu, var_m2, 'b--', 'LineWidth', 3);

h_b5 = loglog(1:nu, bias_m5.^2, 'r-', 'LineWidth', 3);
h_v5 = loglog(1:nu, var_m5, 'r--', 'LineWidth', 3);

title('MSE Breakdown: Squared Bias vs Estimator Variance');
xlabel('$\nu$','Interpreter','latex');
ylabel('Error');
legend([h_b2, h_v2, h_b5, h_v5], ...
       'm=2 (Squared Bias)', 'm=2 (Variance)', ...
       'm=5 (Squared Bias)', 'm=5 (Variance)', ...
       'Location', 'best');
xlim([1, nu]);
ylim([1e-9, 1e-1]);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% ---------------------------------------------------------
% Helper Function: Simulation Runner
% ---------------------------------------------------------
function [all_estimates, squared_errors, avg_error, bias, est_var] = run_simulation(m, nu, N_sims, alpha0, true_r)
    
    all_estimates = zeros(nu, N_sims);
    options = optimset('Display', 'off', 'TolX', 1e-6);
    
    var1 = exp(2 * true_r)/2;
    var2 = (exp(2*true_r) + exp(-2*true_r))/4; 
    mean_d = sqrt(m/2) * alpha0 * (exp(true_r) - exp(-true_r)); 
    
    h_wait = waitbar(0, ['Simulating m=' num2str(m) '...']);
    
    for sim_idx = 1:N_sims
        if mod(sim_idx, 50) == 0
            waitbar(sim_idx / N_sims, h_wait);
        end
        
        if m > 1
            X1 = sqrt(var1) * randn(nu, m - 1); 
        else
            X1 = zeros(nu, 0); 
        end
        X2 = mean_d + sqrt(var2) * randn(nu, 1);
        
        r0 = pi/10; 
        
        for s = 1:nu
            current_X1 = X1(1:s, :);
            current_X2 = X2(1:s);
            S1 = sum(current_X1(:).^2);
            
            cost_func = @(r) neg_log_likelihood(r, current_X2, S1, s, m, alpha0);
            
            if s > 1
                r0 = all_estimates(s-1, sim_idx); 
            end
            
            [r_est, ~] = fminsearch(cost_func, r0, options);
            all_estimates(s, sim_idx) = r_est;
        end
    end
    close(h_wait);
    
    squared_errors = (all_estimates - true_r).^2;
    avg_error = mean(squared_errors, 2);
    bias = mean(all_estimates, 2) - true_r; % Calculate the bias
    est_var = var(all_estimates, 0, 2);     % Calculate estimator variance over simulations
end

% ---------------------------------------------------------
% Helper Function: bootstrap band for the squared bias
% ---------------------------------------------------------
% Resample the N simulations (whole trajectories) with replacement, recompute
% the squared bias for every nu, repeat B times, and return the 2.5% and 97.5%
% percentiles of the B values at each nu.
function [lo, hi] = bootstrap_sqbias(all_estimates, true_r, B)
    [nu, N] = size(all_estimates);
    sqbias_boot = zeros(nu, B);
    for k = 1:B
        idx = randi(N, N, 1);                       % draw N simulations with replacement
        w   = accumarray(idx, 1, [N, 1]) / N;       % weight of each simulation in this resample
        sqbias_boot(:, k) = (all_estimates * w - true_r).^2;
    end
    sqbias_boot = sort(sqbias_boot, 2);
    lo = sqbias_boot(:, max(1, round(0.025 * B)));
    hi = sqbias_boot(:, round(0.975 * B));
end

% ---------------------------------------------------------
% Helper Function: shaded band between two curves
% ---------------------------------------------------------
function h = shade_band(x, lo, hi, col)
    x = x(:); lo = lo(:); hi = hi(:);
    h = fill([x; flipud(x)], [lo; flipud(hi)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% ---------------------------------------------------------
% Helper Function: Negative Log-Likelihood
% ---------------------------------------------------------
function nll = neg_log_likelihood(r, X_m, S1, nu_curr, m, alpha)
    term_inner = sqrt(m/2) * alpha * (exp(r) - exp(-r)) - X_m;
    S2 = sum(term_inner.^2);
    
    L1 = - S1 / exp(2*r);
    L2 = - S2 / cosh(2*r);
    L3 = - nu_curr * (m - 1) * r;
    L4 = - (nu_curr / 2) * log(cosh(2*r));
    
    nll = -(L1 + L2 + L3 + L4);
end
