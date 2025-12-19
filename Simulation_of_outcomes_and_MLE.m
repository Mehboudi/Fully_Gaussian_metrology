% ---------------------------------------------------------
% Parameter Settings
% ---------------------------------------------------------
alpha0 = sqrt(2);      % Tunable parameter alpha (related to input coherent state's amplitude)
nu = 10;                % Total number of samples per simulation
N_sims = 100000;          % Number of independent simulations
true_r = 0.01;            % The true hidden parameter r

% ---------------------------------------------------------
% Run Simulations (Left Column & Bottom Right)
% ---------------------------------------------------------
disp('Running Simulation for m=2...');
[est_m2, sq_err_m2, avg_err_m2] = run_simulation(2, nu, N_sims, alpha0, true_r);

disp('Running Simulation for m=5...');
[est_m5, sq_err_m5, avg_err_m5] = run_simulation(5, nu, N_sims, alpha0, true_r);

% ---------------------------------------------------------
% Prepare Theoretical Bounds
% ---------------------------------------------------------
sample_indices = (1:nu)';

% -- Bounds for m=2 --
m_val = 2;
CRB_m2       = 1 ./ (sample_indices * (2*(m_val-1) + 4 * m_val * alpha0^2));
QCRB_m2      = 1 ./ (sample_indices * (m_val * (4 * alpha0^2 + 2)));
Local_CRB_m2 = 1 ./ (sample_indices * (4 * m_val * alpha0^2));

% -- Bounds for m=5 --
m_val = 5;
CRB_m5       = 1 ./ (sample_indices * (2*(m_val-1) + 4 * m_val * alpha0^2));
QCRB_m5      = 1 ./ (sample_indices * (m_val * (4 * alpha0^2 + 2)));
Local_CRB_m5 = 1 ./ (sample_indices * (4 * m_val * alpha0^2));

% -- Fisher Information vs m (for Panel 3) --
m_range = 1:100;
FI_Global  = (2.*(m_range-1) + 4 .* m_range .* alpha0^2); 
FI_Local   = (4 .* m_range .* alpha0^2);
FI_Quantum = (m_range .* (4 * alpha0^2 + 2));

% ---------------------------------------------------------
% Plotting
% ---------------------------------------------------------
figure('Name', 'MLE Analysis: Convergence and Scaling', 'Color', 'w', 'Position', [100, 100, 1200, 900]);

% --- Panel 1: MLE Trajectory (m=2) ---
subplot(2,2,1);
plot(1:nu, est_m2(:,1), 'b-', 'LineWidth', 3); hold on;
yline(true_r, 'r--', 'LineWidth', 3); % Label removed from graph
title('MLE Estimate (Single Traj, m=2)');
%ylabel('Estimated r');
xlabel('$\nu$','Interpreter','latex');
xlim([1,nu]);
legend('Estimate', 'True r', 'Location', 'best');
grid on;
% Apply formatting: Font 14, Box on, Box LineWidth 1
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% --- Panel 2: Fisher Information vs m ---
subplot(2,2,2);
loglog(m_range, FI_Global./m_range, 'r--o', 'LineWidth', 3, 'MarkerSize', 6); hold on;
loglog(m_range, FI_Local./m_range, 'g-.s', 'LineWidth', 3, 'MarkerSize', 6);
loglog(m_range, FI_Quantum./m_range, 'b:^', 'LineWidth', 3, 'MarkerSize', 6);
xlabel('m');
%ylabel('Fisher Information (per sample)');
legend('${\cal F}_{\rm Global}^C/m$', '${\cal F}_{\rm Global}^C/m$', '${\cal F}^Q/m', 'Interpreter','Latex');
xlim([1, m_range(end)]);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% --- Panel 3: Convergence (m=2) ---
subplot(2,2,3);
% Plot gray trajectories first so they are in background
for i = 1:min(50, N_sims) 
    loglog(1:nu, sq_err_m2(:, i), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); hold on;
end
% Plot major curves with LineWidth 1
h_avg   = loglog(1:nu, avg_err_m2, 'k-', 'LineWidth', 3); 
h_crb   = loglog(1:nu, CRB_m2, 'r--', 'LineWidth', 3);
h_qcrb  = loglog(1:nu, QCRB_m2, 'b:', 'LineWidth', 3);
h_local = loglog(1:nu, Local_CRB_m2, 'g-.', 'LineWidth', 3);

xlabel('$\nu$','Interpreter','latex');
%ylabel('Squared Error');
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
%ylabel('Squared Error');
title('CRB (m=5)');
legend([h_avg, h_crb, h_qcrb, h_local], ...
       'MSE', 'CRB Global', 'QCRB', 'CRB Local', ...
       'Location', 'best');
ylim([5e-4, .1]); 
xlim([1, nu]);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on');

% ---------------------------------------------------------
% Helper Function: Simulation Runner
% ---------------------------------------------------------
function [all_estimates, squared_errors, avg_error] = run_simulation(m, nu, N_sims, alpha0, true_r)
    
    all_estimates = zeros(nu, N_sims);
    options = optimset('Display', 'off', 'TolX', 1e-6);
    
    var1 = exp(2 * true_r);
    var2 = (exp(2*true_r) + exp(-2*true_r))/2; 
    mean_d = sqrt(m) * alpha0 * (exp(true_r) - exp(-true_r)); 
    
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
end

% ---------------------------------------------------------
% Helper Function: Negative Log-Likelihood
% ---------------------------------------------------------
function nll = neg_log_likelihood(r, X_m, S1, nu_curr, m, alpha)
    e_2r = exp(2*r);
    e_neg_2r = exp(-2*r);
    denom = e_2r + e_neg_2r; 
    
    term_inner = sqrt(m) * alpha * (exp(r) - exp(-r)) - X_m;
    S2 = sum(term_inner.^2);
    
    L1 = - S1 / (2 * e_2r);
    L2 = - nu_curr * (m - 1) * r;
    L3 = - S2 / denom;
    L4 = - (nu_curr / 2) * log(denom);
    
    nll = -(L1 + L2 + L3 + L4);
end
