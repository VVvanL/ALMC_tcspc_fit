function [results, final_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data, irf, x0, lb, ub, cost_type, fit_bg, fit_shift, error_type)
% FIT_TCSPC_VARPRO_BGOPTION_WERRORBARS 
%
% PURPOSE:
%   Fits Time-Correlated Single Photon Counting (TCSPC) data using the 
%   Variable Projection (VARPRO) algorithm. 
%
% MATHEMATICAL APPROACH:
%   - Non-linear parameters (lifetimes/taus and IRF shift) are optimized via fmincon.
%   - Linear parameters (amplitudes and background) are solved via linear 
%     least squares (lsqnonneg) inside the kernel for every fmincon iteration.
%   - This reduction in the optimizer's search space leads to faster convergence.
%
% FEATURES:
%   - fit_shift toggle: Removes the shift from the optimizer if false for speed.
%   - Sorting: Returns results with lifetimes in ascending order.
%   - Multi-mode Error Bars: Supports LS (Gaussian) and MLE (Poisson) statistics.

    % --- 1. Pre-processing and Setup ---
    t = t(:); data = data(:); irf = irf(:); 
    n_exp = length(x0) - 1;     % Number of decay components (last x0 element is shift)
    dt = t(2) - t(1);           % Time resolution
    T_pulse = t(end) + dt;      % Effective laser period (used for pile-up math)
    N = length(t);              % Total number of bins
    
    % --- 2. Fourier Domain Preparation ---
    % Pre-calculating the FFT of the IRF allows for extremely precise, 
    % sub-bin shifting using the Fourier Shift Theorem (Shift = multiply by phase).
    irf_f = fft(irf);           
    freq = (0:N-1)';            
    freq(freq > N/2) = freq(freq > N/2) - N; % Shift frequency axis for FFT centering
    freq = freq / T_pulse;      
    
    % --- 3. Optimizer Configuration ---
    % Using 'trust-region-reflective' which requires analytic gradients.
    % We provide these gradients (g) inside the varpro_kernel.
    options = optimoptions('fmincon', ...
        'SpecifyObjectiveGradient', true, ... 
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'off', ... 
        'OptimalityTolerance', 1e-7);

    % --- 4. Parameter Branching: Dynamic Search Space ---
    % To maximize speed when the shift is known, we physically remove it 
    % from the vectors passed to fmincon if fit_shift is false.
    if fit_shift
        % The optimizer handles [tau1, tau2... , shift]
        active_x0 = x0(:)';
        active_lb = lb(:)';
        active_ub = ub(:)';
        % Pass is_shift_active=true so kernel extracts shift from the parameter vector
        obj_fun = @(x) varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg, true, 0);
    else
        % Optimizer handles ONLY [tau1, tau2...]. Search space is reduced by 1 dimension.
        active_x0 = x0(1:n_exp)';
        active_lb = lb(1:n_exp)';
        active_ub = ub(1:n_exp)';
        % Pass is_shift_active=false and provide the fixed shift value as a constant
        fixed_shift_val = x0(end);
        obj_fun = @(x) varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg, false, fixed_shift_val);
    end

    % --- 5. Non-linear Optimization ---
    % fmincon calls the kernel repeatedly to find the optimal lifetimes (and shift).
    [best_nonlin_active, ~] = fmincon(obj_fun, active_x0, [], [], [], [], active_lb, active_ub, [], options);
    
    % --- 6. Final Reconstruction ---
    % Retrieve final amplitudes, model curve, and Jacobian matrix for error estimation.
    [~, ~, lin_params, final_fit, Jacob] = obj_fun(best_nonlin_active);
    
    % --- 7. Error Estimation (Covariance Matrix) ---
    % Calculate degrees of freedom (N - total free parameters).
    p_total = length(best_nonlin_active) + length(lin_params);
    dof = N - p_total; 
    
    if strcmpi(cost_type, 'LS')
        % Least Squares error based on Mean Squared Error (MSE)
        mse = sum((data - final_fit).^2) / dof;
        cov_mat = mse * ( (Jacob' * Jacob) \ eye(size(Jacob,2)) );
    else
        % Poisson MLE error based on the Fisher Information Matrix (Hessian approximation)
        W = 1 ./ max(final_fit, 1e-6); 
        cov_mat = (Jacob' * (W .* Jacob)) \ eye(size(Jacob,2));
    end
    
    % Calculate standard errors (sqrt of diagonal)
    err_vec = sqrt(diag(cov_mat));
    
    % Scaling for 95% Confidence Intervals if requested
    if strcmpi(error_type, '95CI')
        tcrit = 1.96; 
        if exist('tinv', 'file'), tcrit = tinv(0.975, dof); end
        err_vec = err_vec * tcrit;
    end
    
    % --- 8. Sorting and Packaging Results ---
    % To ensure consistent output, we sort by lifetime (taus) in ascending order.
    % We must ensure amplitudes and errors are re-indexed to match.
    unsorted_taus = best_nonlin_active(1:n_exp);
    unsorted_amps = lin_params(1:n_exp);
    unsorted_tau_errs = err_vec(1:n_exp);
    
    % Generate sorting index based on taus
    [sorted_taus, sort_idx] = sort(unsorted_taus);
    
    results.taus = sorted_taus;
    results.amplitudes = unsorted_amps(sort_idx);
    results.err_vals.taus = unsorted_tau_errs(sort_idx);
    
    % Handle Background result
    results.background = 0; 
    if fit_bg, results.background = lin_params(end); end
    
    % Map remaining errors (Shift and Amplitudes) depending on whether shift was fitted
    if fit_shift
        results.shift = best_nonlin_active(end);
        results.err_vals.shift = err_vec(n_exp + 1);
        % Amplitude errors start at column n_exp + 2 in the Jacobian if shift is fitted
        amp_errs_raw = err_vec(n_exp + 2 : 2*n_exp + 1);
        results.err_vals.amplitudes = amp_errs_raw(sort_idx);
    else
        results.shift = fixed_shift_val;
        results.err_vals.shift = 0; % No error for a fixed constant
        % Amplitude errors start at column n_exp + 1 if shift is NOT fitted
        amp_errs_raw = err_vec(n_exp + 1 : 2*n_exp);
        results.err_vals.amplitudes = amp_errs_raw(sort_idx);
    end
    
    if fit_bg, results.err_vals.background = err_vec(end); end
    results.error_type_used = error_type;
end

function [f, g, lin_params, model, Jacob] = varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg, is_shift_active, fixed_shift)
    % INTERNAL FUNCTION: The engine of the fit.
    % Calculates objective function value (f) and analytic gradient (g).
    
    taus = x(1:n_exp); 
    % Select the shift based on whether the optimizer is controlling it
    if is_shift_active, shift = x(end); else, shift = fixed_shift; end
    
    N = length(t); dt = t(2)-t(1);
    
    % --- Step A: IRF Shifting ---
    use_linear = (N < 250);
    if shift == 0
        irf_s = irf;
    else
        if use_linear
            % Use time-domain linear interpolation for small datasets
            irf_s = interp1(t, irf, t - shift, 'linear', 0);
        else
            % Use FFT Shift Theorem: Multiply by complex phase in frequency domain
            shift_kernel = exp(-2j * pi * freq * shift);
            irf_s = real(ifft(irf_f .* shift_kernel));
        end
    end

    % --- Step B: Basis Matrix (Phi) Construction ---
    % Phi contains the convolved decays. VARPRO assumes Model = Phi * Beta.
    n_lin = n_exp + double(fit_bg); 
    Phi = zeros(N, n_lin);          
    dPhi_dtau = zeros(N, n_exp); % Derivatives w.r.t non-linear lifetimes
    
    for i = 1:n_exp
        % Decay with Pulse Pile-up correction (Periodic Boundary Condition)
        decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
        
        if use_linear
            tmp = conv(decay, irf_s);
            Phi(:,i) = tmp(1:N) * dt; 
            % Analytic derivative of decay w.r.t tau
            d_dec = (t/taus(i)^2) .* decay;
            tmp_d = conv(d_dec, irf_s);
            dPhi_dtau(:,i) = tmp_d(1:N) * dt;
        else
            % Fourier domain convolution is multiplication: FFT(A)*FFT(B)
            decay_f = fft(decay);
            % Reuse current shifted IRF in frequency domain
            curr_irf_f = irf_f .* exp(-2j * pi * freq * shift); 
            Phi(:,i) = real(ifft(decay_f .* curr_irf_f));
            
            d_dec = (t/taus(i)^2) .* decay;
            dPhi_dtau(:,i) = real(ifft(fft(d_dec) .* curr_irf_f));
        end
    end
    if fit_bg, Phi(:, end) = 1; end % Ones for constant background
    
    % --- Step C: Linear Solve (Projecting out Amplitudes) ---
    if strcmpi(cost_type, 'LS')
        % Non-negative Least Squares solves the linear amplitudes
        lin_params = lsqnonneg(Phi, data);
        model = Phi * lin_params;
        res = model - data; 
        f = sum(res.^2); % Cost function
        w = 2 * res;     % Weights for gradient
    else
        % Poisson Maximum Likelihood via Iteratively Reweighted Least Squares
        beta = lsqnonneg(Phi, data); 
        for iter = 1:5 
            model = max(Phi * beta, 1e-6);
            W_vec = 1 ./ model;
            beta = (Phi' * (W_vec .* Phi)) \ (Phi' * (W_vec .* data));
            beta(beta < 0) = 0; 
        end
        lin_params = beta; model = max(Phi * lin_params, 1e-6);
        f = sum(model - data .* log(model)); % Poisson Likelihood Cost
        w = 1 - data./model; % Weights for gradient
    end

    % --- Step D: Construct Full Jacobian and Gradient ---
    if is_shift_active
        % If fitting shift, we need the derivative of the model w.r.t. shift
        if use_linear
            dirf_ds = gradient(irf_s) / dt;
            dM_ds = 0;
            for i = 1:n_exp
                 decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
                 tmp_s = conv(decay, dirf_ds);
                 dM_ds = dM_ds - lin_params(i) * tmp_s(1:N) * dt;
            end
        else
            % Derivative of the shift kernel is (-2j*pi*freq)*exp(-2j*pi*freq*shift)
            dirf_ds_f = irf_f .* (-2j*pi*freq) .* exp(-2j * pi * freq * shift);
            dM_ds = 0;
            for i = 1:n_exp
                 decay_f = fft(exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i))));
                 dM_ds = dM_ds + lin_params(i) * real(ifft(decay_f .* dirf_ds_f));
            end
        end
        
        % Jacobian Layout: [dTaus (n_exp), dShift (1), dLinear (n_lin)]
        Jacob = zeros(N, n_exp + 1 + n_lin);
        for i = 1:n_exp, Jacob(:,i) = lin_params(i) * dPhi_dtau(:,i); end
        Jacob(:, n_exp+1) = dM_ds; 
        Jacob(:, n_exp+2:end) = Phi; 
        
        if nargout > 1
            g = zeros(n_exp + 1, 1);
            for i = 1:n_exp, g(i) = sum(w .* Jacob(:, i)); end
            g(end) = sum(w .* dM_ds);
        end
    else
        % If NOT fitting shift, Jacobian is smaller: [dTaus (n_exp), dLinear (n_lin)]
        Jacob = zeros(N, n_exp + n_lin);
        for i = 1:n_exp, Jacob(:,i) = lin_params(i) * dPhi_dtau(:,i); end
        Jacob(:, n_exp+1:end) = Phi; 
        
        if nargout > 1
            g = zeros(n_exp, 1);
            for i = 1:n_exp, g(i) = sum(w .* Jacob(:, i)); end
        end
    end
end

% function [results, final_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data, irf, x0, lb, ub, cost_type, fit_bg, fit_shift, error_type)
% % FIT_TCSPC_VARPRO_BGOPTION_WERRORBARS 
% %
% % PURPOSE:
% %   Fits Time-Correlated Single Photon Counting (TCSPC) data using the 
% %   Variable Projection (VARPRO) algorithm. 
% %
% % INPUT OPTIONS & SETTINGS:
% %   x0        - Initial guess [tau1, tau2, ..., shift]. (length = n_exp + 1).
% %   lb, ub    - Lower and upper bounds for the non-linear parameters.
% %   cost_type - 'LS' (Least Squares) or 'MLE' (Poisson Maximum Likelihood).
% %   fit_bg    - true: include constant background; false: fix to zero.
% %   fit_shift - true: optimize IRF shift; false: fix to x0(end).
% %   error_type- '1sigma' or '95CI'.
% 
%     % --- 1. Pre-processing and Setup ---
%     t = t(:); data = data(:); irf = irf(:); 
%     n_exp = length(x0) - 1;     % Components determined by x0 length
%     dt = t(2) - t(1);           % Time resolution
%     T_pulse = t(end) + dt;      % Laser repetition period for wrap-around math
%     N = length(t);              
% 
%     % --- 2. Fourier Domain Preparation ---
%     % Pre-calculating the FFT of the IRF allows us to use the Fourier Shift Theorem.
%     % Shifting in frequency space is more precise than time-domain interpolation.
%     irf_f = fft(irf);           
%     freq = (0:N-1)';            
%     freq(freq > N/2) = freq(freq > N/2) - N; 
%     freq = freq / T_pulse;      
% 
%     % --- 3. Optimizer Configuration ---
%     % 'trust-region-reflective' requires analytic gradients (which we provide)
%     % but it crashes if LB == UB.
%     options = optimoptions('fmincon', ...
%         'SpecifyObjectiveGradient', true, ... 
%         'Algorithm', 'trust-region-reflective', ...
%         'Display', 'off', ... 
%         'OptimalityTolerance', 1e-7);
% 
%     % --- 4. Parameter Locking Logic (Fix for sfminbx error) ---
%     % To "fix" a parameter in trust-region-reflective, we can't set LB = UB.
%     % Instead, we set UB to be LB + a tiny epsilon (1e-12). This locks the 
%     % parameter effectively while satisfying the algorithm's requirements.
%     local_lb = lb(:)'; 
%     local_ub = ub(:)';
%     if ~fit_shift
%         local_lb(end) = x0(end);
%         local_ub(end) = x0(end) + 1e-12; % Tiny epsilon to prevent LB==UB error
%     end
% 
%     % Create handle for the VARPRO kernel
%     obj_fun = @(x) varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- 5. Non-linear Optimization ---
%     % fmincon searches for the best Lifetimes and Shift.
%     [best_nonlin, ~] = fmincon(obj_fun, x0(:)', [], [], [], [], local_lb, local_ub, [], options);
% 
%     % --- 6. Final Reconstruction ---
%     % Extract the linear parameters (Amplitudes/BG) and Jacobian using the best found values.
%     [~, ~, lin_params, final_fit, Jacob] = varpro_kernel(best_nonlin, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- 7. Error Estimation ---
%     % Calculate effective degrees of freedom. 
%     % Note: even if fit_shift is false, the Jacob matrix still contains the shift column.
%     n_lin = n_exp + double(fit_bg); 
%     p_total = n_exp + n_lin + double(fit_shift); 
%     dof = N - p_total; 
% 
%     if strcmpi(cost_type, 'LS')
%         mse = sum((data - final_fit).^2) / dof;
%         cov_mat = mse * ( (Jacob' * Jacob) \ eye(size(Jacob,2)) );
%     else
%         W = 1 ./ max(final_fit, 1e-6); 
%         cov_mat = (Jacob' * (W .* Jacob)) \ eye(size(Jacob,2));
%     end
% 
%     errs = sqrt(diag(cov_mat));
%     if strcmpi(error_type, '95CI')
%         if exist('tinv', 'file')
%             tcrit = tinv(0.975, dof); 
%         else
%             tcrit = 1.96; 
%         end
%         errs = errs * tcrit;
%     end
% 
%     % --- 8. Package Results ---
%     results.taus = best_nonlin(1:n_exp);
%     results.shift = best_nonlin(end);
%     results.amplitudes = lin_params(1:n_exp);
%     results.background = 0;
%     if fit_bg, results.background = lin_params(end); end
% 
%     results.err_vals.taus = errs(1:n_exp);
%     if fit_shift, results.err_vals.shift = errs(n_exp + 1); else, results.err_vals.shift = 0; end
%     results.err_vals.amplitudes = errs(n_exp + 2 : 2*n_exp + 1);
%     if fit_bg, results.err_vals.background = errs(end); end
%     results.error_type_used = error_type;
% end
% 
% function [f, g, lin_params, model, Jacob] = varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg)
%     % INTERNAL KERNEL: Computes Model, Cost (f), and Gradient (g).
% 
%     taus = x(1:n_exp); 
%     shift = x(end);
%     N = length(t);
%     dt = t(2) - t(1);
% 
%     % --- IRF Shifting ---
%     use_linear = (N < 250);
%     if use_linear
%         % Linear interpolation for short vectors
%         irf_s = interp1(t, irf, t - shift, 'linear', 0);
%     else
%         % Fourier Shift Theorem: Shifting by multiplying by a complex exponential phase.
%         shift_kernel = exp(-2j * pi * freq * shift);
%         irf_s_f = irf_f .* shift_kernel; 
%         irf_s = real(ifft(irf_s_f));
%     end
% 
%     % --- Basis Matrix (Phi) ---
%     % Phi contains the convolved decays. VARPRO assumes: Model = Phi * Amplitudes
%     n_lin = n_exp + double(fit_bg); 
%     Phi = zeros(N, n_lin);          
%     dPhi_dtau = zeros(N, n_exp);    
% 
%     for i = 1:n_exp
%         % Decay formula includes pulse pile-up correction for periodic excitation.
%         decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
% 
%         if use_linear
%             tmp = conv(decay, irf_s);
%             Phi(:,i) = tmp(1:N) * dt; 
%             d_dec = (t/taus(i)^2) .* decay; 
%             tmp_d = conv(d_dec, irf_s);
%             dPhi_dtau(:,i) = tmp_d(1:N) * dt;
%         else
%             % Multiplication in Fourier space = Convolution in Time space.
%             decay_f = fft(decay);
%             Phi(:,i) = real(ifft(decay_f .* irf_s_f));
%             d_dec = (t/taus(i)^2) .* decay;
%             dPhi_dtau(:,i) = real(ifft(fft(d_dec) .* irf_s_f));
%         end
%     end
%     if fit_bg, Phi(:, end) = 1; end 
% 
%     % --- Solve for Linear Parameters ---
%     if strcmpi(cost_type, 'LS')
%         % Solve for non-negative amplitudes (physical constraint)
%         lin_params = lsqnonneg(Phi, data);
%         model = Phi * lin_params;
%         res = model - data;
%         f = sum(res.^2);
%         w = 2 * res; 
%     else
%         % Poisson MLE via 5 iterations of IRLS (Iterative Reweighted Least Squares)
%         beta = lsqnonneg(Phi, data); 
%         for iter = 1:5 
%             model = Phi * beta;
%             model(model < 1e-6) = 1e-6; 
%             W_vec = 1 ./ model;
%             beta = (Phi' * (W_vec .* Phi)) \ (Phi' * (W_vec .* data));
%             beta(beta < 0) = 0; 
%         end
%         lin_params = beta;
%         model = Phi * lin_params;
%         model(model < 1e-6) = 1e-6;
%         f = sum(model - data .* log(model)); % Poisson Likelihood
%         w = 1 - data./model;
%     end
% 
%     % --- Gradient of Shift ---
%     % We compute the analytic derivative of the IRF w.r.t the shift parameter.
%     if use_linear
%         dirf_ds = gradient(irf_s) / dt;
%         dM_ds = 0;
%         for i = 1:n_exp
%              decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
%              tmp_s = conv(decay, dirf_ds);
%              dM_ds = dM_ds - lin_params(i) * tmp_s(1:N) * dt;
%         end
%     else
%         % In Fourier space, the derivative of the shift is just -j*2*pi*f
%         dirf_ds_f = irf_f .* (-2j*pi*freq) .* shift_kernel;
%         dM_ds = 0;
%         for i = 1:n_exp
%              decay_f = fft(exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i))));
%              dM_ds = dM_ds + lin_params(i) * real(ifft(decay_f .* dirf_ds_f));
%         end
%     end
% 
%     % --- Full Jacobian for Error propagation ---
%     Jacob = zeros(N, length(x) + n_lin);
%     for i = 1:n_exp
%         Jacob(:, i) = lin_params(i) * dPhi_dtau(:,i); 
%     end
%     Jacob(:, n_exp + 1) = dM_ds; 
%     Jacob(:, n_exp + 2 : end) = Phi; 
% 
%     % --- Gradient Output ---
%     if nargout > 1
%         g = zeros(size(x));
%         for i = 1:n_exp
%             g(i) = sum(w .* Jacob(:, i));
%         end
%         g(end) = sum(w .* dM_ds);
%     end
% end

% function [results, final_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data, irf, x0, lb, ub, cost_type, fit_bg, error_type)
% % FIT_TCSPC_VARPRO_BGOPTION_WERRORBARS 
% %
% % PURPOSE:
% %   Fits Time-Correlated Single Photon Counting (TCSPC) data using the 
% %   Variable Projection (VARPRO) algorithm. This method separates the 
% %   non-linear parameters (lifetimes, shift) from linear parameters 
% %   (amplitudes, background), reducing the search space for fmincon.
% %
% % INPUT DATA REQUIREMENTS:
% %   t         - [N x 1] Time axis vector (typically ns).
% %   data      - [N x 1] Photon counts (raw data).
% %   irf       - [N x 1] Instrument Response Function.
% %
% % INPUT OPTIONS & SETTINGS:
% %   x0        - Initial guess [tau1, tau2, ..., shift]. (length = n_exp + 1).
% %   lb, ub    - Lower and upper bounds for the non-linear parameters.
% %   cost_type - 'LS'   : Least Squares (Gaussian noise assumption).
% %               'MLE'  : Maximum Likelihood (Poisson noise assumption).
% %   fit_bg    - true   : Includes a constant background parameter in the fit.
% %               false  : Fixes background to zero.
% %   error_type- '1sigma': Standard error (approx. 68% confidence).
% %               '95CI'  : 95% Confidence Interval using Student's t-dist.
% 
%     % --- 1. Pre-processing ---
%     t = t(:); data = data(:); irf = irf(:); 
%     n_exp = length(x0) - 1;     % Number of exponential components
%     dt = t(2) - t(1);           % Time step
%     T_pulse = t(end) + dt;      % Effective pulse period for pile-up correction
%     N = length(t);              % Number of bins
% 
%     % --- 2. Speed Optimization: Fourier Domain Setup ---
%     irf_f = fft(irf);           % FFT of IRF for fast convolution
%     freq = (0:N-1)';            % Frequency axis for Fourier Shift Theorem
%     freq(freq > N/2) = freq(freq > N/2) - N;
%     freq = freq / T_pulse;
% 
%     % --- 3. Optimizer Configuration ---
%     options = optimoptions('fmincon', ...
%         'SpecifyObjectiveGradient', true, ... % Uses analytic gradients from kernel
%         'Algorithm', 'trust-region-reflective', ...
%         'Display', 'off', ... 
%         'OptimalityTolerance', 1e-7);
% 
%     % Create objective function handle (passes x to the kernel)
%     obj_fun = @(x) varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- 4. Non-linear Optimization ---
%     % Optimizes only [taus, shift]
%     [best_nonlin, ~] = fmincon(obj_fun, x0(:)', [], [], [], [], lb(:)', ub(:)', [], options);
% 
%     % --- 5. Final Reconstruction ---
%     % Retrieve optimal linear parameters (amps/bg), model curve, and Jacobian
%     [~, ~, lin_params, final_fit, Jacob] = varpro_kernel(best_nonlin, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- 6. Error Estimation (Covariance Matrix) ---
%     p_total = length(best_nonlin) + length(lin_params);
%     dof = N - p_total; % Degrees of freedom
% 
%     if strcmpi(cost_type, 'LS')
%         % Residual-based error for Least Squares
%         mse = sum((data - final_fit).^2) / dof;
%         cov_mat = mse * ( (Jacob' * Jacob) \ eye(p_total) );
%     else
%         % Information Matrix-based error for Poisson MLE
%         W = 1 ./ max(final_fit, 1e-6); % Weight matrix (inverse of model variance)
%         cov_mat = (Jacob' * (W .* Jacob)) \ eye(p_total);
%     end
% 
%     % Standard Errors (square root of the diagonal of covariance matrix)
%     errs = sqrt(diag(cov_mat));
% 
%     % Apply scaling for 95% Confidence Intervals
%     if strcmpi(error_type, '95CI')
%         if exist('tinv', 'file')
%             tcrit = tinv(0.975, dof); % Exact Student's t critical value
%         else
%             tcrit = 1.96; % Normal distribution approximation fallback
%         end
%         errs = errs * tcrit;
%     end
% 
%     % --- 7. Package Results ---
%     results.taus = best_nonlin(1:n_exp);
%     results.shift = best_nonlin(end);
%     results.amplitudes = lin_params(1:n_exp);
%     results.background = 0;
%     if fit_bg, results.background = lin_params(end); end
% 
%     % Store requested error values
%     results.err_vals.taus = errs(1:n_exp);
%     results.err_vals.shift = errs(n_exp + 1);
%     results.err_vals.amplitudes = errs(n_exp + 2 : 2*n_exp + 1);
%     if fit_bg, results.err_vals.background = errs(end); end
%     results.error_type_used = error_type;
% end
% 
% function [f, g, lin_params, model, Jacob] = varpro_kernel(x, t, data, irf, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg)
%     % Internal Function: Calculates Cost, Gradient, and Jacobian
% 
%     taus = x(1:n_exp); shift = x(end);
%     N = length(t);
%     dt = t(2) - t(1);
% 
%     % --- Determine Convolution Strategy ---
%     % Use linear convolution for small datasets to avoid circular artifacts (waviness)
%     use_linear = (N < 250);
% 
%     if use_linear
%         % Time-domain shifting (interpolation)
%         irf_s = interp1(t, irf, t - shift, 'linear', 0);
%     else
%         % Fourier-domain shifting
%         shift_kernel = exp(-2j * pi * freq * shift);
%         irf_s_f = irf_f .* shift_kernel; 
%         irf_s = real(ifft(irf_s_f));
%     end
% 
%     n_lin = n_exp + double(fit_bg); % Columns for amplitudes + background
%     Phi = zeros(N, n_lin);          % Basis matrix
%     dPhi_dtau = zeros(N, n_exp);    % Basis derivatives for gradient
% 
%     % Construct Basis Functions (Convolved Decays)
%     for i = 1:n_exp
%         % Decay with Periodic Boundary Condition (Pulse Pile-up correction)
%         decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
% 
%         if use_linear
%             tmp = conv(decay, irf_s);
%             Phi(:,i) = tmp(1:N) * dt; % Scale by dt for proper integration
% 
%             d_dec = (t/taus(i)^2) .* decay;
%             tmp_d = conv(d_dec, irf_s);
%             dPhi_dtau(:,i) = tmp_d(1:N) * dt;
%         else
%             decay_f = fft(decay);
%             Phi(:,i) = real(ifft(decay_f .* irf_s_f));
% 
%             d_dec = (t/taus(i)^2) .* decay;
%             dPhi_dtau(:,i) = real(ifft(fft(d_dec) .* irf_s_f));
%         end
%     end
%     if fit_bg, Phi(:, end) = 1; end % Add background column (ones)
% 
%     % Solve for Linear Parameters (Amplitudes/Background)
%     if strcmpi(cost_type, 'LS')
%         lin_params = lsqnonneg(Phi, data);
%         model = Phi * lin_params;
%         res = model - data;
%         f = sum(res.^2);
%         w = 2 * res;
%     else
%         beta = lsqnonneg(Phi, data); 
%         for iter = 1:5 
%             model = Phi * beta;
%             model(model < 1e-6) = 1e-6; 
%             W = 1 ./ model;
%             beta = (Phi' * (W .* Phi)) \ (Phi' * (W .* data));
%             beta(beta < 0) = 0; 
%         end
%         lin_params = beta;
%         model = Phi * lin_params;
%         model(model < 1e-6) = 1e-6;
%         f = sum(model - data .* log(model));
%         w = 1 - data./model;
%     end
% 
%     % Calculate Derivative w.r.t IRF Shift
%     if use_linear
%         % Numeric derivative for shift in time-domain
%         dirf_ds = gradient(irf_s) / dt;
%         dM_ds = 0;
%         for i = 1:n_exp
%              decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
%              tmp_s = conv(decay, dirf_ds);
%              dM_ds = dM_ds - lin_params(i) * tmp_s(1:N) * dt;
%         end
%     else
%         dirf_ds_f = irf_f .* (-2j*pi*freq) .* shift_kernel;
%         dM_ds = 0;
%         for i = 1:n_exp
%              decay_f = fft(exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i))));
%              dM_ds = dM_ds + lin_params(i) * real(ifft(decay_f .* dirf_ds_f));
%         end
%     end
% 
%     % Construct Full Jacobian Matrix
%     Jacob = zeros(N, length(x) + n_lin);
%     for i = 1:n_exp
%         Jacob(:, i) = lin_params(i) * dPhi_dtau(:,i); 
%     end
%     Jacob(:, n_exp + 1) = dM_ds; 
%     Jacob(:, n_exp + 2 : end) = Phi; 
% 
%     if nargout > 1
%         g = zeros(size(x));
%         for i = 1:n_exp
%             g(i) = sum(w .* Jacob(:, i));
%         end
%         g(end) = sum(w .* dM_ds);
%     end
% end


% function [results, final_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data, irf, x0, lb, ub, cost_type, fit_bg, error_type)
% % FIT_TCSPC_VARPRO_BGOPTION_WERRORBARS 
% %
% % PURPOSE:
% %   Fits Time-Correlated Single Photon Counting (TCSPC) data using the 
% %   Variable Projection (VARPRO) algorithm. This method separates the 
% %   non-linear parameters (lifetimes, shift) from linear parameters 
% %   (amplitudes, background), reducing the search space for fmincon.
% %
% % INPUT DATA REQUIREMENTS:
% %   t         - [N x 1] Time axis vector (typically ns).
% %   data      - [N x 1] Photon counts (raw data).
% %   irf       - [N x 1] Instrument Response Function.
% %
% % INPUT OPTIONS & SETTINGS:
% %   x0        - Initial guess [tau1, tau2, ..., shift]. (length = n_exp + 1).
% %   lb, ub    - Lower and upper bounds for the non-linear parameters.
% %   cost_type - 'LS'   : Least Squares (Gaussian noise assumption).
% %               'MLE'  : Maximum Likelihood (Poisson noise assumption).
% %   fit_bg    - true   : Includes a constant background parameter in the fit.
% %               false  : Fixes background to zero.
% %   error_type- '1sigma': Standard error (approx. 68% confidence).
% %               '95CI'  : 95% Confidence Interval using Student's t-dist.
% 
%     % --- 1. Pre-processing ---
%     t = t(:); data = data(:); irf = irf(:); 
%     n_exp = length(x0) - 1;     % Number of exponential components
%     dt = t(2) - t(1);           % Time step
%     T_pulse = t(end) + dt;      % Effective pulse period for pile-up correction
%     N = length(t);              % Number of bins
% 
%     % --- 2. Speed Optimization: Fourier Domain Setup ---
%     irf_f = fft(irf);           % FFT of IRF for fast convolution
%     freq = (0:N-1)';            % Frequency axis for Fourier Shift Theorem
%     freq(freq > N/2) = freq(freq > N/2) - N;
%     freq = freq / T_pulse;
% 
%     % --- 3. Optimizer Configuration ---
%     options = optimoptions('fmincon', ...
%         'SpecifyObjectiveGradient', true, ... % Uses analytic gradients from kernel
%         'Algorithm', 'trust-region-reflective', ...
%         'Display', 'off', ... 
%         'OptimalityTolerance', 1e-7);
% 
%     % Create objective function handle (passes x to the kernel)
%     obj_fun = @(x) varpro_kernel(x, t, data, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- 4. Non-linear Optimization ---
%     % Optimizes only [taus, shift]
%     [best_nonlin, ~] = fmincon(obj_fun, x0(:)', [], [], [], [], lb(:)', ub(:)', [], options);
% 
%     % --- 5. Final Reconstruction ---
%     % Retrieve optimal linear parameters (amps/bg), model curve, and Jacobian
%     [~, ~, lin_params, final_fit, Jacob] = varpro_kernel(best_nonlin, t, data, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- 6. Error Estimation (Covariance Matrix) ---
%     p_total = length(best_nonlin) + length(lin_params);
%     dof = N - p_total; % Degrees of freedom
% 
%     if strcmpi(cost_type, 'LS')
%         % Residual-based error for Least Squares
%         mse = sum((data - final_fit).^2) / dof;
%         cov_mat = mse * ( (Jacob' * Jacob) \ eye(p_total) );
%     else
%         % Information Matrix-based error for Poisson MLE
%         W = 1 ./ max(final_fit, 1e-6); % Weight matrix (inverse of model variance)
%         cov_mat = (Jacob' * (W .* Jacob)) \ eye(p_total);
%     end
% 
%     % Standard Errors (square root of the diagonal of covariance matrix)
%     errs = sqrt(diag(cov_mat));
% 
%     % Apply scaling for 95% Confidence Intervals
%     if strcmpi(error_type, '95CI')
%         if exist('tinv', 'file')
%             tcrit = tinv(0.975, dof); % Exact Student's t critical value
%         else
%             tcrit = 1.96; % Normal distribution approximation fallback
%         end
%         errs = errs * tcrit;
%     end
% 
%     % --- 7. Package Results ---
%     results.taus = best_nonlin(1:n_exp);
%     results.shift = best_nonlin(end);
%     results.amplitudes = lin_params(1:n_exp);
%     results.background = 0;
%     if fit_bg, results.background = lin_params(end); end
% 
%     % Store requested error values
%     results.err_vals.taus = errs(1:n_exp);
%     results.err_vals.shift = errs(n_exp + 1);
%     results.err_vals.amplitudes = errs(n_exp + 2 : 2*n_exp + 1);
%     if fit_bg, results.err_vals.background = errs(end); end
%     results.error_type_used = error_type;
% end
% 
% function [f, g, lin_params, model, Jacob] = varpro_kernel(x, t, data, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg)
%     % Internal Function: Calculates Cost, Gradient, and Jacobian
% 
%     taus = x(1:n_exp); shift = x(end);
%     N = length(t);
% 
%     % Shift IRF using Fourier Shift Theorem (Shift IRF in freq domain)
%     shift_kernel = exp(-2j * pi * freq * shift);
%     irf_s_f = irf_f .* shift_kernel; 
% 
%     n_lin = n_exp + double(fit_bg); % Columns for amplitudes + background
%     Phi = zeros(N, n_lin);          % Basis matrix
%     dPhi_dtau = zeros(N, n_exp);    % Basis derivatives for gradient
% 
%     % Construct Basis Functions (Convolved Decays)
%     for i = 1:n_exp
%         % Decay with Periodic Boundary Condition (Pulse Pile-up correction)
%         decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
%         decay_f = fft(decay);
%         Phi(:,i) = real(ifft(decay_f .* irf_s_f));
% 
%         % Analytic Derivative of the decay w.r.t Tau
%         d_dec = (t/taus(i)^2) .* decay;
%         dPhi_dtau(:,i) = real(ifft(fft(d_dec) .* irf_s_f));
%     end
%     if fit_bg, Phi(:, end) = 1; end % Add background column (ones)
% 
%     % Solve for Linear Parameters (Amplitudes/Background)
%     if strcmpi(cost_type, 'LS')
%         % Non-negative Least Squares
%         lin_params = lsqnonneg(Phi, data);
%         model = Phi * lin_params;
%         res = model - data;
%         f = sum(res.^2);    % Chi-square cost
%         w = 2 * res;        % Residual weights for gradient
%     else
%         % MLE (Poisson) using Iteratively Reweighted Least Squares (IRLS)
%         beta = lsqnonneg(Phi, data); 
%         for iter = 1:5 
%             model = Phi * beta;
%             model(model < 1e-6) = 1e-6; 
%             W = 1 ./ model;
%             beta = (Phi' * (W .* Phi)) \ (Phi' * (W .* data));
%             beta(beta < 0) = 0; % Enforce non-negativity
%         end
%         lin_params = beta;
%         model = Phi * lin_params;
%         model(model < 1e-6) = 1e-6;
%         f = sum(model - data .* log(model)); % Poisson log-likelihood
%         w = 1 - data./model; % Poisson weights for gradient
%     end
% 
%     % Calculate Derivative w.r.t IRF Shift
%     dirf_ds_f = irf_f .* (-2j*pi*freq) .* shift_kernel;
%     dM_ds = 0;
%     for i = 1:n_exp
%          decay_f = fft(exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i))));
%          dM_ds = dM_ds + lin_params(i) * real(ifft(decay_f .* dirf_ds_f));
%     end
% 
%     % Construct Full Jacobian Matrix for Error Estimation
%     % Columns: [dModel/dTaus, dModel/dShift, dModel/dAmps, (dModel/dBG)]
%     Jacob = zeros(N, length(x) + n_lin);
%     for i = 1:n_exp
%         Jacob(:, i) = lin_params(i) * dPhi_dtau(:,i); 
%     end
%     Jacob(:, n_exp + 1) = dM_ds; 
%     Jacob(:, n_exp + 2 : end) = Phi; 
% 
%     % Return analytic gradient if requested by fmincon
%     if nargout > 1
%         g = zeros(size(x));
%         for i = 1:n_exp
%             g(i) = sum(w .* Jacob(:, i));
%         end
%         g(end) = sum(w .* dM_ds);
%     end
% end


% function [results, final_fit] = fit_tcspc_varpro_bgoption_werrorbars(t, data, irf, x0, lb, ub, cost_type, fit_bg, error_type)
%     % Optimized TCSPC Solver using VARPRO with Error Estimation
%     % error_type: '1sigma' for standard error, '95CI' for 95% Confidence Interval
% 
%     t = t(:); data = data(:); irf = irf(:); 
%     n_exp = length(x0) - 1;
%     dt = t(2) - t(1);
%     T_pulse = t(end) + dt;
%     N = length(t);
% 
%     irf_f = fft(irf); 
%     freq = (0:N-1)'; 
%     freq(freq > N/2) = freq(freq > N/2) - N;
%     freq = freq / T_pulse;
% 
%     options = optimoptions('fmincon', ...
%         'SpecifyObjectiveGradient', true, ...
%         'Algorithm', 'trust-region-reflective', ...
%         'Display', 'off', ... 
%         'OptimalityTolerance', 1e-7);
% 
%     obj_fun = @(x) varpro_kernel(x, t, data, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     [best_nonlin, ~] = fmincon(obj_fun, x0(:)', [], [], [], [], lb(:)', ub(:)', [], options);
% 
%     % Final reconstruction and gradient for Jacobian/Hessian calculation
%     [~, ~, lin_params, final_fit, Jacob] = varpro_kernel(best_nonlin, t, data, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg);
% 
%     % --- Error Estimation ---
%     p_total = length(best_nonlin) + length(lin_params);
%     dof = N - p_total;
% 
%     if strcmpi(cost_type, 'LS')
%         mse = sum((data - final_fit).^2) / dof;
%         cov_mat = mse * ( (Jacob' * Jacob) \ eye(p_total) );
%     else
%         % MLE Poisson
%         W = 1 ./ max(final_fit, 1e-6);
%         cov_mat = (Jacob' * (W .* Jacob)) \ eye(p_total);
%     end
% 
%     % Standard Errors (1-sigma)
%     errs = sqrt(diag(cov_mat));
% 
%     % Scale errors if 95% CI is requested
%     if strcmpi(error_type, '95CI')
%         if exist('tinv', 'file') % Requires Statistics Toolbox
%             tcrit = tinv(0.975, dof);
%         else
%             tcrit = 1.96; % Normal approximation fallback
%         end
%         errs = errs * tcrit;
%     end
% 
%     % Map results
%     results.taus = best_nonlin(1:n_exp);
%     results.shift = best_nonlin(end);
%     results.amplitudes = lin_params(1:n_exp);
%     results.background = 0;
%     if fit_bg, results.background = lin_params(end); end
% 
%     % Store requested errors
%     results.err_vals.taus = errs(1:n_exp);
%     results.err_vals.shift = errs(n_exp + 1);
%     results.err_vals.amplitudes = errs(n_exp + 2 : 2*n_exp + 1);
%     if fit_bg, results.err_vals.background = errs(end); end
%     results.error_type_used = error_type;
% end
% 
% function [f, g, lin_params, model, Jacob] = varpro_kernel(x, t, data, irf_f, freq, T_pulse, n_exp, cost_type, fit_bg)
%     taus = x(1:n_exp); shift = x(end);
%     N = length(t);
% 
%     shift_kernel = exp(-2j * pi * freq * shift);
%     irf_s_f = irf_f .* shift_kernel; 
% 
%     n_lin = n_exp + double(fit_bg);
%     Phi = zeros(N, n_lin);
%     dPhi_dtau = zeros(N, n_exp);
% 
%     for i = 1:n_exp
%         decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
%         decay_f = fft(decay);
%         Phi(:,i) = real(ifft(decay_f .* irf_s_f));
%         d_dec = (t/taus(i)^2) .* decay;
%         dPhi_dtau(:,i) = real(ifft(fft(d_dec) .* irf_s_f));
%     end
%     if fit_bg, Phi(:, end) = 1; end
% 
%     % Solve Linear Parameters
%     if strcmpi(cost_type, 'LS')
%         lin_params = lsqnonneg(Phi, data);
%         model = Phi * lin_params;
%         res = model - data;
%         f = sum(res.^2);
%         w = 2 * res;
%     else
%         beta = lsqnonneg(Phi, data); 
%         for iter = 1:5 
%             model = Phi * beta;
%             model(model < 1e-6) = 1e-6; 
%             W = 1 ./ model;
%             beta = (Phi' * (W .* Phi)) \ (Phi' * (W .* data));
%             beta(beta < 0) = 0; 
%         end
%         lin_params = beta;
%         model = Phi * lin_params;
%         model(model < 1e-6) = 1e-6;
%         f = sum(model - data .* log(model));
%         w = 1 - data./model;
%     end
% 
%     % Shift derivative for Jacobian
%     dirf_ds_f = irf_f .* (-2j*pi*freq) .* shift_kernel;
%     dM_ds = 0;
%     for i = 1:n_exp
%          decay_f = fft(exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i))));
%          dM_ds = dM_ds + lin_params(i) * real(ifft(decay_f .* dirf_ds_f));
%     end
% 
%     % Full Jacobian [dM/dtau, dM/dshift, dM/damp, dM/dbg]
%     Jacob = zeros(N, length(x) + n_lin);
%     for i = 1:n_exp
%         Jacob(:, i) = lin_params(i) * dPhi_dtau(:,i); 
%     end
%     Jacob(:, n_exp + 1) = dM_ds; 
%     Jacob(:, n_exp + 2 : end) = Phi; 
% 
%     if nargout > 1
%         g = zeros(size(x));
%         for i = 1:n_exp
%             g(i) = sum(w .* Jacob(:, i));
%         end
%         g(end) = sum(w .* dM_ds);
%     end
% end