function [results, final_fit, normalized_irf] = fit_tcspc_gauss_irf_varpro(t, data, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg, error_type)
% FIT_TCSPC_GAUSS_IRF_VARPRO
%
% PURPOSE:
%   Fits TCSPC data using Variable Projection (VARPRO) while simultaneously
%   modeling the IRF as a sum of Gaussians.
%
% INPUTS:
%   t, data: [N x 1] Vectors.
%   x0, lb, ub: [n_exp x 1] Initial guess/bounds for lifetimes.
%   x0_irf, lb_irf, ub_irf: Parameters for IRF Gaussians.
%       1 Gauss: [pos1, sig1]
%       2 Gauss: [pos1, sig1, pos2, sig2, rel_amp2]
%       3 Gauss: [pos1, sig1, pos2, sig2, rel_amp2, pos3, sig3, rel_amp3]

    % --- 1. Pre-processing ---
    t = t(:);           
    data = data(:);     
    n_exp = length(x0);     
    n_irf_params = length(x0_irf);
    dt = t(2) - t(1);           
    T_pulse = t(end) + dt; 
    N = length(t);              

    % Combine non-linear parameters: [taus; irf_params]
    x0_full = [x0(:); x0_irf(:)];
    lb_full = [lb(:); lb_irf(:)];
    ub_full = [ub(:); ub_irf(:)];

    % --- 2. Optimizer Configuration ---
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ... 
        'MaxIterations', 400, ...
        'OptimalityTolerance', 1e-7,...
        'StepTolerance', 1e-8);

    % Objective Function Handle
    obj_fun = @(x) varpro_kernel_gauss(x, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);

    % --- 3. Non-linear Optimization ---
    [best_nonlin, ~] = fmincon(obj_fun, x0_full, [], [], [], [], lb_full, ub_full, [], options);

    % --- 4. Final Reconstruction & Statistics ---
    [~, lin_params, final_fit, Jacob, normalized_irf] = varpro_kernel_gauss(best_nonlin, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);

    p_total = length(best_nonlin) + length(lin_params);
    dof = N - p_total;

    % Error Estimation (Covariance Matrix)
    if strcmpi(cost_type, 'LS')
        mse = sum((data - final_fit).^2) / dof;
        cov_mat = mse * ( (Jacob' * Jacob) \ eye(p_total) );
    else
        % Use the model variance for Poisson weighting
        W = 1 ./ max(final_fit, 1e-6);
        cov_mat = (Jacob' * (W .* Jacob)) \ eye(p_total);
    end

    errs = sqrt(diag(cov_mat));
    if strcmpi(error_type, '95CI')
        errs = errs * 1.96; 
    end

    % --- 5. Result Packaging ---
    results.taus = best_nonlin(1:n_exp);
    results.irf_gauss_params = best_nonlin(n_exp + 1 : n_exp + n_irf_params);
    results.amplitudes = lin_params(1:n_exp);
    results.background = 0;
    if fit_bg, results.background = lin_params(end); end
    results.estimated_irf = normalized_irf;

    % Map errors back to structure
    results.err_vals.taus = errs(1:n_exp);
    results.err_vals.irf_gauss = errs(n_exp + 1 : n_exp + n_irf_params);
    results.err_vals.amplitudes = errs(n_exp + n_irf_params + 1 : n_exp + n_irf_params + n_exp);
    if fit_bg, results.err_vals.background = errs(end); end
end

function [f, lin_params, model, Jacob, irf_norm] = varpro_kernel_gauss(x, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg)
    % Internal Function for VARPRO Kernel
    N = length(t);
    dt = t(2) - t(1);
    taus = x(1:n_exp);
    gauss_p = x(n_exp + 1 : end);

    % --- A. Synthesize IRF ---
    % First Gaussian: exp(-(t-pos)^2 / (2*sig^2))
    % Ensure sigma (gauss_p(2)) is not zero to avoid NaN
    sig1 = max(gauss_p(2), 1e-9);
    irf_raw = exp(-(t - gauss_p(1)).^2 / (2 * sig1^2));
    
    idx = 3;
    while idx <= n_irf_params
        p_pos = gauss_p(idx);
        p_sig = max(gauss_p(idx+1), 1e-9);
        p_amp = gauss_p(idx+2);
        irf_raw = irf_raw + p_amp * exp(-(t - p_pos).^2 / (2 * p_sig^2));
        idx = idx + 3;
    end
    
    % Normalize Area to 1
    irf_area = sum(irf_raw);% * dt;
    irf_norm = irf_raw / max(irf_area, 1e-12); 
    irf_f = fft(irf_norm);

    % --- B. Build Basis Matrix (Phi) ---
    n_lin = n_exp + double(fit_bg);
    Phi = zeros(N, n_lin);

    for i = 1:n_exp
        % Decay with pile-up correction: exp(-t/tau) / (1 - exp(-T/tau))
        tau_val = max(taus(i), 1e-9);
        decay = exp(-t/tau_val) / (1 - exp(-T_pulse/tau_val));
        
        % Convolution via FFT
        conv_res = real(ifft(fft(decay) .* irf_f)) * dt;
        Phi(:,i) = conv_res(1:N); 
    end
    if fit_bg, Phi(:, end) = 1; end

    % --- C. Linear Parameter Solution ---
    if strcmpi(cost_type, 'LS')
        % Least Squares
        lin_params = lsqnonneg(Phi, data);
        model = Phi * lin_params;
        f = sum((model - data).^2);
    else
        % Poisson MLE (using IRLS)
        beta = lsqnonneg(Phi, data);
        for iter = 1:3
            model = max(Phi * beta, 1e-6);
            W = 1 ./ model;
            beta = (Phi' * (W .* Phi)) \ (Phi' * (W .* data));
            beta(beta < 0) = 0;
        end
        lin_params = beta;
        model = max(Phi * lin_params, 1e-6);
        f = sum(model - data .* log(model));
    end

    % --- D. Numerical Jacobian ---
    if nargout > 3
        total_nonlin = length(x);
        Jacob = zeros(N, total_nonlin + n_lin);
        h = 1e-5; 
        for j = 1:total_nonlin
            x_h = x; x_h(j) = x_h(j) + h;
            % Recursive call (no Jacobian)
            [~, ~, m_h] = varpro_kernel_gauss(x_h, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
            Jacob(:, j) = (m_h - model) / h;
        end
        Jacob(:, total_nonlin + 1 : end) = Phi;
    end
end

% function [results, final_fit, normalized_irf] = fit_tcspc_gauss_irf_varpro(t, data, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg, error_type)
% % FIT_TCSPC_GAUSS_IRF_VARPRO
% %
% % PURPOSE:
% %   Fits TCSPC data using Variable Projection (VARPRO) while simultaneously
% %   modeling the IRF as a sum of Gaussians.
% %
% % INPUTS:
% %   t, data: [N x 1] Vectors.
% %   x0, lb, ub: [n_exp x 1] Initial guess/bounds for lifetimes.
% %   x0_irf, lb_irf, ub_irf: Parameters for IRF Gaussians.
% %       1 Gauss: [pos1, sig1]
% %       2 Gauss: [pos1, sig1, pos2, sig2, rel_amp2] ... and so on.
% 
%     % --- 1. Pre-processing & Dimensionality ---
%     t = t(:);           % Force column vector
%     data = data(:);     % Force column vector
%     n_exp = length(x0);     
%     n_irf_params = length(x0_irf);
%     dt = t(2) - t(1);           
%     T_pulse = t(end) + dt; 
%     N = length(t);              
% 
%     % Combine non-linear parameters
%     x0_full = [x0(:); x0_irf(:)];
%     lb_full = [lb(:); lb_irf(:)];
%     ub_full = [ub(:); ub_irf(:)];
% 
%     % --- 2. Optimizer Configuration ---
%     % Using 'sqp' algorithm which is robust for non-linear constraints/bounds.
%     options = optimoptions('fmincon', ...
%         'Algorithm', 'sqp', ...
%         'Display', 'off', ... 
%         'MaxIterations', 400, ...
%         'OptimalityTolerance', 1e-7);
% 
%     % Objective Function Handle
%     obj_fun = @(x) varpro_kernel_gauss(x, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
% 
%     % --- 3. Non-linear Optimization ---
%     [best_nonlin, ~] = fmincon(obj_fun, x0_full, [], [], [], [], lb_full, ub_full, [], options);
% 
%     % --- 4. Final Reconstruction & Statistics ---
%     [~, lin_params, final_fit, Jacob, normalized_irf] = varpro_kernel_gauss(best_nonlin, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
% 
%     p_total = length(best_nonlin) + length(lin_params);
%     dof = N - p_total;
% 
%     % Error Estimation (Covariance Matrix)
%     if strcmpi(cost_type, 'LS')
%         mse = sum((data - final_fit).^2) / dof;
%         cov_mat = mse * ( (Jacob' * Jacob) \ eye(p_total) );
%     else
%         W = 1 ./ max(final_fit, 1e-6);
%         cov_mat = (Jacob' * (W .* Jacob)) \ eye(p_total);
%     end
% 
%     errs = sqrt(diag(cov_mat));
%     if strcmpi(error_type, '95CI'), errs = errs * 1.96; end
% 
%     % --- 5. Result Packaging ---
%     results.taus = best_nonlin(1:n_exp);
%     results.irf_gauss_params = best_nonlin(n_exp + 1 : n_exp + n_irf_params);
%     results.amplitudes = lin_params(1:n_exp);
%     results.background = 0;
%     if fit_bg, results.background = lin_params(end); end
%     results.estimated_irf = normalized_irf;
% 
%     % Standard error mapping
%     results.err_vals.taus = errs(1:n_exp);
%     results.err_vals.irf_gauss = errs(n_exp + 1 : n_exp + n_irf_params);
%     results.err_vals.amplitudes = errs(n_exp + n_irf_params + 1 : n_exp + n_irf_params + n_exp);
%     if fit_bg, results.err_vals.background = errs(end); end
% end
% 
% function [f, lin_params, model, Jacob, irf_norm] = varpro_kernel_gauss(x, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg)
%     % Internal Function
%     t = t(:); data = data(:); % Ensure column vectors
%     taus = x(1:n_exp);
%     gauss_p = x(n_exp + 1 : end);
%     N = length(t);
%     dt = t(2) - t(1);
% 
%     % --- A. Synthesize IRF ---
%     % First Gaussian
%     irf_raw = exp(-(t - gauss_p(1)).^2 / (2 * gauss_p(2)^2));
%     % Additional Gaussians
%     idx = 3;
%     while idx <= n_irf_params
%         irf_raw = irf_raw + gauss_p(idx+2) * exp(-(t - gauss_p(idx)).^2 / (2 * gauss_p(idx+1)^2));
%         idx = idx + 3;
%     end
%     % Normalize Area to 1
%     irf_norm = irf_raw / (sum(irf_raw) * dt); 
%     irf_f = fft(irf_norm);
% 
%     % --- B. Build Basis Matrix (Phi) ---
%     n_lin = n_exp + double(fit_bg);
%     Phi = zeros(N, n_lin);
% 
%     for i = 1:n_exp
%         % Decay with pile-up correction
%         decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
%         % Convolution via FFT
%         Phi(:,i) = real(ifft(fft(decay) .* irf_f)) * dt; 
%     end
%     if fit_bg, Phi(:, end) = 1; end
% 
%     % --- C. Linear Parameter Solution ---
%     if strcmpi(cost_type, 'LS')
%         lin_params = lsqnonneg(Phi, data);
%         model = Phi * lin_params;
%         f = sum((model - data).^2);
%     else
%         % Poisson MLE (Iterative Reweighted Least Squares)
%         beta = lsqnonneg(Phi, data);
%         for iter = 1:3
%             model = max(Phi * beta, 1e-6);
%             W = 1 ./ model;
%             beta = (Phi' * (W .* Phi)) \ (Phi' * (W .* data));
%             beta(beta < 0) = 0;
%         end
%         lin_params = beta;
%         model = max(Phi * lin_params, 1e-6);
%         f = sum(model - data .* log(model));
%     end
% 
%     % --- D. Numerical Jacobian (Only computed when requested) ---
%     if nargout > 3
%         total_nonlin = length(x);
%         Jacob = zeros(N, total_nonlin + n_lin);
%         h = 1e-5; 
%         for j = 1:total_nonlin
%             x_h = x; x_h(j) = x_h(j) + h;
%             % Recursive call without requesting Jacobian to avoid infinite loop
%             [~, ~, m_h] = varpro_kernel_gauss(x_h, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
%             Jacob(:, j) = (m_h - model) / h;
%         end
%         Jacob(:, total_nonlin + 1 : end) = Phi;
%     end
% end


% function [results, final_fit, normalized_irf] = fit_tcspc_gauss_irf_varpro(t, data, x0, lb, ub, x0_irf, lb_irf, ub_irf, cost_type, fit_bg, error_type)
% % FIT_TCSPC_GAUSS_IRF_VARPRO
% %
% % PURPOSE:
% %   Performs TCSPC lifetime fitting using the Variable Projection (VARPRO) 
% %   algorithm. Unlike standard fits, this version treats the Instrument 
% %   Response Function (IRF) as a sum of N-Gaussians and optimizes their 
% %   shapes simultaneously with the lifetimes.
% %
% % --- INPUT SHAPES & DEFINITIONS ---
% % t         : [N x 1] vector. Time axis (usually nanoseconds).
% % data      : [N x 1] vector. Raw photon counts.
% % x0        : [n_exp x 1] vector. Initial guess for lifetimes [tau1, tau2, ...].
% % lb, ub    : [n_exp x 1] vectors. Lower/Upper bounds for lifetimes.
% % x0_irf    : [Variable x 1] vector. Parameters for the Gaussian IRF:
% %             - If 1 Gaussian: [pos1, sig1]
% %             - If 2 Gaussians: [pos1, sig1, pos2, sig2, rel_amp2]
% %             - If 3 Gaussians: [pos1, sig1, pos2, sig2, rel_amp2, pos3, sig3, rel_amp3]
% %             Note: 'sig' is the standard deviation. FWHM = 2.355 * sig.
% % lb_irf, ub_irf : Bounds for the IRF parameters (must match x0_irf length).
% % cost_type : 'LS' (Least Squares) or 'MLE' (Maximum Likelihood Estimation).
% % fit_bg    : boolean. true = fit a constant background; false = 0 background.
% % error_type: '1sigma' or '95CI'.
% %
% % --- OUTPUTS ---
% % results        : Structure containing fitted parameters and error estimates.
% % final_fit      : [N x 1] vector. The best-fit model curve.
% % normalized_irf : [N x 1] vector. The reconstructed IRF (area normalized to 1).
% 
%     % --- 1. Pre-processing & Dimensionality ---
%     t = t(:); data = data(:);
%     n_exp = length(x0);     
%     n_irf_params = length(x0_irf);
%     dt = t(2) - t(1);           
%     T_pulse = t(end) + dt; % Period of the laser pulses
%     N = length(t);              
% 
%     % Concatenate all non-linear parameters into one vector for fmincon
%     % Structure: [tau1, ..., tauN, irf_p1, ..., irf_pM]
%     x0_full = [x0(:); x0_irf(:)];
%     lb_full = [lb(:); lb_irf(:)];
%     ub_full = [ub(:); ub_irf(:)];
% 
%     % --- 2. Optimizer Configuration ---
%     % Using 'trust-region-reflective' as it handles bounds robustly.
%     options = optimoptions('fmincon', ...
%         'Algorithm', 'sqp', ...
%         'Display', 'off', ... 
%         'MaxIterations', 400, ...
%         'OptimalityTolerance', 1e-7);
% 
%     % Handle for the kernel which returns the cost (f)
%     obj_fun = @(x) varpro_kernel_gauss(x, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
% 
%     % --- 3. Non-linear Optimization ---
%     [best_nonlin, ~] = fmincon(obj_fun, x0_full, [], [], [], [], lb_full, ub_full, [], options);
% 
%     % --- 4. Final Reconstruction & Statistics ---
%     % Run the kernel one last time with optimal parameters to get the model and Jacobian
%     [~, lin_params, final_fit, Jacob, normalized_irf] = varpro_kernel_gauss(best_nonlin, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
% 
%     % Degrees of Freedom
%     p_total = length(best_nonlin) + length(lin_params);
%     dof = N - p_total;
% 
%     % Covariance Matrix calculation
%     if strcmpi(cost_type, 'LS')
%         mse = sum((data - final_fit).^2) / dof;
%         cov_mat = mse * ( (Jacob' * Jacob) \ eye(p_total) );
%     else
%         % For Poisson MLE, weights are the inverse of the model prediction
%         W = 1 ./ max(final_fit, 1e-6);
%         cov_mat = (Jacob' * (W .* Jacob)) \ eye(p_total);
%     end
% 
%     errs = sqrt(diag(cov_mat));
%     if strcmpi(error_type, '95CI'), errs = errs * 1.96; end
% 
%     % --- 5. Result Packaging ---
%     results.taus = best_nonlin(1:n_exp);
%     results.irf_gauss_params = best_nonlin(n_exp + 1 : n_exp + n_irf_params);
%     results.amplitudes = lin_params(1:n_exp);
%     results.background = 0;
%     if fit_bg, results.background = lin_params(end); end
%     results.estimated_irf = normalized_irf;
% 
%     % Error mapping
%     results.err_vals.taus = errs(1:n_exp);
%     results.err_vals.irf_gauss = errs(n_exp + 1 : n_exp + n_irf_params);
%     results.err_vals.amplitudes = errs(n_exp + n_irf_params + 1 : n_exp + n_irf_params + n_exp);
%     if fit_bg, results.err_vals.background = errs(end); end
% end
% 
% function [f, lin_params, model, Jacob, irf_norm] = varpro_kernel_gauss(x, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg)
%     % Internal Function: Maps non-linear parameters to linear ones and calculates cost.
% 
%     taus = x(1:n_exp);
%     gauss_p = x(n_exp + 1 : end);
%     N = length(t);
%     dt = t(2) - t(1);
% 
%     % --- Step A: IRF Synthesis ---
%     % Initialize with the first Gaussian (Amplitude = 1)
%     irf_raw = exp(-(t - gauss_p(1)).^2 / (2 * gauss_p(2)^2));
% 
%     % Add subsequent Gaussians if provided
%     % Parameter logic: [pos1, sig1, pos2, sig2, rel_amp2, pos3, sig3, rel_amp3...]
%     idx = 3;
%     while idx <= n_irf_params
%         irf_raw = irf_raw + gauss_p(idx+2) * exp(-(t - gauss_p(idx)).^2 / (2 * gauss_p(idx+1)^2));
%         idx = idx + 3;
%     end
% 
%     % Normalize: The area of the IRF must be 1 so that the fitted 
%     % amplitudes represent total photon counts accurately.
%     irf_norm = irf_raw / (sum(irf_raw) * dt); 
%     irf_f = fft(irf_norm);
% 
%     % --- Step B: Variable Projection (Linear Solve) ---
%     n_lin = n_exp + double(fit_bg);
%     Phi = zeros(N, n_lin);
% 
%     for i = 1:n_exp
%         % 1. Create decay basis function with pile-up (wrapping) correction
%         decay = exp(-t/taus(i)) / (1 - exp(-T_pulse/taus(i)));
%         % 2. Convolve with the synthetic IRF (Frequency domain multiplication)
%         Phi(:,i) = real(ifft(fft(decay) .* irf_f)) * dt; 
%     end
%     if fit_bg, Phi(:, end) = 1; end % Constant offset
% 
%     % Solve for amplitudes using non-negative least squares
%     if strcmpi(cost_type, 'LS')
%         lin_params = lsqnonneg(Phi, data);
%         model = Phi * lin_params;
%         f = sum((model - data).^2);
%     else
%         % Poisson Likelihood requires Iterative Reweighted Least Squares (IRLS)
%         beta = lsqnonneg(Phi, data);
%         for iter = 1:3 % Typically converges quickly for TCSPC
%             model = max(Phi * beta, 1e-6);
%             W = 1 ./ model;
%             beta = (Phi' * (W .* Phi)) \ (Phi' * (W .* data));
%             beta(beta < 0) = 0;
%         end
%         lin_params = beta;
%         model = max(Phi * lin_params, 1e-6);
%         % Poisson Deviance (simplified MLE cost)
%         f = sum(model - data .* log(model));
%     end
% 
%     % --- Step C: Numerical Jacobian ---
%     % Required for error estimation since the IRF shape is dynamic.
%     if nargout > 3
%         total_nonlin = length(x);
%         Jacob = zeros(N, total_nonlin + n_lin);
%         h = 1e-5; % Finite difference step
%         for j = 1:total_nonlin
%             x_h = x; x_h(j) = x_h(j) + h;
%             [~, ~, m_h] = varpro_kernel_gauss(x_h, t, data, T_pulse, n_exp, n_irf_params, cost_type, fit_bg);
%             Jacob(:, j) = (m_h - model) / h;
%         end
%         % The Jacobian for linear parameters is simply the basis matrix Phi
%         Jacob(:, total_nonlin + 1 : end) = Phi;
%     end
% end