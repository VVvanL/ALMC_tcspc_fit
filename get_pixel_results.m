function pixel_results = get_pixel_results(r_im, i, j)
% GET_PIXEL_RESULTS
%
% PURPOSE:
%   Automatically detects the fit configuration from the image results (r_im)
%   and extracts the parameters for a specific pixel (i, j).
%
% LOGIC:
%   - Uses 'isfield' to check if background was fitted.
%   - Uses 'squeeze' to convert 3D image slices into the 1D vectors 
%     expected by the fitting and plotting functions.
%
% INPUTS:
%   r_im - The image-wide results structure containing 2D/3D arrays.
%   i, j - The row and column indices of the pixel to extract.

    % --- 1. Non-Linear Parameters ---
    % 'squeeze' is critical here. r_im.taus(i,j,:) is a 1x1xN array;
    % plotting functions expect a simple Nx1 or 1xN vector.
    pixel_results.taus = squeeze(r_im.taus(i,j,:));
    pixel_results.shift = r_im.shift(i,j);
    
    % --- 2. Linear Parameters ---
    pixel_results.amplitudes = squeeze(r_im.amplitudes(i,j,:));
    
    % --- 3. Background Detection (Self-Sensing) ---
    % We check the structure itself to see if background data exists.
    if isfield(r_im, 'background')
        pixel_results.background = r_im.background(i,j);
    else
        % If the field doesn't exist, it wasn't fitted; we default to 0.
        pixel_results.background = 0;
    end

    % --- 4. Error Values ---
    % We replicate the sub-structure 'err_vals' for the single pixel.
    pixel_results.err_vals.taus = squeeze(r_im.err_vals.taus(i,j,:));
    pixel_results.err_vals.amplitudes = squeeze(r_im.err_vals.amplitudes(i,j,:));
    pixel_results.err_vals.shift = r_im.err_vals.shift(i,j);
    
    % Check for background error field existence
    if isfield(r_im.err_vals, 'background')
        pixel_results.err_vals.background = r_im.err_vals.background(i,j);
    end

    % --- 5. Error Metadata ---
    % Carrying over the type of error (1sigma vs 95CI) if it was stored.
    if isfield(r_im, 'error_type_used')
        pixel_results.error_type_used = r_im.error_type_used;
    else
        pixel_results.error_type_used = 'Unknown';
    end
end