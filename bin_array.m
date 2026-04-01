function binned_data = bin_array(data, bin_size, bin_dims)
% BIN_ARRAY Bins an array along specified dimensions.
%
%   binned_data = BIN_ARRAY(data, bin_size, bin_dims)
%
%   Inputs:
%     data      - Input array of any dimension.
%     bin_size  - Scalar or vector of bin sizes for each bin_dim.
%     bin_dims  - Vector of dimensions to apply binning.
%
%   Output:
%     binned_data - Binned array, with reduced size along binned dimensions.

    % Input validation
    if isscalar(bin_size)
        bin_size = repmat(bin_size, size(bin_dims));
    elseif numel(bin_size) ~= numel(bin_dims)
        error('Length of bin_size must match length of bin_dims');
    end

    original_size = size(data);
    nd = ndims(data);

    % Permute the data so that binned dimensions come first
    all_dims = 1:max(nd, max(bin_dims));
    non_bin_dims = setdiff(all_dims, bin_dims, 'stable');
    permute_order = [bin_dims, non_bin_dims];
    data = permute(data, permute_order);

    % Get new shape after permutation
    new_size = size(data);
    bin_sizes = bin_size(:)'; % row vector

    % Calculate reshaped size for binning
    reshaped_size = [];
    binned_shape = [];
    for i = 1:length(bin_dims)
        d_len = new_size(i);
        bsz = bin_sizes(i);
        if mod(d_len, bsz) ~= 0
            warning('Dimension %d is not divisible by bin size %d. Truncating data.', bin_dims(i), bsz);
            d_len = floor(d_len / bsz) * bsz;
            idx = repmat({':'}, 1, nd);
            idx{i} = 1:d_len;
            data = data(idx{:});
        end
        reshaped_size = [reshaped_size, bsz, d_len/bsz];
        binned_shape = [binned_shape, d_len/bsz];
    end

    % Append non-binned dimensions
    reshaped_size = [reshaped_size, new_size(length(bin_dims)+1:end)];
    binned_shape = [binned_shape, new_size(length(bin_dims)+1:end)];

    % Reshape and bin (mean over bin axes)
    data = reshape(data, reshaped_size);

    % Compute mean over every second dimension (bin axes)
    for i = 1:length(bin_dims)
        data = sum(data, 2*i - 1); % 2*i-1 corresponds to bin axis
    end

    % Reshape to final size
    binned_data = reshape(data, binned_shape);

    % Inverse permute to original dimension order
    [~, inv_order] = sort(permute_order);
    binned_data = ipermute(binned_data, permute_order);
end
