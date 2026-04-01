function save_fit_data_as_tiff(r_im, data, filepath)

%% data output to be saved in ome.tiff files
% total counts, fit mask, taus, amplitudes, shift if fitted, background if fitted, error vals of
% all fitted parameters ratio taus, ratio amplitudes, ratio total_counts

n_exp = size(r_im.taus,3);
if max(r_im.shift,[],'all') == min(r_im.shift,[],'all')
    fit_shift = false;
else 
    fit_shift = true;
end
if isfield(r_im, 'background')
    fit_bg = true;
else
    fit_bg = false;
end


switch n_exp
    case 1
        n_datalayers = 7 + n_exp*4 ;
    case 2
        n_datalayers = 10 + n_exp*4;
    case 3
        n_datalayers = 10 + n_exp*4;
end
channel_names = cell(1,n_datalayers);

data_2_tiff = zeros([size(r_im.taus,[1:2]),n_datalayers]);

cl = 1; % keep track of current layer

% total counts

data_2_tiff(:,:,cl) = sum(data,3);
channel_names{cl} = "total counts";
cl = cl + 1;

data_2_tiff(:,:,cl) = sum(r_im.raw_data,3);
channel_names{cl} = "total counts binned xy";
cl = cl + 1;

% fit mask
data_2_tiff(:,:,cl) = r_im.mask;
channel_names{cl} = "fit mask";
cl = cl + 1;

% taus
for i = 1:n_exp
    data_2_tiff(:,:,cl) = r_im.taus(:,:,i);
    channel_names{cl} = strcat("lifetime tau",num2str(i)," (ns)");
    cl = cl + 1;
end

% amplitudes
for i = 1:n_exp
    data_2_tiff(:,:,cl) = r_im.amplitudes(:,:,i);
    channel_names{cl} = strcat("amplitude ",num2str(i));
    cl = cl + 1;
end

% background if fitted
if fit_bg
    data_2_tiff(:,:,cl) = r_im.background;
    channel_names{cl} = "background";
else
    data_2_tiff(:,:,cl) = zeros(size(r_im.mask));
    channel_names{cl} = "background (not fitted)";
end
cl = cl+1;

% shift, if fitted
if fit_shift 
    data_2_tiff(:,:,cl) = r_im.shift;
    channel_names{cl} = "t0 shift";
else
    data_2_tiff(:,:,cl) = r_im.shift;
    channel_names{cl} = "t0 shift (not fitted)";
end
cl = cl+1;

% tau error vals
for i = 1:n_exp
    data_2_tiff(:,:,cl) = r_im.err_vals.taus(:,:,i);
    channel_names{cl} = strcat("lifetime error tau",num2str(i)," (ns)");
    cl = cl + 1;
end

% amplitdue error vals
for i = 1:n_exp
    data_2_tiff(:,:,cl) = r_im.err_vals.amplitudes(:,:,i);
    channel_names{cl} = strcat("amplitude error ",num2str(i));
    cl = cl + 1;
end

% background if fitted
if fit_bg
    data_2_tiff(:,:,cl) = r_im.err_vals.background;
    channel_names{cl} = "background error";
    cl = cl+1;
end

% shift, if fitted
if fit_shift 
    data_2_tiff(:,:,cl) = r_im.err_vals.shift;
    channel_names{cl} = "t0 shift error";
    cl = cl+1;
end

% ratios
if n_exp>1
    ratio = zeros(size(r_im.taus,1:2));
    num = r_im.taus(:,:,1);
    den = r_im.taus(:,:,2);
    ratio(r_im.mask) = num(r_im.mask)./den(r_im.mask);
    data_2_tiff(:,:,cl) = ratio;
    channel_names{cl} = "ratio liftetime 1 to lifetime 2";
    cl = cl + 1;

    num = r_im.amplitudes(:,:,1);
    den = r_im.amplitudes(:,:,2);
    ratio(r_im.mask) = num(r_im.mask)./den(r_im.mask);
    data_2_tiff(:,:,cl) = ratio;
    channel_names{cl} = "ratio amplitude 1 to amplitude 2";
    cl = cl + 1;

    num = r_im.amplitudes(:,:,1).*r_im.taus(:,:,1);
    den = r_im.amplitudes(:,:,2).*r_im.taus(:,:,2);
    ratio(r_im.mask) = num(r_im.mask)./den(r_im.mask);
    data_2_tiff(:,:,cl) = ratio;
    channel_names{cl} = "ratio photon counts fit component 1 to fit component 2";
    cl = cl + 1;
end

%% save the data to tiff


data_32bit = single(data_2_tiff); 
[h, w, numC] = size(data_32bit); 

metadata = createMinimalOMEXMLMetadata(data_32bit, 'XYCZT');

for i = 0:numC-1
    chanID = sprintf('Channel:0:%d', i);
    metadata.setChannelID(chanID, 0, i);
    metadata.setChannelName(channel_names{i+1}, 0, i);
end

bfsave(data_32bit, filepath, 'metadata', metadata);