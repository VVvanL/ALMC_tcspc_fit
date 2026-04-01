function [irf,fp] = extract_IRF_2(data, min_prominence, center_pos, plot_fits)
% extracts an IRF estimate 
% data is the TCSPC data, 
% min_prominence is the minimum porminence relative to the maximum value 
% that a peak has to exceed to be included in the IRF
% center_pos is a number between 0 and 1 that describes how the 2*FWHM fit 
% range is placed around a the peak (0 means peak position is the start of 
% the range, 1 means peak position is the end of the range. plot_fits is 
% true/false

data = data(:);
x = 0:length(data)-1;

%
data(1)=data(2);
data(end)=data(end-1);

% smooth
data_sm = smooth(data);
data_sm(1:2)=data_sm(3);
data_sm(end-1:end)=data_sm(end-2);

% find peaks with prominence larger than min_prominence
[pks,locs,w,p] = findpeaks(data_sm/max(data_sm));
[pks,I] = sort(pks,'descend');

locs = locs(I);
w = w(I);
p = p(I);

sel_vec = p>=min_prominence;
n_peaks = sum(sel_vec);
% pks = pks(sel_vec);
% locs = locs(sel_vec);
% w = w(sel_vec);
% p = p(sel_vec);

% gradient
data_sm_grad = gradient(data_sm,x);

% find peaks
[pks,locs,w,p] = findpeaks(data_sm_grad);
[pks,I] = sort(pks,'descend');

locs = locs(I(1:n_peaks));
w = w(I(1:n_peaks));
p = p(I(1:n_peaks));

% get range around maxium
range = round([locs - 2*w*center_pos, locs + 2*w*(1-center_pos)]);

% make sure range does not exceed data
range(range<1)=1;
range(range>length(data))=length(data);

% get data within the ranges, fit it with double gaussian and generate IRf
irf = zeros(length(x),1);
fp = cell(n_peaks,1);
for i=1:n_peaks
    dmy = data_sm_grad(range(i,1):range(i,2));
    dmy = dmy - min(dmy);
    x_dmy = (range(i,1):range(i,2))';
    f = fit(x_dmy,dmy,'gauss2');
    fp{i} = f;
    if plot_fits
        figure;
        plot(f,x_dmy,dmy);
    end
    irf = irf + f(x);
end

% normalize irf
irf = round(irf/max(irf)*100000);
irf = irf/sum(irf);

if plot_fits
    figure;
    semilogy(data/max(data));
    hold on
    semilogy(irf/max(irf));
    hold off
end