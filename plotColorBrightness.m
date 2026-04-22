function [rgbImage, hcb] = plotColorBrightness(colorData, brightData, cmapName, gamma, plotTitle, cbLabel, prctileGate, varargin)
    % PLOTCOLORBRIGHTNESS Maps colorData to a colormap and brightData to intensity.
    % 
    % Inputs:
    %   colorData:   2D array for color (values <= 0 are forced to black)
    %   brightData:  2D array for brightness/intensity
    %   cmapName:    (Optional) String name of a colormap. Default 'parula'.
    %   gamma:       (Optional) Gamma for brightness. Default 0.5.
    %   plotTitle:   (Optional) String for the main figure title.
    %   cbLabel:     (Optional) String for the colorbar side label.
    %   prctileGate: (Optional) Scalar [0-50). Gates extreme values.
    %   'ScaleBar', pos, lenPx, lenUm: (Optional) Adds a scale bar.

    % Handle optional inputs
    if nargin < 3 || isempty(cmapName), cmapName = 'parula'; end 
    if nargin < 4 || isempty(gamma), gamma = 0.5; end 
    if nargin < 5 || isempty(plotTitle), plotTitle = ''; end
    if nargin < 6 || isempty(cbLabel), cbLabel = 'Color Data Scale'; end
    if nargin < 7 || isempty(prctileGate), prctileGate = 0; end

    % 1. Create a mask for values > 0
    validMask = colorData > 0;
    normC = zeros(size(colorData));
    
    % 2. Scale Color Data (with gating)
    minV = 0; maxV = 1; 
    if any(validMask(:))
        vals = colorData(validMask);
        minV = prctile(vals, prctileGate);
        maxV = prctile(vals, 100 - prctileGate);
        
        if maxV > minV
            clippedC = colorData;
            clippedC(validMask & (colorData < minV)) = minV;
            clippedC(validMask & (colorData > maxV)) = maxV;
            normC(validMask) = (clippedC(validMask) - minV) / (maxV - minV);
        else
            normC(validMask) = 1; 
        end
    end

    % 3. Normalize Brightness Data (with gating)
    bVals = brightData(:);
    bMin = prctile(bVals, prctileGate);
    bMax = prctile(bVals, 100 - prctileGate);
    if bMax > bMin
        clippedB = brightData;
        clippedB(brightData < bMin) = bMin;
        clippedB(brightData > bMax) = bMax;
        normB = (clippedB - bMin) / (bMax - bMin);
    else
        normB = ones(size(brightData)); 
    end
    normB = normB .^ gamma;

    % 4. Map to RGB
    numColors = 256;
    cmap = feval(cmapName, numColors);
    idx = round(normC * (numColors - 1)) + 1;
    
    r = reshape(cmap(idx, 1), size(colorData)) .* normB .* validMask;
    g = reshape(cmap(idx, 2), size(colorData)) .* normB .* validMask;
    b = reshape(cmap(idx, 3), size(colorData)) .* normB .* validMask;
    rgbImage = cat(3, r, g, b);

    % 5. Scale Bar Logic
    addScaleBar = false;
    for i = 1:length(varargin)
        if ischar(varargin{i}) && strcmpi(varargin{i}, 'ScaleBar')
            sbPos = varargin{i+1};
            sbLenPx = varargin{i+2};
            sbLenUm = varargin{i+3};
            addScaleBar = true;
            break;
        end
    end

    if addScaleBar
        [h, w, ~] = size(rgbImage);
        thickness = max(2, round(h/120)); % Slightly thinner bar for aesthetics
        margin = round(w/25);
        
        % Determine pixel coordinates for the bar
        switch lower(sbPos)
            case 'tl'
                rRange = margin:(margin + thickness);
                cRange = margin:(margin + sbLenPx);
            case 'tr'
                rRange = margin:(margin + thickness);
                cRange = (w - margin - sbLenPx):(w - margin);
            case 'bl'
                rRange = (h - margin - thickness):(h - margin);
                cRange = margin:(margin + sbLenPx);
            case 'br'
                rRange = (h - margin - thickness):(h - margin);
                cRange = (w - margin - sbLenPx):(w - margin);
        end
        
        % Draw white bar into image matrix
        rgbImage(rRange, cRange, :) = 1;
    end

    % 6. Display and Colorbar
    imshow(rgbImage);
    axis image;
    if ~isempty(plotTitle), title(plotTitle); end
    
    colormap(gca, cmap); 
    hcb = colorbar;
    if maxV > minV, clim([minV, maxV]); else, clim([minV - 0.1, minV + 0.1]); end
    ylabel(hcb, cbLabel);

    % 7. Add Scale Bar Text (Closer positioning)
    if addScaleBar
        hold on;
        txt = [num2str(sbLenUm), ' \mum'];
        textX = mean(cRange);
        
        % Vertical offset set to a very small constant to keep it tight
        if contains(lower(sbPos), 'b')
            % Text slightly above the bar
            textY = min(rRange) - (thickness * 0.5); 
            valign = 'bottom';
        else
            % Text slightly below the bar
            textY = max(rRange) + (thickness * 0.5); 
            valign = 'top';
        end
        
        text(textX, textY, txt, 'Color', 'white', 'FontWeight', 'bold', ...
            'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', valign);
        hold off;
    end
end

% function [rgbImage, hcb] = plotColorBrightness(colorData, brightData, cmapName, gamma, plotTitle, cbLabel, prctileGate)
%     % PLOTCOLORBRIGHTNESS Maps colorData to a colormap and brightData to intensity.
%     % 
%     % Inputs:
%     %   colorData:   2D array for color (values <= 0 are forced to black)
%     %   brightData:  2D array for brightness/intensity
%     %   cmapName:    (Optional) String name of a colormap. Default 'parula'.
%     %   gamma:       (Optional) Gamma for brightness. Default 0.5.
%     %   plotTitle:   (Optional) String for the main figure title.
%     %   cbLabel:     (Optional) String for the colorbar side label.
%     %   prctileGate: (Optional) Scalar [0-50). Gates out the lowest and highest 
%     %                x percent of values. Default is 0 (no gating).
% 
%     % Handle optional inputs
%     if nargin < 3 || isempty(cmapName), cmapName = 'parula'; end 
%     if nargin < 4 || isempty(gamma), gamma = 0.5; end 
%     if nargin < 5 || isempty(plotTitle), plotTitle = ''; end
%     if nargin < 6 || isempty(cbLabel), cbLabel = 'Color Data Scale'; end
%     if nargin < 7 || isempty(prctileGate), prctileGate = 0; end
% 
%     % 1. Create a mask for values > 0 (ignore background/zeros)
%     validMask = colorData > 0;
%     normC = zeros(size(colorData));
% 
%     % Initialize min/max for colorbar scaling
%     minV = 0; maxV = 1; 
% 
%     % 2. Scale Color Data (excluding zeros and applying percentile gating)
%     if any(validMask(:))
%         vals = colorData(validMask);
% 
%         % Calculate gated thresholds
%         minV = prctile(vals, prctileGate);
%         maxV = prctile(vals, 100 - prctileGate);
% 
%         if maxV > minV
%             % Clip values to the gated range before normalizing
%             clippedC = colorData;
%             clippedC(validMask & (colorData < minV)) = minV;
%             clippedC(validMask & (colorData > maxV)) = maxV;
% 
%             normC(validMask) = (clippedC(validMask) - minV) / (maxV - minV);
%         else
%             normC(validMask) = 1; 
%         end
%     end
% 
%     % 3. Normalize and apply Gamma to Brightness Data (with gating)
%     bVals = brightData(:);
%     bMin = prctile(bVals, prctileGate);
%     bMax = prctile(bVals, 100 - prctileGate);
% 
%     if bMax > bMin
%         % Clip brightness values
%         clippedB = brightData;
%         clippedB(brightData < bMin) = bMin;
%         clippedB(brightData > bMax) = bMax;
% 
%         normB = (clippedB - bMin) / (bMax - bMin);
%     else
%         normB = ones(size(brightData)); 
%     end
%     normB = normB .^ gamma; % Apply Gamma Correction
% 
%     % 4. Map to RGB using the chosen colormap
%     numColors = 256;
%     cmap = feval(cmapName, numColors);
%     idx = round(normC * (numColors - 1)) + 1;
% 
%     % Combine components
%     r = reshape(cmap(idx, 1), size(colorData)) .* normB .* validMask;
%     g = reshape(cmap(idx, 2), size(colorData)) .* normB .* validMask;
%     b = reshape(cmap(idx, 3), size(colorData)) .* normB .* validMask;
%     rgbImage = cat(3, r, g, b);
% 
%     % 5. Display the result
%     imshow(rgbImage);
%     axis image; 
% 
%     % 6. Apply custom Title
%     if ~isempty(plotTitle)
%         title(plotTitle);
%     end
% 
%     % 7. Add and configure the Colorbar
%     colormap(gca, cmap); 
%     hcb = colorbar;
% 
%     % Set color limits based on gated thresholds
%     if maxV > minV
%         clim([minV, maxV]); 
%     else
%         clim([minV - 0.1, minV + 0.1]); 
%     end
% 
%     ylabel(hcb, cbLabel);
% end

% function [rgbImage, hcb] = plotColorBrightness(colorData, brightData, cmapName, gamma, plotTitle, cbLabel)
%     % PLOTCOLORBRIGHTNESS Maps colorData to a colormap and brightData to intensity.
%     % 
%     % Inputs:
%     %   colorData:  2D array for color (values <= 0 are forced to black)
%     %   brightData: 2D array for brightness/intensity (must be >= 0)
%     %   cmapName:   (Optional) String name of a MATLAB colormap (e.g., 'turbo')
%     %   gamma:      (Optional) Gamma for brightness. Default is 0.5.
%     %   plotTitle:  (Optional) String for the main figure title
%     %   cbLabel:    (Optional) String for the colorbar side label
% 
%     % Handle optional inputs
%     if nargin < 3 || isempty(cmapName), cmapName = 'parula'; end 
%     if nargin < 4 || isempty(gamma), gamma = 0.5; end 
%     if nargin < 5 || isempty(plotTitle), plotTitle = ''; end
%     if nargin < 6 || isempty(cbLabel), cbLabel = 'Color Data Scale'; end
% 
%     % 1. Create a mask for values > 0 (ignore background/zeros)
%     validMask = colorData > 0;
%     normC = zeros(size(colorData));
% 
%     % Initialize min/max for colorbar scaling
%     minV = 0; maxV = 1; 
% 
%     % 2. Scale Color Data (excluding zeros)
%     if any(validMask(:))
%         vals = colorData(validMask);
%         minV = min(vals);
%         maxV = max(vals);
% 
%         if maxV > minV
%             normC(validMask) = (colorData(validMask) - minV) / (maxV - minV);
%         else
%             normC(validMask) = 1; % Handle flat non-zero data
%         end
%     end
% 
%     % 3. Normalize and apply Gamma to Brightness Data
%     bMin = min(brightData(:));
%     bMax = max(brightData(:));
%     if bMax > bMin
%         normB = (brightData - bMin) / (bMax - bMin);
%     else
%         normB = ones(size(brightData)); 
%     end
%     normB = normB .^ gamma; % Apply Gamma Correction
% 
%     % 4. Map to RGB using the chosen colormap
%     numColors = 256;
%     cmap = feval(cmapName, numColors);
%     idx = round(normC * (numColors - 1)) + 1;
% 
%     % Combine components
%     r = reshape(cmap(idx, 1), size(colorData)) .* normB .* validMask;
%     g = reshape(cmap(idx, 2), size(colorData)) .* normB .* validMask;
%     b = reshape(cmap(idx, 3), size(colorData)) .* normB .* validMask;
%     rgbImage = cat(3, r, g, b);
% 
%     % 5. Display the result with fixed aspect ratio
%     imshow(rgbImage);
%     axis image; % CRITICAL: Keeps aspect ratio without stretching
% 
%     % 6. Apply custom Title
%     if ~isempty(plotTitle)
%         title(plotTitle);
%     end
% 
%     % 7. Add and configure the Colorbar
%     colormap(gca, cmap); 
%     hcb = colorbar;
% 
%     if maxV > minV
%         clim([minV, maxV]); 
%     else
%         clim([minV - 1, minV + 1]); 
%     end
% 
%     % Apply custom colorbar description
%     ylabel(hcb, cbLabel);
% end


% function [rgbImage, hcb] = plotColorBrightness(colorData, brightData, cmapName, gamma)
%     % PLOTCOLORBRIGHTNESS Maps colorData to a colormap and brightData to intensity.
%     % 
%     % Inputs:
%     %   colorData: 2D array for color (values <= 0 are forced to black)
%     %   brightData: 2D array for brightness/intensity (must be >= 0)
%     %   cmapName: (Optional) String name of a MATLAB colormap (e.g., 'turbo')
%     %   gamma: (Optional) Gamma correction for brightness. Default is 0.5 (brighter).
%     %          Use 1.0 for linear (no change), < 1 to brighten, > 1 to darken.
% 
%     if nargin < 3 || isempty(cmapName), cmapName = 'parula'; end 
%     if nargin < 4 || isempty(gamma), gamma = 0.5; end % Default gamma to brighten
% 
%     % 1. Create a mask for values > 0 (to ignore background/zeros)
%     validMask = colorData > 0;
%     normC = zeros(size(colorData));
% 
%     % Initialize min/max for colorbar scaling
%     minV = 0; 
%     maxV = 1; 
% 
%     % 2. Scale Color Data (excluding zeros)
%     if any(validMask(:))
%         vals = colorData(validMask);
%         minV = min(vals);
%         maxV = max(vals);
% 
%         if maxV > minV
%             normC(validMask) = (colorData(validMask) - minV) / (maxV - minV);
%         else
%             normC(validMask) = 1; % If all non-zero values are identical
%         end
%     end
% 
%     % 3. Normalize and apply Gamma to Brightness Data
%     bMin = min(brightData(:));
%     bMax = max(brightData(:));
% 
%     if bMax > bMin
%         normB = (brightData - bMin) / (bMax - bMin);
%     else
%         normB = ones(size(brightData)); 
%     end
% 
%     % Apply Gamma Correction
%     normB = normB .^ gamma;
% 
%     % 4. Map to RGB using the chosen colormap
%     numColors = 256;
%     cmap = feval(cmapName, numColors);
%     idx = round(normC * (numColors - 1)) + 1; % Convert 0-1 to 1-256 indices
% 
%     % Extract RGB components and apply brightness + zero mask
%     r = reshape(cmap(idx, 1), size(colorData)) .* normB .* validMask;
%     g = reshape(cmap(idx, 2), size(colorData)) .* normB .* validMask;
%     b = reshape(cmap(idx, 3), size(colorData)) .* normB .* validMask;
% 
%     rgbImage = cat(3, r, g, b);
% 
%     % 5. Display the result
%     imshow(rgbImage);
%     axis image;
% 
%     % 6. Add and configure the Colorbar
%     colormap(gca, cmap); 
%     hcb = colorbar;
% 
%     if maxV > minV
%         clim([minV, maxV]); 
%     else
%         clim([minV - 1, minV + 1]); 
%     end
% 
%     ylabel(hcb, 'Color Data Scale (Excluding 0)');
%     title(['Gamma = ', num2str(gamma)]);
% end


% function [rgbImage, hcb] = plotColorBrightness(colorData, brightData, cmapName)
%     % PLOTCOLORBRIGHTNESS Maps colorData to a colormap and brightData to intensity.
%     % 
%     % Inputs:
%     %   colorData: 2D array for color (values <= 0 are forced to black)
%     %   brightData: 2D array for brightness/intensity
%     %   cmapName: (Optional) String name of a MATLAB colormap (e.g., 'turbo', 'parula')
%     %
%     % Outputs:
%     %   rgbImage: The final processed RGB image
%     %   hcb: Handle to the colorbar for further customization
% 
%     if nargin < 3, cmapName = 'parula'; end 
% 
%     % 1. Create a mask for values > 0 (to ignore background/zeros)
%     validMask = colorData > 0;
%     normC = zeros(size(colorData));
% 
%     % Initialize min/max for colorbar scaling
%     minV = 0; 
%     maxV = 1; 
% 
%     % 2. Scale Color Data (excluding zeros)
%     if any(validMask(:))
%         vals = colorData(validMask);
%         minV = min(vals);
%         maxV = max(vals);
% 
%         if maxV > minV
%             normC(validMask) = (colorData(validMask) - minV) / (maxV - minV);
%         else
%             normC(validMask) = 1; % If all non-zero values are identical
%         end
%     end
% 
%     % 3. Normalize Brightness Data (0 to 1)
%     bMin = min(brightData(:));
%     bMax = max(brightData(:));
%     if bMax > bMin
%         normB = (brightData - bMin) / (bMax - bMin);
%     else
%         normB = ones(size(brightData)); 
%     end
% 
%     % 4. Map to RGB using the chosen colormap
%     numColors = 256;
%     cmap = feval(cmapName, numColors);
%     idx = round(normC * (numColors - 1)) + 1; % Convert 0-1 to 1-256 indices
% 
%     % Extract RGB components and apply brightness + zero mask
%     r = reshape(cmap(idx, 1), size(colorData)) .* normB .* validMask;
%     g = reshape(cmap(idx, 2), size(colorData)) .* normB .* validMask;
%     b = reshape(cmap(idx, 3), size(colorData)) .* normB .* validMask;
% 
%     rgbImage = cat(3, r, g, b);
% 
%     % 5. Display the result
%     imshow(rgbImage);
%     axis image;
% 
%     % 6. Add and configure the Colorbar
%     colormap(gca, cmap); % Set colormap for this specific axes
%     hcb = colorbar;
% 
%     if maxV > minV
%         clim([minV, maxV]); % Updates the colorbar ticks to match your data
%     else
%         % If data is flat, just label the single value
%         clim([minV - 1, minV + 1]); 
%     end
% 
%     ylabel(hcb, 'Color Data Scale (Excluding 0)');
% end


% function rgbImage = plotColorBrightness(colorData, brightData, cmapName)
%     % Default colormap if not provided
%     if nargin < 3, cmapName = 'parula'; end 
% 
%     % 1. Create a mask for values > 0
%     validMask = colorData > 0;
% 
%     % 2. Initialize normalized color array (all zeros/black by default)
%     normC = zeros(size(colorData));
% 
%     % 3. Scale only the non-zero values
%     if any(validMask)
%         vals = colorData(validMask);
%         minV = min(vals);
%         maxV = max(vals);
% 
%         % Avoid division by zero if all non-zero values are the same
%         if maxV > minV
%             normC(validMask) = (colorData(validMask) - minV) / (maxV - minV);
%         else
%             normC(validMask) = 1; % Default to top of colormap
%         end
%     end
% 
%     % 4. Normalize brightness (0 to 1)
%     normB = (brightData - min(brightData(:))) / (max(brightData(:)) - min(brightData(:)));
% 
%     % 5. Map to RGB using the chosen colormap
%     cmap = feval(cmapName, 256);
%     idx = round(normC * 255) + 1;
% 
%     r = reshape(cmap(idx, 1), size(colorData));
%     g = reshape(cmap(idx, 2), size(colorData));
%     b = reshape(cmap(idx, 3), size(colorData));
% 
%     % 6. Apply brightness AND the zero-mask (forcing colorData == 0 to black)
%     rgbImage = cat(3, r .* normB .* validMask, ...
%                       g .* normB .* validMask, ...
%                       b .* normB .* validMask);
% 
%     % Display
%     imshow(rgbImage);
%     axis image;
% end
