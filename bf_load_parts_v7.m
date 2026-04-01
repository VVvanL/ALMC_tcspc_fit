function [data,omeMeta,seriesMetadata]=bf_load_parts_v7(file, dataseries, xrange, yrange, zslices, timeseries, channels, showWaitbar)
% Uses the bioformats plugin to read in a data from a dataseries within the
% specified data file. All indizes start at 1! zslices, timeseries and
% channels require 1d arrays that contain the indices of the zslices,
% timeseries, and channels that are supposed to be read in. If xrange, yrange
% zslices, timeseries or channels are set -1 then all the data in that
% dimension is read in. xrange and yrange are specified as [start pixel,
% width]
% for example: [data,omeMeta]=bf_load_parts_v7(strcat(path,file),1,-1,-1,1:100,1,[1 3])
% reads in the file specified by file in the directory path (ideally use
% uigetfile to open a dialogbox to get file and path). It reads in data
% from the first dataseries. It reads the entire xy plane. 
% It reads the first to the hundredth zplane, at
% the first timepoint of channel 1 and 3. Setting the 
% The output data is a matrix with up to 5 dimensions with the bioformats
% order XYZTC (1. dimension: x, 2.dimension: y, 3. dimension: z, 4.
% dimension: T, 5. dimension: C). If the last dimension only has one entry
% then the dimension is not kept, i.e. if Z,T,and C only have one entry
% then only output matrix has 2 dimensions. If Z and C would only have one
% entry but T would have 2 the matrix will have 4 dimensions [X Y 1 2]
% (note, only the C dimension is cut)
% omeMeta is the standardized bioformats metadata. For accessing the
% metadata infromation have a look at http://static.javadoc.io/org.openmicroscopy/ome-xml/5.5.4/ome/xml/meta/MetadataRetrieve.html
% for example: omeMeta.getPixelsSizeX(dataseries-1) returns the number of pixels in
% x-direction of dataseries (in bioformats all the indizes start at 0)
% omeMeta.getChannelName(0,3) returns the Channel name of Channel 3 for in
% dataseries 1
if nargin<8
    showWaitbar=false;
end

[XYZTC,meta]=bf_file_info(file);
[filepath,filename,fileext] = fileparts(file);
filename = strrep(filename,'_','\_');
bit_depth = str2num(char(meta.getPixelsSignificantBits(dataseries-1)));

if bit_depth == 8
    data_type = 'uint8';
elseif bit_depth == 16
    data_type = 'uint16';
elseif bit_depth == 32
    data_type = 'single';
elseif bit_depth == 64
    data_type = 'double';
else
    data_type = 'uint16';
end

% check if entire XYZTC of the selected timeseries should be loaded
if all(xrange == -1)
    xrange = [1 XYZTC(dataseries,1)];
end
if all(yrange == -1)
    yrange = [1, XYZTC(dataseries,2)];
end
if all(zslices == -1)
    zslices = 1:XYZTC(dataseries,3);
end
if all(timeseries == -1)
    timeseries = 1:XYZTC(dataseries,4);
end
if all(channels == -1)
    channels = 1:XYZTC(dataseries,5);
end

% bioformats can only load up to 2^31 pixels in one layer so if image is
% larger it has to be read in in multiple steps
n_tiles = ceil(xrange(2)*yrange(2)/2/10^9);
if n_tiles > 1
    dmy = zeros(n_tiles,2);
    dmy(1,1) = 1;
    dmy(1,2) = round(XYZTC(dataseries,1)/n_tiles);
    for i = 2:n_tiles
        dmy(i,1) = dmy(i-1,1) + dmy(i-1,2);
        dmy(i,2) = round(XYZTC(dataseries,1)*i/n_tiles);
    end
    xrange = dmy;
end

if dataseries>0 && dataseries<=size(XYZTC,1)
    if min(zslices)>0 && max(zslices)<=XYZTC(dataseries,3)
        if min(timeseries)>0 && max(timeseries)<=XYZTC(dataseries,4)
            if min(channels)>0 && max(channels)<=XYZTC(dataseries,5)

                if showWaitbar
                    f = waitbar(0, ['loading ',filename,fileext]);
                end
                
                % make data variable to fill
                data = zeros(yrange(2),xrange(2),max(size(zslices)),max(size(timeseries)),max(size(channels)),data_type);
                
                % load reader
                reader = bfGetReader(file);                
                reader.setSeries(dataseries-1);
                globalMetadata = reader.getGlobalMetadata();
                omeMeta=reader.getMetadataStore;
                seriesMetadata = reader.getSeriesMetadata();
                javaMethod('merge', 'loci.formats.MetadataTools', globalMetadata, seriesMetadata, 'Global ');
                
                n_steps = max(size(zslices))*max(size(timeseries))*max(size(channels)*n_tiles);
                current_step = 0;
                last_disp_val = 0;
                for z=1:max(size(zslices))
                    for t=1:max(size(timeseries))
                        for c=1:max(size(channels))
                            iPlane = reader.getIndex(zslices(z)-1, channels(c)-1, timeseries(t) - 1) + 1;
                            if n_tiles == 1
                                data(:,:,z,t,c) = (bfGetPlane(reader, iPlane,xrange(1),yrange(1),xrange(2),yrange(2)));
                                if showWaitbar
                                    if floor(current_step/n_steps *100)/100 - last_disp_val >= 0.05 
                                        last_disp_val = floor(current_step/n_steps *100)/100;
                                        waitbar(last_disp_val, f, ['loading ',filename,fileext])
                                        drawnow;
                                    end
                                end
                                current_step = current_step + 1;
                            else
                            % load the image in several tiles
                                for tile = 1 : n_tiles
                                    data(:,xrange(tile,1):xrange(tile,2),z,t,c) = (bfGetPlane(reader, iPlane,xrange(tile,1),yrange(1),xrange(tile,2)-xrange(tile,1)+1,yrange(2)));
                                    if showWaitbar
                                        if floor(current_step/n_steps *100)/100 - last_disp_val >= 0.05 
                                            last_disp_val = floor(current_step/n_steps *100)/100;
                                            waitbar(last_disp_val, f, ['loading ',filename,fileext])
                                            drawnow;
                                        end
                                    end
                                    current_step = current_step + 1;    
                                end
                            end
                        end
                    end
                end

                if showWaitbar
                    close(f);
                end
                reader.close();
            else    
                disp(['At least some of the channels do not exist in the current dataseries. No data was read in!']);
            end
        else
            disp(['At least some of the time points do not exist in the current dataseries. No data was read in!']);
        end
    else
        disp(['At least some of the z-planes do not exist in the current dataseries. No data was read in!']);
    end    
else    
    disp(['Dataseries ',num2str(dataseries), ' does not exist in the current dataset. No data was read in!']);
end

end