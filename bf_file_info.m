function [XYZTC,meta,image_name]=bf_file_info(filepath)
% uses

reader = bfGetReader(filepath);
n_series=reader.getSeriesCount;

meta=reader.getMetadataStore;

XYZTC=zeros(n_series,5);
series_names = cell(n_series,1);
image_name = [];

for i=1:n_series
    image_name = [image_name, meta.getImageName(i-1)];
    reader.setSeries(i-1);

    XYZTC(i,1)=reader.getSizeX;
    XYZTC(i,2)=reader.getSizeY;
    XYZTC(i,3)=reader.getSizeZ;
    XYZTC(i,4)=reader.getSizeT;
    XYZTC(i,5)=reader.getSizeC;

end
reader.close;
image_name = string(image_name);
end