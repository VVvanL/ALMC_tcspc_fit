function h = plot_intensity_image(im_data)

h = figure; 
imagesc(im_data);
colormap('parula');
c = colorbar;
c.Label.String = 'Intensity Counts';
axis equal