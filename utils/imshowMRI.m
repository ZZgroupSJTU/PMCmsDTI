function imshowMRI(image_matrix, scale, tile_shape, titles,cmap)
%
%  ismrm_imshow(image_matrix, [low high], [image_tile_rows
%  image_tile_columns], titles )
%
%  Displays a collection of images in a tiled figure
%
%  INPUT:
%    - image_matrix    : 2-D, displays single image
%                        3-D, third dimension indexes subimages
%                        More then 3D, format into 3D
%    - [low high]      : low displays as black, high as white
%    - [image_tile_rows image_tile_columns] : specifies how to tile
%                                             subimages in display
%    - titles  : {'title1', 'title2',...}
%
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%--------------------------------------------------------
%  Updated to enable setting the colormap  -Zhiyong Zhang
%  Updated with fig pos setting -Zhiyong Zhang
%  (zhiyongxmu@gmail.com)
%---------------------------------------------------------
%%
% Validate input

image_matrix=reshape(image_matrix,size(image_matrix,1),size(image_matrix,2),[]);

switch nargin
    case 1
        scale = [];
        tile_shape = [];
        titles = [];
        cmap='gray';
    case 2
        tile_shape = [];
        titles = [];
        cmap='gray';
    case 3
        titles = [];
        cmap='gray';
    case 4
        cmap='gray';
    case 5

    otherwise
        error('valid number of arguments: 1-4');
end

assert( ndims(image_matrix) > 1 && ndims(image_matrix) < 5, 'image_matrix must have 2 or 3 dimensions')

if isempty(scale),
    scale = [min(image_matrix(:)), max(image_matrix(:))];
end

if isempty(tile_shape)
    tile_shape = [1 size(image_matrix,3)];
end

if (prod(tile_shape)>=size(image_matrix,3))
    image_matrix=cat(3,image_matrix,zeros(size(image_matrix,1),size(image_matrix,2),prod(tile_shape)-size(image_matrix,3)));
else
    error('image tile rows x columns must larger than the 3rd dim extent of image_matrix');
end

if ~isempty(titles),
    assert(numel(titles) == size(image_matrix,3), 'number of titles must equal 3rd dim extent of image_matrix');
end

%%
% Set figure shapes and strides
border = 4; % pixels
size_scale = 1;

num_rows = tile_shape(1);
num_cols = tile_shape(2);

im_shape = [size(image_matrix,1) size(image_matrix,2)] * size_scale;
im_stride = im_shape + border;
frame_shape = im_stride .* tile_shape + border ;

screensize = get( groot, 'Screensize' );

size_scale=min(0.9*screensize([4,3])./frame_shape);
if size_scale<1
    im_shape=im_shape*size_scale;
    im_stride=im_stride*size_scale;
    frame_shape=frame_shape*size_scale;
    border=border*size_scale;
end

% normalized shape and strides for subplot function
im_norm_shape = im_shape ./ frame_shape;
im_norm_stride = im_stride ./ frame_shape;
border_norm=[border border]./frame_shape;

if isempty(titles),
    title_margin = 0;
else
    title_margin = num_rows*50;
end


figpos=[screensize(3)/2-frame_shape(2)/2,screensize(4)/2-frame_shape(1)/2,frame_shape(2),frame_shape(1)+title_margin];


figure_handle = gcf;
set(figure_handle,'Units','Pixels');
set(figure_handle,'Position',figpos);
set(figure_handle,'Color',[1 1 1]);
% set(figure_handle,'Color','none'); 

%%
% Traverse & display images
curr_index = 1;
for row_index = 1:num_rows,
    for col_index = 1:num_cols,
        subplot('Position',[border_norm(2)+im_norm_stride(2) * (col_index-1) border_norm(1)+im_norm_stride(1) * (num_rows - row_index)  im_norm_shape(2) im_norm_shape(1)]);
        if scale ==0
            imshow(image_matrix(:,:,curr_index),[]); 
            colormap(gca,cmap);
        else
            imshow(image_matrix(:,:,curr_index), scale); 
            colormap(gca,cmap);
        end
        if ~isempty(titles),
            title(titles(curr_index),'color','y');
        end
        curr_index = curr_index + 1;
    end
end

