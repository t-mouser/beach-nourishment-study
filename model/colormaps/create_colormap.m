%% COLORMAP - automatically creates a matrix of RGB values 
% (expressed as percentages) for use in a colormap

%% Input an array of RGB values below
% There should be at least three colors
colors=[230,97,0;255,255,255;93,58,155]; % enter colors here

%% Caculate and propogate the colormaps
arrnum = 256; % determines the number of rows in the matrix
numberofcolors = size(colors,1)-1; % determines the number of colors in the map
genmap = zeros(arrnum,3); % creates an empty array for faster processing
colorbuckets = arrnum./(numberofcolors); % is used in automation calculations

for i=1:arrnum % loop to create the colormap matrix
    if i == 1
        genmap(1,:) = colors(1,:); % sets first row to be the initial color values
    else 
        genmap(i,1) = genmap(i-1,1) + (colors(floor((i-1)/colorbuckets)+2,1) - colors(floor((i-1)/colorbuckets)+1,1))/colorbuckets;
        genmap(i,2) = genmap(i-1,2) + (colors(floor((i-1)/colorbuckets)+2,2) - colors(floor((i-1)/colorbuckets)+1,2))/colorbuckets;
        genmap(i,3) = genmap(i-1,3) + (colors(floor((i-1)/colorbuckets)+2,3) - colors(floor((i-1)/colorbuckets)+1,3))/colorbuckets;
         
    end
    genmap(:)=abs(genmap(:));
end

%% These Variables Need Unique Names for Each Colormap!!
op_equal_cmap = genmap / arrnum; %converts RGB values in the percentages for MatLab to interpret
op_equal_cmap_rev = flip(op_equal_cmap,1); %flip the colormap in case a reverse color order is necessary

%% Save the Colormap
% IMPORTANT: make these unique values when saving multiple colormaps
save('constrast_equal.mat','op_equal_cmap','op_equal_cmap_rev') %save colormap for use in other files
