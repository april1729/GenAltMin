listing = dir('*.png');
[~, rect] = imcrop(imread(listing(7).name));
for im=1:length(listing)
I = imread(listing(im).name);
I2 = imcrop(I,rect);
imwrite(I2,listing(im).name);
end