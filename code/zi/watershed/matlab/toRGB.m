function [ toRGB ] = toRGB( img )

size(img)
toRGB = zeros([size(img) 3]);
toRGB(:,:,1) = img;
toRGB(:,:,2) = img;
toRGB(:,:,3) = img;

end

