function pngExport( im, fname, skip, border )

I = im(skip+1:end-skip,skip+1:end-skip,:);
I(1:border,:,:) = 255;
I(:,1:border,:) = 255;
I(end+1-border:end,:,:) = 255;
I(:,end+1-border:end,:) = 255;
imwrite(I, fname);

end

