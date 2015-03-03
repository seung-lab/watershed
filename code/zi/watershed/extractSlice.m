function [ ] = extractSlice( filename, size, slice, cmap )

    f = fopen(filename);
    A = fread(f, 'uint32');
    A = reshape(A, size, size, size);
    I = toColormapped(squeeze(A(slice,:,:)), cmap);
    imwrite(I, [filename '.png']);

end
