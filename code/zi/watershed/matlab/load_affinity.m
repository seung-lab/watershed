function [ A ] = load_affinity( filename, size )

    f = fopen(filename);
    A = fread(f, 'single');
    A = reshape(A, size, size, size, 3);

end
