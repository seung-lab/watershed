function [ V ] = load_volume( filename, size )

f = fopen(filename);
V = fread(f, 'uint32');
V = reshape(V, size, size, size);

end

