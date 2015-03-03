function [ cmap ] = makeCmap( sz, intensity, zval )

cmap = rand(sz,3);

cmap(1,:) = zval;

end

