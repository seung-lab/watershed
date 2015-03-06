function [ P, R ] = loadPrecisionRecall( fname )

f = fopen(fname, 'r');
d = fread(f, 'double');
P = d(1:2:end-1);
R = d(2:2:end);

end

