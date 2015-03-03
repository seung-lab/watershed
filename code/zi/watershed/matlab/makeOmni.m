function makeOmni( outf, segf, sz, dendf )

vol = load_volume(segf, sz);
f = fopen(dendf, 'r');
d_raw = fread(f, 'single');
fclose(f);
s = size(d_raw, 1) / 3;
d_raw = reshape(d_raw, 3, s);
d_raw = d_raw';
dend = uint32(d_raw(:,1:2));
dendV = d_raw(:,3);

hdf5write(outf, '/main', uint32(vol), '/dend', uint32(dend), '/dendValues', single(dendV));

end

