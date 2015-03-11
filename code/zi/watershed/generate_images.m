%load mycmap.mat;
chan = hdf5read('../../data/x10y50z30_im_s2483_11475_6931_e2738_11730_7186.h5', '/main');
aff  = load_affinity('../../data/ws_test_256.raw', 256);

pngExport(toRGB(squeeze(chan(128,:,:))), '../../../paper/Figures/chan.png', 16, 10);
pngExport(squeeze(aff(128,:,:,:)), '../../../paper/Figures/aff.png', 16, 10);

caff = max(aff,[],4);

all = hdf5read('../../data/x10y50z30_s2483_11475_6931_e2738_11730_7186.val.h5', '/main');
im = toColormappedChannel(squeeze(all(128,:,:)), cmap, squeeze(caff(128,:,:)), 0.5);
pngExport(im, '../../../paper/Figures/gt.png', 16, 10);

all = load_volume('./zished/experiments/felzenszwalb/10.000000.dat', 256);
im = toColormappedChannel(squeeze(all(128,:,:)), cmap, squeeze(caff(128,:,:)), 0.5);
pngExport(im, '../../../paper/Figures/fe1.png', 16, 10);

all = load_volume('./zished/experiments/felzenszwalb/0.500000.dat', 256);
im = toColormappedChannel(squeeze(all(128,:,:)), cmap, squeeze(caff(128,:,:)), 0.5);
pngExport(im, '../../../paper/Figures/fe2.png', 16, 10);


all = load_volume('./zished/experiments/square/3000.dat', 256);
im = toColormappedChannel(squeeze(all(128,:,:)), cmap, squeeze(caff(128,:,:)), 0.5);
pngExport(im, '../../../paper/Figures/all.png', 16, 10);

all = load_volume('./zished/experiments/watershed/basic.out', 256);
im = toColormappedChannel(squeeze(all(128,:,:)), cmap, squeeze(caff(128,:,:)), 0.5);
pngExport(im, '../../../paper/Figures/raw.png', 16, 10);

all = load_volume('./zished/experiments/watershed/minmax.out', 256);
im = toColormappedChannel(squeeze(all(128,:,:)), cmap, squeeze(caff(128,:,:)), 0.5);
pngExport(im, '../../../paper/Figures/minmax.png', 16, 10);

