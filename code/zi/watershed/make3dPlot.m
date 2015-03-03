chan = h5read('../../data/x10y50z30_im_s2483_11475_6931_e2738_11730_7186.h5', '/main');
Ix = toColormappedChannel(squeeze(s(128,:,:)), cmap, squeeze(chan(128,:,:)), 0.2); imshow(I);
Iy = toColormappedChannel(squeeze(s(:,128,:)), cmap, squeeze(chan(:,128,:)), 0.2); imshow(I);
Iz = toColormappedChannel(squeeze(s(:,:,128)), cmap, squeeze(chan(:,:,128)), 0.2); imshow(I);
do3Dplot(Ix,Iy,Iz);
light('Position',[10 12 14],'Style','local')
axis off