function make3dPlotAff(s, cmap)

aff = load_affinity('../../data/ws_test_256.raw', 256);
chan = max(aff, [], 4);
%chan = squeeze(chan(:,:,:,1));
Ix = toColormappedChannel(squeeze(s(128,:,:)), cmap, squeeze(chan(128,:,:)), 0.5); imshow(Ix);
Iy = toColormappedChannel(squeeze(s(:,128,:)), cmap, squeeze(chan(:,128,:)), 0.5); imshow(Iy);
Iz = toColormappedChannel(squeeze(s(:,:,128)), cmap, squeeze(chan(:,:,128)), 0.5); imshow(Iz);

%Ix = imadjust(Ix,[.2 .2 .2; .7 .7 .7],[]);
%Iy = imadjust(Iy,[.2 .2 .2; .7 .7 .7],[]);
%Iz = imadjust(Iz,[.2 .2 .2; .7 .7 .7],[]);

do3Dplot(Ix,Iy,Iz);
view([90 80 40]);
light('Position',[5 8 12],'Style','local')
axis off

end