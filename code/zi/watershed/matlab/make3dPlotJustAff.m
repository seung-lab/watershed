function make3dPlotJustAff()

aff = load_affinity('../../data/ws_test_256.raw', 256);

%chan = squeeze(chan(:,:,:,1));
Ix = squeeze(aff(128,:,:,:));
Iy = squeeze(aff(:,128,:,:));
Iz = squeeze(aff(:,:,128,:));

Ix = uint8(Ix * 255);
Iy = uint8(Iy * 255);
Iz = uint8(Iz * 255);

%Ix = imadjust(Ix,[.2 .2 .2; .7 .7 .7],[]);
%Iy = imadjust(Iy,[.2 .2 .2; .7 .7 .7],[]);
%Iz = imadjust(Iz,[.2 .2 .2; .7 .7 .7],[]);

do3Dplot(Ix,Iy,Iz);
view([90 80 40]);
light('Position',[5 8 12],'Style','local')
axis off

end