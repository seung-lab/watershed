function make3dPlotChan()

s = h5read('../../data/x10y50z30_im_s2483_11475_6931_e2738_11730_7186.h5', '/main');
Ix = toRGB(squeeze(s(128,:,:)));
Iy = toRGB(squeeze(s(:,128,:)));
Iz = toRGB(squeeze(s(:,:,128)));

%Ix = imadjust(Ix,[.2 .2 .2; .7 .7 .7],[]);
%Iy = imadjust(Iy,[.2 .2 .2; .7 .7 .7],[]);
%Iz = imadjust(Iz,[.2 .2 .2; .7 .7 .7],[]);

do3Dplot(Ix,Iy,Iz);
view([90 80 40]);
light('Position',[5 8 12],'Style','local')
axis off

end