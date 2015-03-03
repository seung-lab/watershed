function [  ] = do3Dplot( xs, ys, zs )

figure;
hold on;

xlabel('x');
ylabel('y');
zlabel('z');

xImage = [ 0    0  ;  0    0  ];
yImage = [-0.5  0.5; -0.5  0.5];
zImage = [ 0.5  0.5; -0.5 -0.5];

surf(xImage,yImage,zImage, 'CData', xs, 'FaceColor','texturemap');

xImage = [-0.5 -0.5;  0.5  0.5];
yImage = [-0.5  0.5; -0.5  0.5];
zImage = [ 0    0  ;  0    0  ];

surf(xImage,yImage,zImage, 'CData', ys, 'FaceColor','texturemap');

xImage = [-0.5 -0.5;  0.5  0.5];
yImage = [ 0    0  ;  0    0  ];
zImage = [ 0.5 -0.5;  0.5 -0.5];

surf(xImage,yImage,zImage, 'CData', zs, 'FaceColor','texturemap');

end