function flyOver(qq)

[d1 d2 d3] = size(qq);

close all;
figure;


for x = [1:d3 d3:1]
    imagesc(reshape(rem(qq(:, :, x)', 128), d1, d2));
    colormap jet;
    drawnow
    pause(0.02)
end

for x = [1:d2 d2:1]
    imagesc(reshape(rem(qq(:, x, :), 128), d1, d3)');
    colormap jet;
    drawnow
    pause(0.02)
end

for x = [1:d1 d1:1]
    imagesc(reshape(rem(qq(x, :, :), 128), d2, d3));
    colormap jet;
    drawnow
    pause(0.02)
end
