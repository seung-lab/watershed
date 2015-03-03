function flyOver(qq)

[d1 d2 d3] = size(qq);

close all;
figure;

cmap = rand(1024,3);
cmap(1,:) = 0;

colormap(cmap);

for x = [1:d3 d3:1]
    imagesc(reshape(rem(qq(:, :, x)', 1024), d1, d2));
    drawnow
    pause(0.02)
end

for x = [1:d2 d2:1]
    imagesc(reshape(rem(qq(:, x, :), 1024), d1, d3)');
    drawnow
    pause(0.02)
end

for x = [1:d1 d1:1]
    imagesc(reshape(rem(qq(x, :, :), 1024), d2, d3));
    drawnow
    pause(0.02)
end
