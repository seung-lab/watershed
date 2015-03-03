function [ Out ] = toColormapped( Slice, ColorMap )

Modulo = size(ColorMap,1);
s = size(Slice);

Out = zeros([s 3]);

for x = 1:s(1)
    for y = 1:s(2)
        for c = 1:3
            Out(x,y,c) = ColorMap(rem(Slice(x,y),Modulo-1)+2,c);
            if Slice(x,y) == 0
                Out(x,y,c) = ColorMap(1,c);
            end;
        end;
    end;
end;

Out = uint8(Out * 255);

end

