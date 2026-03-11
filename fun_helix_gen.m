function stack2 = fun_helix_gen(stack, paras, tifnm, sp, num1, rot1)

% ================ Parameter setting ================
sz = paras(1:3);           % Stack size
numTurns = paras(2+2);            % Turns of helix
radius = paras(3+2);              % Radius of helix
pitch = paras(4+2);               % Pitch of helix
center = [paras(5+2) paras(6+2) paras(7+2)];  % Center coordinates
lineThickness = paras(8+2);       % Thickness of helix
intensity = 255;                  % Intensity of helix
% ===================================================

%% Generating helix point
t = linspace(0, 2*pi*numTurns, 5000);
if rot1
    x = radius * sin(t) + center(1);  % Rotate
    y = radius * cos(t) + center(2);
else
    x = radius * cos(t) + center(1);
    y = radius * sin(t) + center(2);
end
    z = (pitch/(2*pi)) * t + center(3);

% Constraining generating range
x = min(max(round(x), 1), sz(1));
y = min(max(round(y), 1), sz(2));
z = min(max(round(z), 1), sz(3));

%% Generating helix point
heliceStackm = stack;

% Linewith procession
r = floor(lineThickness/2);
[dx, dy, dz] = ndgrid(-r:r, -r:r, -r:r);
distMask = dx.^2 + dy.^2 + dz.^2 <= r^2;

% Plotting all points
for i = 1:length(t)
    cx = x(i); cy = y(i); cz = z(i);

    for k = 1:numel(dx)
        if distMask(k)
            nx = cx + dx(k);
            ny = cy + dy(k);
            nz = cz + dz(k);

            if nx >= 1 && nx <= sz(1) && ...
                    ny >= 1 && ny <= sz(2) && ...
                    nz >= 1 && nz <= sz(3)
                heliceStackm(nx, ny, nz) = intensity;
            end
        end
    end
end

stack2 = heliceStackm;
end