%% Code for generating 3D fibers in helix morphology
clear; close all; 
savepath = '.\helix_m\'; mkdir(savepath)
sp = 2;                      % Number of helix fibers
sz0 = [80 80 185];           % Stack size of all groups
sz = [240 80 185];           % Stack size of one group

%% Parameter setting
numTurns = 4;              % Turns of helix
radius = 6;                % Radius of helix
pitch = 32;                % Pitch of helix
center = [40, 40, 30];     % Center coordinates
lineThickness = 4;         % Thickness of helix

stack00 = zeros(sz0, 'uint8');
%% Elastin channel generation
pth = '.\helix_tif_e\';  tifnm = pth(3:end-1);  mkdir(pth)

stack01 = fun_helix_gen(stack00, [sz0 numTurns radius pitch center+[0 0 pitch/4] lineThickness],  tifnm, sp, 1, 0);
stack02 = fun_helix_gen(stack00, [sz0 numTurns radius pitch center lineThickness],  tifnm, sp, 1, 0);
stack03 = fun_helix_gen(stack00, [sz0 numTurns radius pitch center lineThickness],  tifnm, sp, 1, 0);

stack3 = zeros(sz, 'uint8');
stack3 (1:80,:,:) = stack01;
stack3 (81:160,:,:) = stack02;
stack3 (161:240,:,:) = stack03;

heliceStackm1 = stack3;
save([savepath,'heliceStackm1'])
for i = 1 : sz(3)
    imwrite(heliceStackm1 (:,:,i),[pth,tifnm,num2str(i),'.tif'])
end

%% Collagen channel generation
pth = '.\helix_tif_c\';  tifnm = pth(3:end-1);  mkdir(pth)

stack01 = fun_helix_gen(stack00, [sz0 numTurns radius pitch center+[0 0 -pitch/4] lineThickness],  tifnm, sp, 1, 0);
stack02 = fun_helix_gen(stack00, [sz0 numTurns radius pitch center lineThickness],  tifnm, sp, 1, 0);
stack03 = fun_helix_gen(stack00, [sz0 numTurns radius pitch center lineThickness],  tifnm, sp, 1, 1);

stack3 = zeros(sz, 'uint8');
stack3 (1:80,:,:) = stack01; imrotate3(stack01,60,[0 0 1]);
stack3 (81:160,:,:) = imrotate3(stack02,90,[0 0 1]);
stack3 (161:240,:,:) = stack03;

heliceStackm2 = stack3;
save([savepath,'heliceStackm2'])
for i = 1 : sz(3)
    imwrite(heliceStackm2 (:,:,i),[pth,tifnm,num2str(i),'.tif'])
end

%% Save elastin and collagen fiber pictures
pth = '.\helix_tif_ec\';  tifnm = pth(3:end-1);   mkdir(pth)
heliceStackm3 = heliceStackm1 + heliceStackm2; clear heliceStackm1 heliceStackm2
save([savepath,'heliceStackm3'])

for i = 1 : sz(3)
    imwrite(heliceStackm3 (:,:,i),[pth,tifnm,num2str(i),'.tif'])
end
