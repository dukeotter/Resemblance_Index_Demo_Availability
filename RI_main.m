%% Here to load the 3D image stack of fiber-like structure and perform segmentation.
clear

tic

sz  = 512;
sz2 = 512;
totpic = 58;  % modify 1: stack depth
root_read = '.\Example_cancer\cancer\cancer'; % modify 2: path of reading
Duplim = 23; % modify 3: thickness-related distance parameter

stack = zeros(sz, sz2, totpic, 3);
for picnum = 1:totpic
    stack(:,:,picnum,:) = im2double( imread([root_read,num2str(picnum),'.tif']) );
end

stack_ela = stack(:,:,:,1);
stack_col = stack(:,:,:,2);    clear stack

fm_threshint = 0.05;  % modify 4: raw background intensity acquired by averaging intensity of several regions identified as background
fm_ratio = 1.1;  % modify 5: to provide an inmproved intensity threshold, and the factor '0.45' can be adjusted according to different sample

finalmask_ela = fm_gen(stack_ela, fm_threshint, fm_ratio);
% 'finalmask_ela' is the final binary mask selecting the elastin-fiber-only region
finalmask_col = fm_gen(stack_col, fm_threshint, fm_ratio);
% 'finalmask_col' is the final binary mask selecting the collagen-fiber-only region

%% Multi-parametric morpho-structural quantification loading
load(['.\datum_cancer\aSDEe.mat']);      load(['.\datum_cancer\aSDEc.mat'])
load(['.\datum_cancer\LDmatr3De.mat']);  load(['.\datum_cancer\LDmatr3Dc.mat'])
load(['.\datum_cancer\wavmatr3De.mat']); load(['.\datum_cancer\wavmatr3Dc.mat'])
load(['.\datum_cancer\Vmatre.mat']);     load(['.\datum_cancer\Vmatrc.mat'])

%% Closest inter-fiber distance searching: elastin to collagen (e2c)
pdsz_ = 11;    % modify 6: Median filtering parameter of scope on orientation theta (aSDE)
sigma_ = 1.2;  % modify 7: Median filtering parameter of sigma on orientation theta (aSDE)

aSDE_colmed = aSDE_medfilt(aSDEc,pdsz_);
[De2c, aSDEe2c, LDmatr3De2c, Vmatre2c, wavmatr3De2c] = closestsearching...
    (finalmask_ela,finalmask_col,Duplim,aSDE_colmed,LDmatr3Dc,Vmatrc,wavmatr3Dc);clear aSDE_colmed

De2c_ = De2c .* finalmask_ela;
Duplimit_ratio = 1;
De2c_(De2c_ > Duplim) = Duplim * Duplimit_ratio;

% Elastin to collagen difference map generation
Ae2c = abs((aSDEe - aSDEe2c) .* finalmask_ela);
Ae2c = (Ae2c > 90).*(180 - Ae2c) + (Ae2c <= 90).*(Ae2c);

LDe2c_ = abs((LDmatr3De - LDmatr3De2c) .* finalmask_ela);
Ve2c_ = abs((Vmatre - Vmatre2c) .* finalmask_ela);
We2c_ = abs((wavmatr3De - wavmatr3De2c) .* finalmask_ela);

%% Closest inter-fiber distance searching: collagen to elastin (c2e)
aSDE_elamed = aSDE_medfilt(aSDEe,pdsz_);
[Dc2e, aSDEc2e, LDmatr3Dc2e, Vmatrc2e, wavmatr3Dc2e] = closestsearching...
    (finalmask_col,finalmask_ela,Duplim,aSDE_elamed,LDmatr3De,Vmatre,wavmatr3De);clear aSDE_elamed

Dc2e_ = Dc2e.*finalmask_col;
Dc2e_(Dc2e_ > Duplim) = Duplim * Duplimit_ratio;

% Collagen to elastin difference map generation
Ac2e = abs((aSDEc - aSDEc2e) .* finalmask_col);
Ac2e = (Ac2e > 90).*(180 - Ac2e) + (Ac2e <= 90).*(Ac2e);

LDc2e_ = abs((LDmatr3Dc - LDmatr3Dc2e) .* finalmask_col);
Vc2e_ = abs((Vmatrc - Vmatrc2e) .* finalmask_col);
Wc2e_ = abs((wavmatr3Dc - wavmatr3Dc2e) .* finalmask_col);

%% Comprehensive RI parameter constructing
% Auxiliary matrix generation
fm_or = finalmask_ela|finalmask_col;           fm_and = finalmask_ela & finalmask_col;

% Inter-fiber difference map fusing
D_fused = (finalmask_ela-fm_and).*De2c_ + (finalmask_col-fm_and).*Dc2e_ + fm_and.*( (De2c_ + Dc2e_)./2 );      clear De2c_ Dc2e_
A_fused = (finalmask_ela-fm_and).*Ae2c  + (finalmask_col-fm_and).*Ac2e +  fm_and.*( (Ae2c  + Ac2e )./2 );      clear Ae2c Ac2e
V_fused = (finalmask_ela-fm_and).*Ve2c_ + (finalmask_col-fm_and).*Vc2e_ + fm_and.*( (Ve2c_ + Vc2e_)./2 );      clear Ve2c_ Vc2e_
W_fused = (finalmask_ela-fm_and).*We2c_ + (finalmask_col-fm_and).*Wc2e_ + fm_and.*( (We2c_ + Wc2e_)./2 );      clear We2c_ Wc2e_
LD_fused = (finalmask_ela-fm_and).*LDe2c_ + (finalmask_col-fm_and).*LDc2e_ + fm_and.*( (LDe2c_ + LDc2e_)./2 ); clear LDe2c_ LDc2e_

% Normalization scoring to generate resulting resemblance index matrix RI
VW_scoring_ratio = 4; % modify 8: to adjust directional variance and waviness scoring level
RI = ...
    sqrt( (1/(1 + (D_fused/Duplim)*(exp(1) - 1) )) ).*...
    sqrt( (1/(1 + LD_fused*        (exp(1) - 1) )) ).*...
    cos(A_fused/90*pi/2).*...
    1/2.*( (exp(-VW_scoring_ratio.*(V_fused-0.5))-exp(VW_scoring_ratio.*(V_fused-0.5)))./ (exp(-VW_scoring_ratio.*(V_fused-0.5))+exp(VW_scoring_ratio.*(V_fused-0.5)))  + 1).*...
    1/2.*( (exp(-VW_scoring_ratio.*(W_fused-0.5))-exp(VW_scoring_ratio.*(W_fused-0.5)))./ (exp(-VW_scoring_ratio.*(W_fused-0.5))+exp(VW_scoring_ratio.*(W_fused-0.5)))  + 1);

% Background control and sigular region control
RI(isnan(RI)) = 0;
RI = RI.*fm_or;
mean((RI(fm_or)))
%% Here to do post-processing.
% For post-processing, we prepare the 'pretty' images of orientation and
% waviness. In these 'pretty' images, the fused intensity image 
% is used to provide the contrast of fiber features, and the resemblance
% index maps are labeled by different colors to show the resemblance
% index information

% First, prepare the fused intensity image
stack_ela_re = stack_ela./(max(stack_ela(:)));
stack_col_re = stack_col./(max(stack_col(:)));
stackre = stack_ela_re + stack_col_re;            clear stack_ela_re stack_col_re

% Then, output pretty resemblance index images
root_outputpics = ('.\RI_cancer\') ; % modify 9:  path of outputing
mkdir(root_outputpics)

uplim = 1;  botlim = 0; bright = 0.99; dark = 0.01;
for mm = 1:totpic
    intensityima = stackre(:,:,mm);        RIima = RI(:,:,mm);
    prettyima = prettymap(imgaussfilt(RIima,sigma_),intensityima,[root_outputpics,'cancer_RI',num2str(mm)],jet(64),uplim,botlim,bright,dark);
    % Here 'prettyima' is the RI map.
    % 'prettymap.m' is the function used to create 'pretty' images.
    % 'RIima' is the resemblance index map, 'intensityima' is the fused
    % inter-channel image. 'jet(64)' designates the color scheme. 'uplim' 
    % and 'botlim' are the upper and bottom limits of the resemblance index, 
    % the range of which is from 0 to 180. 'bright' and 'dark' are used 
    % to enhance the contrast of the image
end
toc
