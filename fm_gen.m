function finalmask = fm_gen (stack, threshint, ratio1)
[sz,sz2,totpic] = size(stack);

maska3d = stack > threshint;
% maska3d is a raw binary mask based on the raw background intensity
meanint = mean(stack(maska3d));
% calculate the mean intensity within the raw binary mask
threshinta = ratio1*meanint;

maska3da = stack > threshinta;
% 'maska3da' is an improved binary mask
finalmaska = maska3d.*maska3da;
% 'finalmaska' is the binary mask acquired based on the above two masks

finalmask = zeros(sz,sz2,totpic);
for j = 1:totpic
    Maskf1 = finalmaska(:,:,j);
    Maskf2 = miprmdebrisn(Maskf1,15);
    % 'miprmdebrisn.m' is a function which is able to remove very tiny structures
    finalmask(:,:,j) = Maskf2;
end
clear maska3d maska3da finalmaska
finalmask = logical(finalmask);
end
% 'finalmask' is the final binary mask selecting the fiber-only region