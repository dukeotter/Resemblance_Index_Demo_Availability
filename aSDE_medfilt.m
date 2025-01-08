function aSDE_medfilted = aSDE_medfilt(aSDE,pdsz)
% This function performs median filtering on orientation theta matrix.
[sz,sz2,totpic] = size(aSDE);
aSDEmedtmp = zeros(sz+2*pdsz,sz2+2*pdsz,totpic);
for i = 1:totpic
    aSDEmedtmp(:,:,i) = medfilt2(padarray(aSDE(:,:,i),[pdsz pdsz]),[pdsz pdsz]);
end
aSDE_medfilted = aSDEmedtmp(1+pdsz:sz+pdsz,1+pdsz:sz+pdsz,:);
end