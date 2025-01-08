function [D, aSDE_, LDmatr3D_, Vmatr_, wavmatr3D_] = closestsearching(finalmask1,finalmask2,Duplim,aSDE2,LDmatr3D2,Vmatr2,wavmatr3D2)
% Find closest pixel to generate the minimum distance matrix and
% closest-pixel parameter matrice, from finalmask1 to finalmask2.
% If fm2s is zero matrix, let D be Duplim and let other parameter matrix be
% NaN to be set zero after final scoring.
[sz,sz2,totpic] = size(finalmask1);

% Output minimum distance map initialization
D  = zeros(sz,sz2,totpic);
ind_ = zeros(sz,sz2,totpic);

aSDE_ = zeros(sz,sz2,totpic);
LDmatr3D_ = zeros(sz,sz2,totpic);
Vmatr_ = zeros(sz,sz2,totpic);
wavmatr3D_ = zeros(sz,sz2,totpic);

%% calculation process
for kk = 1 : totpic
    % 2D slicing of input and output matrix
    % Input matrix slicing
    fm1s = finalmask1(:,:,kk);
    fm2s = finalmask2(:,:,kk);

    % Output matrix slicing
    Ds = zeros(sz,sz2); 
    inds = zeros(sz,sz2);

    aSDE_s= zeros(sz,sz2);
    LDmatr3D_s= zeros(sz,sz2);
    Vmatr_s= zeros(sz,sz2);
    wavmatr3D_s= zeros(sz,sz2);

    if any(fm2s(:))
    % 2D slicing of destination parameter maxtrix 
    aSDEs= aSDE2(:,:,kk);
    LDmatr3Ds= LDmatr3D2(:,:,kk);
    Vmatrs= Vmatr2(:,:,kk);
    wavmatr3Ds= wavmatr3D2(:,:,kk);

    [Prow,   Pcol] = find(fm2s);  P  = [Prow   Pcol]; %  P: two column vectors
    [PQrow, PQcol] = find(fm1s);  PQ = [PQrow PQcol]; % PQ: two column vectors

    [k,dist] = dsearchn(P,PQ);  % Search closest pixel from PQ to P.

    % Distance map
    Ds(sub2ind([sz sz2],PQrow,PQcol)) = dist;
    D(:,:,kk) = Ds;

    % Index map
    indseq = sub2ind([sz sz2], P(k,1), P(k,2)); % Index vector of closest pixel from PQ to P.
    inds(fm1s) = indseq; % 2D index matrix of closest pixel from PQ to P.
    ind_(:,:,kk) = inds; % Index matrix of closest pixel from PQ to P.

    % Orientation map's substraction item
    aSDE_s(fm1s) = aSDEs(indseq);
    aSDE_(:,:,kk) = aSDE_s;

    % LD map's substraction item
    LDmatr3D_s(fm1s) = LDmatr3Ds(indseq);
    LDmatr3D_(:,:,kk) = LDmatr3D_s;

    % Variance map's substraction item
    Vmatr_s(fm1s) = Vmatrs(indseq);
    Vmatr_(:,:,kk) = Vmatr_s;

    % Waviness map's substraction item
    wavmatr3D_s(fm1s) = wavmatr3Ds(indseq);
    wavmatr3D_(:,:,kk) = wavmatr3D_s;

    else
    %% esle
    % Distance map
    D(:,:,kk) = fm1s*Duplim;

    % Index map
    inds(fm1s) = NaN;
    ind_(:,:,kk) = inds; % Index matrix of closest pixel from PQ to P.

    % Orientation map's substraction item
    aSDE_s(fm1s) = NaN;
    aSDE_(:,:,kk) = aSDE_s;

    % L map's substraction item
    LDmatr3D_s(fm1s) = NaN;
    LDmatr3D_(:,:,kk) = LDmatr3D_s;

    % Variance map's substraction item
    Vmatr_s(fm1s) = NaN;
    Vmatr_(:,:,kk) = Vmatr_s;

    % Waviness map's substraction item
    wavmatr3D_s(fm1s) = NaN;
    wavmatr3D_(:,:,kk) = wavmatr3D_s;
    end
end
end