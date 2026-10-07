function sk = prepare_spline_fast(spl, W, D)
%PREPARE_SPLINE_FAST Prepare spline coefficients for fast kernel evaluation.
%
%   sk = prepare_spline_fast(spl, W, D)
%
%   Inputs
%   ------
%   spl : spline pair
%       Two spline structures containing the real and imaginary
%       coefficient tables.
%
%   W : scalar
%       Kernel half-width.
%
%   D : scalar
%       Discretization factor.
%
%   Output
%   ------
%   sk : structure
%       Precomputed spline coefficient data for fast evaluation.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Kernel dimensions

    K = 2*W;

    if D <= 0 || D ~= round(D)
        error('D must be a positive integer.');
    end

    if W <= 0 || K ~= round(K)
        error('W must be a positive integer.');
    end

    if numel(spl) < 2
        error('spl must contain two spline structures.');
    end


    %% Spline coefficients

    Cr = spl(1).coefs;
    Ci = spl(2).coefs;

    nPieces = size(Cr, 1);

    if ~isequal(size(Cr), size(Ci))
        error('Real and imaginary spline coefficient tables must have the same size.');
    end

    if size(Cr, 2) ~= 4
        error('Expected cubic spline coefficients with 4 columns.');
    end


    %% Spline interval indices

    q = (0:D).';
    k = 1:K;

    idxRaw = k*D + 1 - q;
    idxRaw(1,:) = idxRaw(end,:);

    if any(idxRaw(:) < 1) || any(idxRaw(:) > nPieces)
        error('prepare_spline_fast:tableRange', ...
            ['W=%d, D=%d need spline pieces 1..%d, but spl only has %d pieces.\n' ...
             'Rebuild spl over a wider domain before increasing W.'], ...
            W, D, K*D, nPieces);
    end

    idx = idxRaw;
    idx = min(max(idx, 1), nPieces);   % now provably a no-op, kept as a safety net

    idxT = idx.';   % K x (D+1): k fastest, q slowest -- matches the layout below


    %% Precomputed coefficient tables

    sk.A_r = Cr(idx,1);
    sk.B_r = Cr(idx,2);
    sk.C_r = Cr(idx,3);
    sk.D_r = Cr(idx,4);

    sk.A_i = Ci(idx,1);
    sk.B_i = Ci(idx,2);
    sk.C_i = Ci(idx,3);
    sk.D_i = Ci(idx,4);


    %% Contiguous coefficient blocks

    % For a fixed q, all K taps now sit in one contiguous 4*K block.
    ABCD_r = zeros(4, K, D+1, 'like', Cr);
    ABCD_i = zeros(4, K, D+1, 'like', Ci);

    ABCD_r(1,:,:) = reshape(Cr(idxT,1), 1, K, D+1);
    ABCD_r(2,:,:) = reshape(Cr(idxT,2), 1, K, D+1);
    ABCD_r(3,:,:) = reshape(Cr(idxT,3), 1, K, D+1);
    ABCD_r(4,:,:) = reshape(Cr(idxT,4), 1, K, D+1);

    ABCD_i(1,:,:) = reshape(Ci(idxT,1), 1, K, D+1);
    ABCD_i(2,:,:) = reshape(Ci(idxT,2), 1, K, D+1);
    ABCD_i(3,:,:) = reshape(Ci(idxT,3), 1, K, D+1);
    ABCD_i(4,:,:) = reshape(Ci(idxT,4), 1, K, D+1);

    sk.ABCD_r = ABCD_r;
    sk.ABCD_i = ABCD_i;


    %% Store parameters

    sk.Dn = D;
    sk.K  = K;

end