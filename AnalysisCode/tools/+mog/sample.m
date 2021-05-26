function x = sample(mog, n)
%MOG.SAMPLE sample from a 1d mixture of gaussians.
%
%See MOG.PDF for format.

if nargin < 2, n = 1; end

[M, NN] = size(mog);

cumul_modes = cumsum(mog(:,3:3:end), 2);

for i=n:-1:1
    % In each row, cumul_modes is the CMF of mode weights. rand > CMF results in a vector of 1s
    % followed by 0s where the first zero is the chosen index. For example, if it is [1 1 1 0 0 0 0 0]
    % that indicates that the 4th item should be chosen. Hence, the index of the chosen mode is
    % sum(rand <= CMF) + 1.
	mode(:, i) = sum(rand(M, 1) > cumul_modes, 2) + 1;
    
    % Get mu and sigma for each mode and draw a sample
    iMode = 3*(mode(:, i)-1)';
    iMu = sub2ind([M, NN], 1:M, iMode+1)';
    iSig = sub2ind([M, NN], 1:M, iMode+2)';
    
    x(:, i) = mog(iMu) + randn(M, 1) .* mog(iSig);
end
end