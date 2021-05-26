function L = logpdf(x, mog)
%MOG.LOGPDF compute log pdf of 1d mixture-of-gaussians more stably than
%log(mog.pdf())
%
%See mog.pdf for more information.
%
%To use use vectorization to evaluate many PDFs at once, mog can be a [M x 3N] matrix specifying M
%different MoG distributions. In this case, x should be [M x 1]
%
%Note: unexpected behavior may occur if any mode's sigma is too small, as in a delta distribution

mus = mog(:, 1:3:end);
sigmas = mog(:, 2:3:end);
pis = mog(:, 3:3:end);

log_mode_probs = log(pis) - (x - mus).^2 ./ sigmas.^2 / 2 + log(2*pi*sigmas);
L = logsumexp(log_mode_probs, 2);

end

function s = logsumexp(a, dim)
% Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
% Default is dim = 1 (columns).
% logsumexp(a, 2) will sum across rows instead of columns.
% Unlike matlab's "sum", it will not switch the summing direction
% if you provide a row vector.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

if nargin < 2
  dim = 1;
end

% subtract the largest in each column
[y, ~] = max(a,[],dim);
dims = ones(1,ndims(a));
dims(dim) = size(a,dim);
a = a - repmat(y, dims);
s = y + log(sum(exp(a),dim));
i = find(~isfinite(y));
if ~isempty(i)
  s(i) = y(i);
end
end