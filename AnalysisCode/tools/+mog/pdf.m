function P = pdf(x, mog_p, discretize)
%MOG.PDF compute pdf of 1d mixture-of-gaussians.
%
%A MoG (in 1d only) is specified by a mean, standard deviation, and weight
%at each mode. A distribution with N modes is represented by a [1 x 3N] row vector:
%
%   mog = [mu_1 sigma_1 pi_1, ..., mu_n sigma_n, pi_n]
%
%To use use vectorization to evaluate many PDFs at once, mog can be a [M x 3N] matrix specifying M
%different MoG distributions. In this case, x should be [M x 1]
%
%Note: unexpected behavior may occur if any mode's sigma is too small, as in a delta distribution

P = exp(mog.logpdf(x, mog_p));
if nargin > 2 && discretize, P = P / sum(P); end
P = reshape(P, size(x));
end