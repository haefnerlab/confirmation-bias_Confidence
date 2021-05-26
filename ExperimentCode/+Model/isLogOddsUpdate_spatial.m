function [lpo, x_samples, weights] = isLogOddsUpdate_spatial(params, e, lpo)
%MODEL.ISLOGODDSUPDATE compute update to log odds of C using importance sampling model.

trials = size(e, 1);
oz = ones(trials, 1);

updates = params.updates;
noise = params.noise;
gamma = params.gamma;
sig_s = sqrt(params.var_s);
sig_x = sqrt(params.var_x);
p_match = params.p_match;
samples = params.samples;
frames = params.frames;

% Create two distributions representing p(x|C=+1) and p(x|C=-1) in the generative model
p_x_Cp = mog.create([+1 -1], [sig_x sig_x], [p_match 1-p_match]);
p_x_Cm = mog.create([-1 +1], [sig_x sig_x], [p_match 1-p_match]);

x_samples = zeros(trials, samples, frames, updates);
weights = zeros(trials, samples, frames, updates);

for n=1:updates
    % Convert from lpo (log odds) to the probability that C=+1
    pC = 1 ./ (1 + exp(-lpo));
    
    % Create likelihoods. Format is a matrix where each row specifies triples of [mu, sigma, pi] of
    % a mixture of Gaussians. Only one mode in the likelihood, but it's useful to use the MoG
    % format. See @mog.create
    for fr=1:frames
        likelihoods(:,:,:, fr) = [e(:,fr) sig_s*oz oz];
    end
    
    % Create the prior on x by marginalizing over the current posterior of C. The prior is also a
    % mixture of gaussians, but with 2 modes corresponding to C = +/-1
    p = p_match * pC + (1 - p_match) * (1 - pC);
%     sig_x = sig_x + 1;
    priors = [+oz, sig_x*oz, p, -oz, sig_x*oz, 1-p];
        
    % Q is the distribution from which samples of x are drawn; it is the current estimate of the
    % posterior over x using lpo as the prior over C
    for fr=1:frames
        Q(:,:,fr) = mog.prod(squeeze(likelihoods(:,:,:,fr)), priors);
    end
    
    % Draw samples from Q
    for fr=1:frames
        samp = mog.sample(squeeze(Q(:,:,fr)), samples);
        x_samples(:, :, fr, n) = samp;
        flat_samples(fr,:) = samp(:);
    end
    
    % Get unnormalized importance weights for each sample (note that this is vectorized over trials,
    % but we have to loop over samples. Typically number of samples << number of trials)
    for fr=1:frames
        for s=1:samples
            weights(:, s, fr, n) = 1 ./ mog.pdf(x_samples(:, s, fr, n), priors);
        end
    end

    % Normalize importance-sampling weights
    for fr=1:frames
        if params.importance_norm
            weights(:, :, fr, n) = weights(:, :, fr, n) ./ sum(weights(:, :, fr, n), 2);
        end
    end

    % Compute p(x|C=+1) and p(x|C=-1), then take weighted sum for each trial.
    for fr=1:frames
        pCp(:,fr) = sum(reshape(mog.pdf(squeeze(flat_samples(fr,:))', p_x_Cp), [trials samples]) .* weights(:, :, fr, n), 2);
        pCm(:,fr) = sum(reshape(mog.pdf(squeeze(flat_samples(fr,:))', p_x_Cm), [trials samples]) .* weights(:, :, fr, n), 2);
    end
    % Log likelihood odds is log(pCp/pCm)
    for fr=1:frames
        llo(:,fr) = (log(pCp(:,fr)) - log(pCm(:,fr)));
    end
    
    lpo = lpo * (1 - gamma / updates) + sum(llo,2) / (updates);
    
    % Add zero-mean additive noise.
    lpo = lpo + randn(trials, 1) * noise;
end

end