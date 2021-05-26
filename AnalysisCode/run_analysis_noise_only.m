function [params_boot,sobl,abbl,best_hprs,trials,bin_centers, means, stderrs,actual_lo,actual_hi,data,sobl_time_locked] = run_analysis_noise_only(subjectID, expt_type, time, threshold, boot_n,threshold_trials,hpr1,hpr2,dir)
% [alpha bias kernel] = maxLikelihood('gaborV2-subject03', 'Contrast', [0 10], 20)

% initialize
datadir = fullfile(pwd, dir);
lo = threshold(1);
hi = threshold(2);

data = LoadAllSubjectData(subjectID,expt_type );
disp('Data loaded! Data: ');

signal_raw = [];
choice_raw = [];
noises = [];
if threshold_trials==0
    actual_lo = lo;
    actual_hi = hi;
else
    
    [floor, thresh] = getThresholdWrapper(subjectID);
    trials = sum(data.noise <= thresh & data.noise >= floor);
    actual_lo = floor;
    actual_hi = thresh;
end
for k = 1:size(data.ideal_frame_signals, 1)
    noises = [noises data.noise(k)];
    if data.noise(k) >= actual_lo && data.noise(k) <= actual_hi
        signal_raw = [signal_raw; data.ideal_frame_signals(k, :)];
        choice_raw = [choice_raw data.choice(k)];
    end
end

trials = size(choice_raw, 2);


[best_hprs, ~] = xValidatePK_with_lapse(signal_raw, choice_raw, hpr1, 0, hpr2, 0, 5);

disp('best params found')
for j = 1:boot_n
    [signal, choice] = bootstrap(signal_raw, choice_raw,trials);
    [sobl(j,:), ~] = LinearPK_with_lapse(signal, choice, 0);
    [sobl_time_locked(j,:), ~] = LinearPK_with_lapse_time_locked(signal, choice, time, 0);
    disp(j)
    [params_boot(j,:), ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(signal, choice, best_hprs(1), 0, best_hprs(3), 0);%,hprs, 0, hprs, 1);
    [abbl(j,:), ~, ~] = ExponentialPK_with_lapse(signal, choice, 0);
end

    function [signals, choices] = bootstrap(signals_raw, choices_raw,trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for i = 1:trials
            trial_num = sample_nums(i);
            signals = [signals; signals_raw(trial_num, :)];
            choices = [choices choices_raw(trial_num)];
        end
    end
    function [floor, thresh] = getThresholdWrapper(subjectID)
        window_low = 0.5;
        window_high = 0.7;
        [floor, thresh,~] = getThresholdWindow(data,subjectID, 2, window_low, window_high, datadir);
    end

bin_edges = linspace(min(mean(signal_raw,2)), max(mean(signal_raw,2)), 11);
bin_halfwidth = (bin_edges(2) - bin_edges(1)) / 2;
bin_centers = bin_edges(1:end-1) + bin_halfwidth;
means = zeros(size(bin_centers));
stderrs = zeros(size(bin_centers));
for b=1:length(bin_centers)
    % Select all points for which bin i is closest.
    bin_dists = abs(mean(signal_raw,2) - bin_centers(b));
    indices = bin_dists <= bin_halfwidth;
    means(b) = mean(data.accuracy(indices));
    stderrs(b) = std(data.accuracy(indices)) / sqrt(sum(indices));
end

end