function [params_boot,sobl,abbl,best_hprs,trials,bin_centers,means,stderrs,...
    confidence,conf_th,conf_trials,subj_resp,subj_resp_err,confidence_resp,confidence_resp_err,ntrial_subj,...
    conf_hist, conf_rt,choice_rt,direct_acc_conf, data, acc_m_l, acc_h_m, log_bernoulli] = run_analysis_both_ratio_and_noise(subjectID, expt_type,boot_n,best_hprs,standardize,dir)

% initialize
datadir = fullfile(pwd, dir);
data = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded for the subject! ');
data.conf = data.conf + 1; % because values are 0, 1 and 2 for pressing keys 1, 2 and 3. Convert confidence values to 1, 2 and 3.
signal_raw = [];
choice_raw = [];
for k = 1:size(data.ideal_frame_signals, 1)
    signal_raw = [signal_raw; data.ideal_frame_signals(k, :)];
    choice_raw = [choice_raw data.choice(k)];
end
disp('Data preprocessing complete!!');
disp(max(data.conf));
disp(unique(data.conf));
trials = size(choice_raw, 2);
[confidence, conf_th, conf_trials, acc_m_l, acc_h_m] = ConfidenceCompute(data,expt_type);
[subj_resp,subj_resp_err,confidence_resp,confidence_resp_err,ntrial_subj,conf_hist, conf_rt, choice_rt,direct_acc_conf] = computeConfidenceStatistics(data,expt_type);
disp('Confidence analysis complete!!');
disp('Starting kernel analysis ...');
for j = 1:boot_n
    if (j==1 || (mod(j,100)==0))
        disp(['Completed ' num2str(j) '/' num2str(boot_n) ' steps..']);
    end
    [signal, choice] = bootstrap(signal_raw, choice_raw, trials);
    if expt_type==1
        signal = sign(signal);
    end
    [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal, choice, standardize);
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal, choice, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal, choice, standardize);
end

all_frames = size(signal_raw,2);
temporal_kernel = prctile(params_boot(:, 1:all_frames), 50);
bias =  prctile(params_boot(:, end-1), 50);
disp('Getting log odds...');
[log_bernoulli] = compute_log_odds(signal_raw, temporal_kernel, bias);

    function [logits] = compute_log_odds(data,weights,bias_computed)
        logits = data * weights(:) + bias_computed;
    end

    function [signal, choice] = bootstrap(signal_raw, choice_raw,trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signal = [];
        choice = [];
        for i = 1:trials
            trial_num = sample_nums(i);
            signal = [signal; signal_raw(trial_num, :)];
            choice = [choice choice_raw(trial_num)];
        end
    end


    function [confidence_mean, thresh_return, conf_tr, acc_m_l, acc_h_m] = ConfidenceCompute(d,exp)
        [flr, th] = getThresholdWrapper(d, exp);
        if exp==1
            thresh_return = 0.5 + (abs(th-0.5) + abs(0.5-flr))/2;
        elseif exp==2
            thresh_return = 0.0 + (abs(th-0.0) + abs(0.0-flr))/2;
        end
        trl = size(d.ideal_frame_signals, 1);
        noise = [];
        true_ratio = [];
        conf = [];
        if exp==2
            for i = 1:trl
                if d.sign_noise(i)<=th && d.sign_noise(i)>=flr
                    noise = [noise; d.noise(i)];
                    conf = [conf; d.conf(i)];
                end
            end
            conf_tr = length(noise);
        elseif exp==1
            for i = 1:trl
                if (d.true_ratio(i)<=th && d.true_ratio(i)>=flr)
                    true_ratio = [true_ratio; d.true_ratio(i)];
                    conf = [conf; d.conf(i)];
                end
            end
            conf_tr = length(true_ratio);
        end
        for bt=1:10000
            sample_nums = randsample(conf_tr, conf_tr, true);
            conf_temp = conf(sample_nums);
            confidence_mean(bt) = mean(conf_temp);
            temp_acc = d.accuracy(sample_nums);
            acc_m_l(bt) = sum(temp_acc(conf(sample_nums)==2))/sum((conf(sample_nums)==2)) - sum(temp_acc(conf(sample_nums)==1))/sum((conf(sample_nums)==1));
            acc_h_m(bt) = sum(temp_acc(conf(sample_nums)==3))/sum((conf(sample_nums)==3)) - sum(temp_acc(conf(sample_nums)==2))/sum((conf(sample_nums)==2));
        end
    end


    function [floor, thresh] = getThresholdWrapper(data_sub,phase)
        perf_low = 0.3;
        perf_high = 0.7;
        if phase==1
            phase_thresh = 1;
        else
            phase_thresh = -2;
        end
        [fit_result, ~, ~, ~] = GaborPsychometric(data_sub, phase_thresh);
        floor = getThreshold(fit_result, perf_low, false);
        thresh = getThreshold(fit_result, perf_high, false);
%         [floor, thresh,~] = getThresholdWindow(data_sub, phase, window_low, window_high);
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