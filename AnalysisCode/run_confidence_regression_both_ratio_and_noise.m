function [conf_analysis,thresh,count] = run_confidence_regression_both_ratio_and_noise(subjectID,expt_type,perf_lo, perf_hi,boot_n,hpr_ridge,hpr_ar1,hpr_curvature,standardize,folds,dir,thresh_trials_only)

datadir = fullfile(pwd, dir);
data = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded for the subject! ');
data.conf = data.conf + 1; % because values are 0, 1 and 2 for pressing keys 1, 2 and 3. Convert confidence values to 1, 2 and 3.
signal_raw = [];
choice_raw = [];
conf_mid_low = [];
conf_high_mid = [];
conf_high_low = [];
conf_low_rest = [];
conf_mid_rest = [];
conf_high_rest = [];
ind_high_mid = [];
ind_high_low = [];
ind_mid_low = [];

if thresh_trials_only==1
[thresh_lo, thresh_hi, ~] = getThresholdWindow(data, expt_type, perf_lo, perf_hi);
else
    if expt_type==2
        thresh_lo = -0.8;
        thresh_hi = 0.8;
    elseif expt_type==1
        thresh_lo = 0.0;
        thresh_hi = 1.0;
    end
end
thresh = [thresh_lo thresh_hi];
count = 0;
ind_thresh = [];
for k = 1:size(data.ideal_frame_signals, 1)
    flag = 0;
    if expt_type==1
        if (data.true_ratio(k)>=thresh_lo && data.true_ratio(k)<=thresh_hi)
            flag = 1;
            count = count + 1;
        end
    elseif expt_type==2
        if (data.sign_noise(k)>=thresh_lo && data.sign_noise(k)<=thresh_hi)
            flag = 1;
            count = count + 1;
        end
    end
    
    signal_raw = [signal_raw; data.ideal_frame_signals(k, :)];
    choice_raw = [choice_raw data.choice(k)];
    if flag==1
        ind_thresh = [ind_thresh k];
        if data.conf(k)==1
            if ~isnan(data.conf(k))
                conf_low_rest = [conf_low_rest 1];
                conf_mid_rest = [conf_mid_rest 0];
                conf_high_rest = [conf_high_rest 0];
                
                conf_high_low = [conf_high_low 0];
                conf_mid_low = [conf_mid_low 0];
                ind_high_low = [ind_high_low k];
                ind_mid_low = [ind_mid_low k];
            end
        elseif data.conf(k)==2
            if ~isnan(data.conf(k))
                conf_mid_rest = [conf_mid_rest 1];
                conf_low_rest = [conf_low_rest 0];
                conf_high_rest = [conf_high_rest 0];
                
                conf_high_mid = [conf_high_mid 0];
                conf_mid_low = [conf_mid_low 1];
                ind_high_mid = [ind_high_mid k];
                ind_mid_low = [ind_mid_low k];
            end
        elseif data.conf(k)==3
            if ~isnan(data.conf(k))
                conf_high_rest = [conf_high_rest 1];
                conf_low_rest = [conf_low_rest 0];
                conf_mid_rest = [conf_mid_rest 0];
                
                conf_high_mid = [conf_high_mid 1];
                conf_high_low = [conf_high_low 1];
                ind_high_low = [ind_high_low k];
                ind_high_mid = [ind_high_mid k];
            end
        end
    end
end
disp('Data preprocessing complete!!');
disp(unique(data.conf));
trials = size(choice_raw, 2);
disp('Starting confidence kernel analysis ...');
for conf_compare_cases=1:6
    signal_chosen = [];
    reg_chosen = [];
    if conf_compare_cases==1
        conf_analysis{conf_compare_cases}.case='mid-low';
        signal_chosen = signal_raw(ind_mid_low,:);
        trials_chosen = size(conf_mid_low,2);
        conf_analysis{conf_compare_cases}.trials = trials_chosen;
        reg_chosen = conf_mid_low;
    elseif conf_compare_cases==2
        conf_analysis{conf_compare_cases}.case='high-mid';
        signal_chosen = signal_raw(ind_high_mid,:);
        trials_chosen = size(conf_high_mid,2);
        conf_analysis{conf_compare_cases}.trials = trials_chosen;
        reg_chosen = conf_high_mid;
    elseif conf_compare_cases==3
        conf_analysis{conf_compare_cases}.case='high-low';
        signal_chosen = signal_raw(ind_high_low,:);
        trials_chosen = size(conf_high_low,2);
        conf_analysis{conf_compare_cases}.trials = trials_chosen;
        reg_chosen = conf_high_low;
    elseif conf_compare_cases==4
        conf_analysis{conf_compare_cases}.case='low-rest';
        signal_chosen = signal_raw(~isnan(data.conf(ind_thresh)),:);
        trials_chosen = sum(~isnan(data.conf(ind_thresh)));
        conf_analysis{conf_compare_cases}.trials = trials_chosen;
        reg_chosen = conf_low_rest;
    elseif conf_compare_cases==5
        conf_analysis{conf_compare_cases}.case='mid-rest';
        signal_chosen = signal_raw(~isnan(data.conf(ind_thresh)),:);
        trials_chosen = sum(~isnan(data.conf(ind_thresh)));
        conf_analysis{conf_compare_cases}.trials = trials_chosen;
        reg_chosen = conf_mid_rest;
    elseif conf_compare_cases==6
        conf_analysis{conf_compare_cases}.case='high-rest';
        signal_chosen = signal_raw(~isnan(data.conf(ind_thresh)),:);
        trials_chosen = sum(~isnan(data.conf(ind_thresh)));
        conf_analysis{conf_compare_cases}.trials = trials_chosen;
        reg_chosen = conf_high_rest;
    end
    [best_hprs, ~] = CustomRegression.xValidatePKwithlapseSabya(signal_chosen, reg_chosen, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds);
    conf_analysis{conf_compare_cases}.hprs_used = best_hprs;
    disp(['Hyperparameters for Case ' num2str(conf_compare_cases)]);
    disp(best_hprs);
    for j = 1:boot_n
        if (j==1 || (mod(j,100)==0))
            disp(['Completed ' num2str(j) '/' num2str(boot_n) ' steps..']);
        end
        [signal, choice] = bootstrap(signal_chosen, reg_chosen, trials_chosen);
        if expt_type==1
            signal = sign(signal);
        end
        %         [conf_analysis{conf_compare_cases}.sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal, choice, standardize);
        %         [conf_analysis{conf_compare_cases}.abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal, choice, standardize);
        [conf_analysis{conf_compare_cases}.params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal, choice, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    end
    
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

end