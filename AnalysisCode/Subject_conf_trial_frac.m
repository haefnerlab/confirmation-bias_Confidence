load('/Users/rupamacharyya/Ankani_Projects/CB/confirmation-bias_Confidence/AnalysisCode/SavedWorkspace/ConfidenceAnalysis_21-Jun-2021.mat');
perf_lo = 0.3;
perf_hi = 0.7;
for i=1:10
    for j=1:2
        total_trials(i,j) = length(data_sub{i,j}.choice);
        rr(i,j) = sum(data_sub{i,j}.correct_answer)/length(data_sub{i,j}.ideal_answer);
        all_trials_conf_count(i,j,1) = sum(data_sub{i,j}.conf==1)/total_trials(i,j);
        all_trials_conf_count(i,j,2) = sum(data_sub{i,j}.conf==2)/total_trials(i,j);
        all_trials_conf_count(i,j,3) = sum(data_sub{i,j}.conf==3)/total_trials(i,j);
        [thresh_lo(i,j), thresh_hi(i,j), ~] = getThresholdWindow(data_sub{i,j}, 3-j, perf_lo, perf_hi);
        if j==1
            ind = data_sub{i,j}.sign_noise>=thresh_lo(i,j) & data_sub{i,j}.sign_noise<=thresh_hi(i,j);
            total_thresh_trials(i,j) = sum(ind);
            thresh_trials_conf_count(i,j,1) = sum(data_sub{i,j}.conf(ind)==1)/total_thresh_trials(i,j);
            thresh_trials_conf_count(i,j,2) = sum(data_sub{i,j}.conf(ind)==2)/total_thresh_trials(i,j);
            thresh_trials_conf_count(i,j,3) = sum(data_sub{i,j}.conf(ind)==3)/total_thresh_trials(i,j);
            unique_conf{i,j} = unique(data_sub{i,j}.conf(data_sub{i,j}.sign_noise>=thresh_lo(i,j) & data_sub{i,j}.sign_noise<=thresh_hi(i,j)));
        else
            ind = data_sub{i,j}.true_ratio>=thresh_lo(i,j) & data_sub{i,j}.true_ratio<=thresh_hi(i,j);
            total_thresh_trials(i,j) = sum(ind);
            thresh_trials_conf_count(i,j,1) = sum(data_sub{i,j}.conf(ind)==1)/total_thresh_trials(i,j);
            thresh_trials_conf_count(i,j,2) = sum(data_sub{i,j}.conf(ind)==2)/total_thresh_trials(i,j);
            thresh_trials_conf_count(i,j,3) = sum(data_sub{i,j}.conf(ind)==3)/total_thresh_trials(i,j);
            unique_conf{i,j} = unique(data_sub{i,j}.conf(data_sub{i,j}.true_ratio>=thresh_lo(i,j) & data_sub{i,j}.true_ratio<=thresh_hi(i,j)));
        end
    end
end

across_sub_conf(1,1) = mean(thresh_trials_conf_count(:,1,1));
across_sub_conf(1,2) = mean(thresh_trials_conf_count(:,1,2));
across_sub_conf(1,3) = mean(thresh_trials_conf_count(:,1,3));

across_sub_conf(2,1) = mean(thresh_trials_conf_count(:,2,1));
across_sub_conf(2,2) = mean(thresh_trials_conf_count(:,2,2));
across_sub_conf(2,3) = mean(thresh_trials_conf_count(:,2,3));