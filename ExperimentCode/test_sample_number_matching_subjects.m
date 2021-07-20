subj_info = load('Subject_conf_trial_frac_20samples.mat');
% Range of values for category and sensory information
ps = 0.51:0.02:0.99;
THRESHOLD = 0.7;
pts = 50;
% Sampling model with gamma = 0.0
params = Model.newModelParams('frames',10, 'model', 'is_temporal', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 20);
%%
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)
disp('Loading/Running sampling model');
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
%%
% comparable 
% hslc = sens_cat_pts(end-37,:);
% lshc = sens_cat_pts(12,:);

hslc = sens_cat_pts(end,:);
lshc = sens_cat_pts(5,:);


params_hslc = params;
params_hslc.category_info = hslc(2);
params_hslc.sensory_info = hslc(1); 
params_hslc.p_match = hslc(2);
params_hslc.var_s = Model.getEvidenceVariance(hslc(1));
params_lshc = params;
params_lshc.category_info = lshc(2);
params_lshc.sensory_info = lshc(1);
params_lshc.p_match = lshc(2);
params_lshc.var_s = Model.getEvidenceVariance(lshc(1));


results_lshc = Model.runVectorized(params_lshc);
results_hslc = Model.runVectorized(params_hslc);
[~, answer_lshc] = Model.genDataWithParams(params_lshc);
[~, answer_hslc] = Model.genDataWithParams(params_hslc);
perf_lshc = mean(results_lshc.choices == answer_lshc)
perf_hslc = mean(results_hslc.choices == answer_hslc)



% test part
model_lshc_mean_lpo = mean(abs(results_lshc.lpo(:,end)))
model_hslc_mean_lpo = mean(abs(results_hslc.lpo(:,end)))
subj_lshc_mean_conf = subj_info.across_sub_conf(1,1)+subj_info.across_sub_conf(1,2)*2+subj_info.across_sub_conf(1,3)*3
subj_hslc_mean_conf = subj_info.across_sub_conf(2,1)+subj_info.across_sub_conf(2,2)*2+subj_info.across_sub_conf(2,3)*3
