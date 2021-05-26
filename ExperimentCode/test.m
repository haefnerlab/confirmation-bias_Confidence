% Range of values for category and sensory information
ps = 0.51:0.02:0.99;
THRESHOLD = 0.7;

% Sampling model with gamma = 0.1
params = Model.newModelParams('model', 'is', 'frames', 1,'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 1);
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)
disp('Loading/Running sampling model, getting threshold points for PKs');
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
hslc = sens_cat_pts(end,:);
lshc = sens_cat_pts(1,:);
sens_cat_pts_new = [hslc;lshc];
[correct, fig, fig1] = Model.plotCategorySensorySpace(ps, ps, params);
% [cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts_new);
% [lpo_fig, ~] = Model.plotCSPKLPO(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts_new);
% disp('Loading/Running sampling model, getting slopes over CS-Space');
% [~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts_new);


%%

% params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 1);
hslc = sens_cat_pts(end,:); %[0.9100    0.6300];
lshc = sens_cat_pts(1,:); %[0.7500    0.8500];
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
%%
results_lshc = Model.runVectorized(params_lshc);
results_hslc = Model.runVectorized(params_hslc);
[~, answer_lshc] = Model.genDataWithParams(params_lshc);
[~, answer_hslc] = Model.genDataWithParams(params_hslc);
perf_lshc = mean(results_lshc.choices == answer_lshc)
perf_hslc = mean(results_hslc.choices == answer_hslc)
%%
n_trajectories = 100;
ax1 = subplot(1,3,1);
plot((0:1),results_lshc.lpo(1:n_trajectories, :)', 'r');
hold on;
plot((0:1),zeros(1,1),'-.k','LineWidth',2);
xlabel('Frames');
ylabel('LPO');
title('LSHC Condition');

ax2 = subplot(1,3,2);
plot((0:1),results_hslc.lpo(1:n_trajectories, :)', 'b');
hold on;
plot((0:1),zeros(1,1),'-.k','LineWidth',2);
xlabel('Frames');
ylabel('LPO');
title('HSLC Condition');
[conf_lshc(1), conf_lshc(2), conf_lshc(3)] = meanci(abs(results_lshc.lpo(:, end)), .67);
[conf_hslc(1), conf_hslc(2), conf_hslc(3)] = meanci(abs(results_hslc.lpo(:, end)), .67);

ax3 = subplot(1,3,3); hold on;
b_lshc = bar([1], [conf_lshc(1)]);
b_hlsc = bar([2], [conf_hslc(1)]);
set(b_lshc,'FaceColor','r');
set(b_hlsc,'FaceColor','b');
errorbar([1 2], [conf_lshc(1) conf_hslc(1)], [conf_lshc(1)-conf_lshc(2) conf_hslc(1)-conf_hslc(2)], [conf_lshc(3)-conf_lshc(1) conf_hslc(3)-conf_hslc(1)], 'ok');
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel',{'LSHC' 'HSLC'});
ylabel('Confidence as absolute LPO')
title('Confidence Compare')

linkaxes([ax1, ax2], 'y')

%%

ps = 0.51:0.02:0.99;


