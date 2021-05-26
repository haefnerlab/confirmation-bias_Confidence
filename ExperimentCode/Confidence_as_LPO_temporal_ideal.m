% Range of values for category and sensory information
ps = 0.51:0.02:0.99;
THRESHOLD = 0.7;

% Sampling model with gamma = 0.1
params = Model.newModelParams('frames',10, 'model', 'ideal', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
%%
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)
disp('Loading/Running sampling model, getting threshold points for PKs');
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, 5);
hslc = sens_cat_pts(end,:);
lshc = sens_cat_pts(1,:);
sens_cat_pts_new = [hslc;lshc];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts_new);
% [lpo_fig, ~] = Model.plotCSPKLPO(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts_new);
% disp('Loading/Running sampling model, getting slopes over CS-Space');
% [~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts_new);


%%
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
ax1 = subplot(2,2,1);
plot((0:params.frames),results_lshc.lpo(1:n_trajectories, :)', 'r');
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',2);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Number of stimulus frames','FontSize', 20);
ylabel('Log Posterior Odds (LPO)','FontSize', 20);
% title('LSHC Condition');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax2 = subplot(2,2,2);
plot((0:params.frames),results_hslc.lpo(1:n_trajectories, :)', 'b');
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',1.5);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-.k','LineWidth',2);
xlabel('Number of stimulus frames','FontSize', 20);
ylabel('Log Posterior Odds (LPO)','FontSize', 20);
% title('HSLC Condition');
[conf_lshc(1), conf_lshc(2), conf_lshc(3)] = meanci(abs(results_lshc.lpo(:, end)), .67);
[conf_hslc(1), conf_hslc(2), conf_hslc(3)] = meanci(abs(results_hslc.lpo(:, end)), .67);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax3 = subplot(2,2,3); hold on;
b_lshc = bar([1], [conf_lshc(1)],'LineWidth',0.75);
b_hlsc = bar([2], [conf_hslc(1)],'LineWidth',0.75);
set(b_lshc,'FaceColor','r');
set(b_hlsc,'FaceColor','b');
errorbar([1 2], [conf_lshc(1) conf_hslc(1)], [conf_lshc(1)-conf_lshc(2) conf_hslc(1)-conf_hslc(2)], [conf_lshc(3)-conf_lshc(1) conf_hslc(3)-conf_hslc(1)], 'ok','LineWidth', 1.5);
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel',{'LSHC' 'HSLC'});
ylabel('Confidence as absolute LPO')
% title('Confidence Compare')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

linkaxes([ax1, ax2], 'y')

% temp_samples_lshc = reshape(results_lshc.samples,params.trials,params.samples*params.frames*params.updates);
% temp_samples_hslc = reshape(results_hslc.samples,params.trials,params.samples*params.frames*params.updates);
% [cer_lshc(1), cer_lshc(2), cer_lshc(3)] = meanci(var(temp_samples_lshc,0,2),0.67);
% [cer_hslc(1), cer_hslc(2), cer_hslc(3)] = meanci(var(temp_samples_hslc,0,2),0.67);
% ax4 = subplot(1,4,4); hold on;
% b_lshc = bar([1], [cer_lshc(1)]);
% b_hlsc = bar([2], [cer_hslc(1)]);
% set(b_lshc,'FaceColor','r');
% set(b_hlsc,'FaceColor','b');
% errorbar([1 2], [cer_lshc(1) cer_hslc(1)], [cer_lshc(1)-cer_lshc(2) cer_hslc(1)-cer_hslc(2)], [cer_lshc(3)-cer_lshc(1) cer_hslc(3)-cer_hslc(1)], 'ok');
% set(gca, 'XTick', [1 2])
% set(gca, 'XTickLabel',{'LSHC' 'HSLC'});
% ylabel('Uncertainty as absolute Sample Variance')
% title('Uncertainty Compare')


