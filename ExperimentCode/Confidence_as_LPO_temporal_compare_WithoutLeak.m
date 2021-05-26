% Range of values for category and sensory information
ps = 0.51:0.02:0.99;
THRESHOLD = 0.7;
pts = 4;
% Sampling model with gamma = 0.0
params = Model.newModelParams('frames',10, 'model', 'is_temporal', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
params_ideal = Model.newModelParams('frames',10, 'model', 'ideal', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
%%
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)
disp('Loading/Running sampling model');
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
hslc = sens_cat_pts(end,:);
lshc = sens_cat_pts(1,:);
sens_cat_pts_new = [hslc;lshc];
% [fig, fig1] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts_new);
[~, pc, conf] = Model.plotCategorySensorySpace(ps, ps, params); 
% [lpo_fig, ~] = Model.plotCSPKLPO(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts_new);
% disp('Loading/Running sampling model, getting slopes over CS-Space');
% [~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts_new);


disp('Loading/Running ideal model');
sens_cat_pts_ideal = Model.getThresholdPoints(ps, params_ideal, THRESHOLD, pts);
hslc_ideal = sens_cat_pts_ideal(end,:);
lshc_ideal = sens_cat_pts_ideal(1,:);
sens_cat_pts_new_ideal = [hslc_ideal;lshc_ideal];
% [cs_fig_ideal, pk_fig] = Model.plotCSPK(ps, ps, params_ideal, [0 0 0], 'beta', beta_range, sens_cat_pts_new_ideal);
[~, pc_ideal, conf_ideal] = Model.plotCategorySensorySpace(ps, ps, params_ideal);
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

% figure(pc);
% hold on;
% scatter(find(ps==hslc(1)), find(ps==hslc(2)), 200, 'b', 'filled');
% hold on;
% scatter(find(ps==lshc(1)), find(ps==lshc(2)), 200, 'r', 'filled');
% figure(conf);
% hold on;
% scatter(find(ps==hslc(1)), find(ps==hslc(2)), 200, 'b', 'filled');
% hold on;
% scatter(find(ps==lshc(1)), find(ps==lshc(2)), 200, 'r', 'filled');

hslc_ideal = sens_cat_pts_ideal(end,:); %[0.9100    0.6300];
lshc_ideal = sens_cat_pts_ideal(1,:); %[0.7500    0.8500];
params_hslc_ideal = params_ideal;
params_hslc_ideal.category_info = hslc_ideal(2);
params_hslc_ideal.sensory_info = hslc_ideal(1);
params_hslc_ideal.p_match = hslc_ideal(2);
params_hslc_ideal.var_s = Model.getEvidenceVariance(hslc_ideal(1));
params_lshc_ideal = params_ideal;
params_lshc_ideal.category_info = lshc_ideal(2);
params_lshc_ideal.sensory_info = lshc_ideal(1);
params_lshc_ideal.p_match = lshc_ideal(2);
params_lshc_ideal.var_s = Model.getEvidenceVariance(lshc_ideal(1));

% figure(pc_ideal);
% hold on;
% scatter(find(ps==hslc(1)), find(ps==hslc(2)), 200, 'b', 'filled');
% hold on;
% scatter(find(ps==lshc(1)), find(ps==lshc(2)), 200, 'r', 'filled');
% figure(conf_ideal);
% hold on;
% scatter(find(ps==hslc(1)), find(ps==hslc(2)), 200, 'b', 'filled');
% hold on;
% scatter(find(ps==lshc(1)), find(ps==lshc(2)), 200, 'r', 'filled');
%%
results_lshc = Model.runVectorized(params_lshc);
results_hslc = Model.runVectorized(params_hslc);
[~, answer_lshc] = Model.genDataWithParams(params_lshc);
[~, answer_hslc] = Model.genDataWithParams(params_hslc);
perf_lshc = mean(results_lshc.choices == answer_lshc)
perf_hslc = mean(results_hslc.choices == answer_hslc)


results_lshc_ideal = Model.runVectorized(params_lshc_ideal);
results_hslc_ideal = Model.runVectorized(params_hslc_ideal);
[~, answer_lshc_ideal] = Model.genDataWithParams(params_lshc_ideal);
[~, answer_hslc_ideal] = Model.genDataWithParams(params_hslc_ideal);
perf_lshc_ideal = mean(results_lshc_ideal.choices == answer_lshc_ideal)
perf_hslc_ideal = mean(results_hslc_ideal.choices == answer_hslc_ideal)

%%
[weights_hslc_ideal, errors_hslc_ideal] = Model.plotPK(params_hslc_ideal);
[weights_lshc_ideal, errors_lshc_ideal] = Model.plotPK(params_lshc_ideal);

[weights_lshc, errors_lshc] = Model.plotPK(params_lshc);
[weights_hslc, errors_hslc] = Model.plotPK(params_hslc);


%%
figure()
n_trajectories = 100;
ax1 = subplot(2,2,1);
plot((0:params.frames),results_lshc.lpo(1:n_trajectories, :)', 'r');
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',2);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Number of stimulus frames','FontSize', 20);
ylabel('Log Posterior Odds (LPO)','FontSize', 20);
title('LSHC Condition in model','FontSize', 20);
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
title('HSLC Condition in model','FontSize', 20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax3 = subplot(2,2,3);
plot((0:params_ideal.frames),results_lshc_ideal.lpo(1:n_trajectories, :)', 'color',[0.5 0 0]);
hold on;
plot((0:params_ideal.frames),zeros(1,params_ideal.frames+1),'-k','LineWidth',2);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Number of stimulus frames','FontSize', 20);
ylabel('Log Posterior Odds (LPO)','FontSize', 20);
title('LSHC Condition Ideal observer','FontSize', 20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax4 = subplot(2,2,4);
plot((0:params_ideal.frames),results_hslc_ideal.lpo(1:n_trajectories, :)', 'color',[0 0 0.5]);
hold on;
plot((0:params_ideal.frames),zeros(1,params_ideal.frames+1),'-k','LineWidth',1.5);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-.k','LineWidth',2);
xlabel('Number of stimulus frames','FontSize', 20);
ylabel('Log Posterior Odds (LPO)','FontSize', 20);
title('HSLC Condition Ideal observer','FontSize', 20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

linkaxes([ax1, ax2, ax3, ax4], 'y')

%%
load('ModelPredictionsWithoutLeak.mat')
wl = load('ModelPredictionsWithLeak.mat');
load('ModelPredictionsWithoutLeak.mat');

figure()
subplot(2,2,3)

conf_lshc = mean(abs(results_lshc.lpo(:, end)));
conf_hslc = mean(abs(results_hslc.lpo(:, end)));
err_lshc = sqrt(var(squeeze(abs(results_lshc.lpo(:, end))))/params_lshc.trials);
err_hslc = sqrt(var(squeeze(abs(results_hslc.lpo(:, end))))/params_hslc.trials);

wlconf_lshc = mean(abs(wl.results_lshc.lpo(:, end)));
wlconf_hslc = mean(abs(wl.results_hslc.lpo(:, end)));
wlerr_lshc = sqrt(var(squeeze(abs(wl.results_lshc.lpo(:, end))))/wl.params_lshc.trials);
wlerr_hslc = sqrt(var(squeeze(abs(wl.results_hslc.lpo(:, end))))/wl.params_hslc.trials);

conf_lshc_ideal = mean(abs(results_lshc_ideal.lpo(:, end)));
conf_hslc_ideal = mean(abs(results_hslc_ideal.lpo(:, end)));
err_lshc_ideal = sqrt(var(abs(results_lshc_ideal.lpo(:, end)))/params_lshc_ideal.trials);
err_hslc_ideal = sqrt(var(abs(results_hslc_ideal.lpo(:, end)))/params_hslc_ideal.trials);


model_series = [conf_lshc conf_hslc; wlconf_lshc wlconf_hslc; conf_lshc_ideal  conf_hslc_ideal];
model_error = [err_lshc err_hslc; wlerr_lshc wlerr_hslc; err_lshc_ideal  err_hslc_ideal];

labels = {{'Approx. Inf';'  model'}, {'Approx. Inf model'; '          with leak'},'Ideal observer'};

b = bar([1 2 3],model_series, 'grouped','FaceColor','flat');
hold on;

b(1).CData(1,:) = [1 0 0]; % group 1 1st bar
b(1).CData(2,:) = [0.8500 0.3250 0.0980]; % group 1 2nd bar
b(1).CData(3,:) = [0.5 0 0];

b(2).CData(1,:) = [0 0 1]; % group 1 1st bar
b(2).CData(2,:) = [0.5843 0.8157 0.9882]; % group 1 2nd bar
b(2).CData(3,:) = [0 0 0.5];

% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     text(x-0.075,model_series(:,i)+1,labels{i},'FontSize',15,'FontWeight','bold');
    errorbar(x, model_series(:,i), model_error(:,i),model_error(:,i),'ok', 'linestyle', 'none','linewidth', 2);
end
text([0.55 1.25 2.5],[6 4.5 1.5],labels,'FontSize',20,'FontWeight','bold');

ylabel({'Mean of absolute';' LPO as confidence'});
set(gca, 'XTick', [0.75 1.25 1.75 2.25])
set(gca, 'XTickLabel',{'LSHC' 'HSLC' 'LSHC' 'HSLC'});
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylim([0 7]);

subplot(2,2,4)
normalized = 1;
% LH(1)=errorbar(1:params_lshc.frames, weights_lshc(1:end-1), errors_lshc(1:end-1),'linewidth',5,'color','r');
% L{1} = 'LHSC';
% hold on;
% LH(2)=errorbar(1:params_hslc.frames, weights_hslc(1:end-1), errors_hslc(1:end-1),'linewidth',5,'color','b');
% L{2} = 'HSLC';
if normalized==1
    LH(1)=plot(1:params_lshc.frames, weights_lshc(1:end-1)/mean(weights_lshc(1:end-1)),'-o','linewidth',5,'color','r');
    L{1} = 'Approx. Inf model (LSHC)';
    hold on;
    LH(2)=plot(1:params_hslc_ideal.frames, weights_hslc(1:end-1)/mean(weights_hslc(1:end-1)),'-o','linewidth',5,'color','b');
    L{2} = 'Approx. Inf model (HSLC)';
    ylim([0.0 2]);
    ylabel('Temporal weights');
    hold on;
    LH(3)=plot(1:params_lshc.frames, wl.weights_lshc(1:end-1)/mean(wl.weights_lshc(1:end-1)),'-o','linewidth',5,'color',[0.8500, 0.3250, 0.0980]);
    L{3} = 'Approx. Inf model with leak (LSHC)';
    hold on;
    LH(4)=plot(1:params_hslc_ideal.frames, wl.weights_hslc(1:end-1)/mean(wl.weights_hslc(1:end-1)),'-o','linewidth',5,'color',[0.5843 0.8157 0.9882]);
    L{4} = 'Approx. Inf model with leak (HSLC)';
    ylim([0.0 2]);
    ylabel('Temporal weights');
    hold on;
else
    LH(1)=plot(1:params_lshc_ideal.frames, weights_lshc(1:end-1),'-o','linewidth',5,'color','r');
    L{1} = 'Model (LHSC)';
    hold on;
    LH(2)=plot(1:params_hslc_ideal.frames, weights_hslc(1:end-1),'-o','linewidth',5,'color','b');
    L{2} = 'Model (HSLC)';
    ylim([0. 10]);
    ylabel('Temporal weights');
end
xlabel('Stimulus frames');
% title('Model predictions','FontSize', 20)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
% legend(LH,L,'box','off','fontsize',20);
xlim([1 10]);


% subplot(2,2,4)
% LH(1)=errorbar(1:params_lshc_ideal.frames, weights_lshc_ideal(1:end-1), errors_lshc_ideal(1:end-1),'linewidth',5,'color','r');
% L{1} = 'LHSC';
% hold on;
% LH(2)=errorbar(1:params_hslc_ideal.frames, weights_hslc_ideal(1:end-1), errors_hslc_ideal(1:end-1),'linewidth',5,'color','b');
% L{2} = 'HSLC';
if normalized==1
    w1 = weights_lshc_ideal(1:end-1)/mean(weights_lshc_ideal(1:end-1));
    w2 = weights_hslc_ideal(1:end-1)/mean(weights_lshc_ideal(1:end-1));
%     LH(3)=plot(1:params_lshc_ideal.frames, weights_lshc_ideal(1:end-1)/mean(weights_lshc_ideal(1:end-1)),'--o','linewidth',5,'color',[0.5 0 0]);
%     L{3} = 'Ideal(LHSC)';
%     hold on;
%     LH(4)=plot(1:params_hslc_ideal.frames, weights_hslc_ideal(1:end-1)/mean(weights_lshc_ideal(1:end-1)),'--o','linewidth',5,'color',[0 0 0.5]);
%     L{4} = 'Ideal(HSLC)';
    
%     LH(3)=plot(1:params_lshc_ideal.frames, (w1+w2)/2,'--o','linewidth',5,'color',[0.5 0 0.5]);
%     L{3} = 'Ideal observer(both LHSC and HSLC)';
%     ylim([0.5 2.0]);
%     ylabel({'Temporal weights'});
    LH(5)=plot(1:params_lshc_ideal.frames, (w1+w2)/2,'--o','linewidth',5,'color',[0.5 0 0.5]);
    L{5} = 'Ideal observer (both LHSC and HSLC)';
    ylim([0.5 2.0]);
    ylabel({'Temporal weights'});
else
    LH(3)=plot(1:params_lshc_ideal.frames, weights_lshc_ideal(1:end-1),'-o','linewidth',5, 'color',[0.5 0 0]);
    L{3} = 'LHSC';
    hold on;
    LH(4)=plot(1:params_hslc_ideal.frames, weights_hslc_ideal(1:end-1),'-o','linewidth',5, 'color',[0 0 0.5]);
    L{4} = 'HSLC';
%     ylim([0.5 1.5]);
    ylabel('Temporal weights');
end
xlabel('Stimulus frames');

% title('Ideal Observer','FontSize', 20)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
legend(LH,L,'box','off','fontsize',20);
xlim([1 10]);
ylim([0.25 3]);


n_trajectories = 100;
ax1 = subplot(2,8,1);
plot((0:params.frames),results_lshc.lpo(1:n_trajectories, :)', 'r');
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Stimulus frames','FontSize', 20);
ylabel({'Log Posterior Odds';'(LPO)'},'FontSize', 20);
title('Model(LSHC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax2 = subplot(2,8,3);
plot((0:params.frames),results_hslc.lpo(1:n_trajectories, :)', 'b');
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-.k','LineWidth',2);
xlabel('Stimulus frames','FontSize', 20);
% ylabel('Log Posterior Odds (LPO)','FontSize', 20);
title('Model(HSLC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax3 = subplot(2,8,5);
plot((0:params_ideal.frames),results_lshc_ideal.lpo(1:n_trajectories, :)', 'color',[0.5 0 0]);
hold on;
plot((0:params_ideal.frames),zeros(1,params_ideal.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Stimulus frames','FontSize', 20);
% ylabel('Log Posterior Odds (LPO)','FontSize', 20);
title('Ideal(LHSC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax4 = subplot(2,8,7);
plot((0:params_ideal.frames),results_hslc_ideal.lpo(1:n_trajectories, :)', 'color',[0 0 0.5]);
hold on;
plot((0:params_ideal.frames),zeros(1,params_ideal.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-.k','LineWidth',2);
xlabel('Stimulus frames','FontSize', 20);
% ylabel('Log Posterior Odds (LPO)','FontSize', 20);
title('Ideal(HSLC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

linkaxes([ax1, ax2, ax3, ax4], 'y')


n_trajectories1 = 500;
bins1 = 100;
bins2 = 50;
subplot(2,8,2)
histogram(results_lshc.lpo(1:n_trajectories1, end),bins1,'Normalization','probability','FaceColor',[1 0 0],'Orientation','Horizontal','EdgeColor', [1 0 0],'FaceAlpha',0.25);
ylim([-10 10]);
ax = gca;
ax.Visible = 'off';

subplot(2,8,4)
histogram(results_hslc.lpo(1:n_trajectories1, end),bins2,'Normalization','probability','FaceColor',[0 0 1],'Orientation','Horizontal','EdgeColor', [0 0 1],'FaceAlpha',0.25);
ylim([-10 10]);
ax = gca;
ax.Visible = 'off';

subplot(2,8,6)
histogram(results_lshc_ideal.lpo(1:n_trajectories1, end),bins2,'Normalization','probability','FaceColor',[0.5 0 0],'Orientation','Horizontal','EdgeColor', [0.5 0 0],'FaceAlpha',0.25);
ylim([-10 10]);
ax = gca;
ax.Visible = 'off';

subplot(2,8,8)
histogram(results_hslc_ideal.lpo(1:n_trajectories1, end),bins2,'Normalization','probability','FaceColor',[0 0 0.5],'Orientation','Horizontal','EdgeColor', [0 0 0.5],'FaceAlpha',0.25);
ylim([-10 10]);
ax = gca;
ax.Visible = 'off';


%%
figure()

n_trajectories = 100;
ax1 = subplot(2,4,1);
plot((0:params.frames),results_lshc.lpo(1:n_trajectories, :)', 'r','LineWidth',0.01);
hold on;
scatter(ones(1,n_trajectories)*(params_ideal.frames+.25), results_lshc.lpo(1:n_trajectories, end),20,'r','filled','MarkerFaceAlpha',0.25 );
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Stimulus frames','FontSize', 20);
ylabel({'Log Posterior Odds';'(LPO)'},'FontSize', 20);
title('Model(LSHC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax2 = subplot(2,4,2);
plot((0:params.frames),results_hslc.lpo(1:n_trajectories, :)', 'b','LineWidth',0.01);
hold on;
scatter(ones(1,n_trajectories)*(params_ideal.frames+.25), results_hslc.lpo(1:n_trajectories, end),20,'b','filled','MarkerFaceAlpha',0.25 );
hold on;
plot((0:params.frames),zeros(1,params.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-.k','LineWidth',2);
xlabel('Stimulus frames','FontSize', 20);
ylabel({'Log Posterior Odds';'(LPO)'},'FontSize', 20);
title('Model(HSLC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax3 = subplot(2,4,3);
plot((0:params_ideal.frames),results_lshc_ideal.lpo(1:n_trajectories, :)', 'color',[0.5 0 0],'LineWidth',0.01);
hold on;
scatter(ones(1,n_trajectories)*(params_ideal.frames+.25), results_lshc_ideal.lpo(1:n_trajectories, end),20,[0.5 0 0],'filled','MarkerFaceAlpha',0.25 );
hold on;
plot((0:params_ideal.frames),zeros(1,params_ideal.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-k','LineWidth',1.5);
xlabel('Stimulus frames','FontSize', 20);
ylabel({'Log Posterior Odds';'(LPO)'},'FontSize', 20);
title('Ideal(LHSC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

ax4 = subplot(2,4,4);
plot((0:params_ideal.frames),results_hslc_ideal.lpo(1:n_trajectories, :)', 'color',[0 0 0.5],'LineWidth',0.01);
hold on;
scatter(ones(1,n_trajectories)*(params_ideal.frames+.25), results_hslc_ideal.lpo(1:n_trajectories, end),20,[0 0 0.5],'filled','MarkerFaceAlpha',0.25 );
plot((0:params_ideal.frames),zeros(1,params_ideal.frames+1),'-k','LineWidth',3);
hold on;
% plot(ones(1,params.frames+1),linspace(-30,30,params.frames+1),'-.k','LineWidth',2);
xlabel('Stimulus frames','FontSize', 20);
ylabel({'Log Posterior Odds';'(LPO)'},'FontSize', 20);
title('Ideal(HSLC)','FontSize', 15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

linkaxes([ax1, ax2, ax3, ax4], 'y')


n_trajectories1 = 10000;

subplot(2,4,5)
[f,xi] = ksdensity(results_lshc.lpo(1:n_trajectories1, end)) ;
bins1 = 30;%ceil((max(results_lshc.lpo(1:n_trajectories1, end))-min(results_lshc.lpo(1:n_trajectories1, end)))/bw);%50;
bins2 = 20;
histogram(results_lshc.lpo(1:n_trajectories1, end),29,'Normalization','probability','FaceColor',[1 0 0],'EdgeColor', 'none','FaceAlpha',0.5);%,'Orientation','Horizontal');
hold on;
plot(xi,f,'color',[1 0 0],'linewidth',2)
% set(gca,'YDir','reverse');
xlim([-15 15]);
ax = gca;
ax.Visible = 'off';
camroll(-90);

subplot(2,4,6)
[f,xi] = ksdensity(results_hslc.lpo(1:n_trajectories1, end)) ;
histogram(results_hslc.lpo(1:n_trajectories1, end),8,'Normalization','probability','FaceColor',[0 0 1],'EdgeColor', 'none','FaceAlpha',0.5);
hold on;
plot(xi,f,'color',[0 0 1],'linewidth',2)
xlim([-15 15]);
ax = gca;
ax.Visible = 'off';
camroll(-90);

subplot(2,4,7)
[f,xi] = ksdensity(results_lshc_ideal.lpo(1:n_trajectories1, end)) ;
histogram(results_lshc_ideal.lpo(1:n_trajectories1, end),9,'Normalization','probability','FaceColor',[0.5 0 0],'EdgeColor', 'none','FaceAlpha',0.5);
hold on;
plot(xi,f,'color',[0.5 0 0],'linewidth',2)
xlim([-15 15]);
ax = gca;
ax.Visible = 'off';
camroll(-90);

subplot(2,4,8)
[f,xi,bw] = ksdensity(results_hslc_ideal.lpo(1:n_trajectories1, end)) ;
histogram(results_hslc_ideal.lpo(1:n_trajectories1, end),9,'Normalization','probability','FaceColor',[0 0 0.5],'EdgeColor', 'none','FaceAlpha',0.5);
hold on;
plot(xi,f,'color',[0 0 0.5],'linewidth',2)
xlim([-15 15]);
ax = gca;
ax.Visible = 'off';
camroll(-90);