clear all; close all; clc;

load('Model_test_20samples.mat');
subj_info = load('Subject_conf_trial_frac.mat');

results_lshc = Model.runVectorized(params_lshc);
results_hslc = Model.runVectorized(params_hslc);
% lower = min(min(abs(results_lshc.lpo(:,end))),min(abs(results_hslc.lpo(:,end))));
% upper = max(max(abs(results_lshc.lpo(:,end))),max(abs(results_hslc.lpo(:,end))));
% bn_edges = linspace(lower, upper,4);
bn_edges_lshc(1) = 0; bn_edges_lshc(2) = prctile(abs(results_lshc.lpo(:,end)),subj_info.across_sub_conf(1,1)*100);
bn_edges_lshc(3) = prctile(abs(results_lshc.lpo(:,end)),sum(subj_info.across_sub_conf(1,1:2))*100);
bn_edges_lshc(4) = max(abs(results_lshc.lpo(:,end)));
bn_edges_hslc(1) = 0; bn_edges_hslc(2) = prctile(abs(results_hslc.lpo(:,end)),subj_info.across_sub_conf(2,1)*100);
bn_edges_hslc(3) = prctile(abs(results_hslc.lpo(:,end)),sum(subj_info.across_sub_conf(2,1:2))*100);
bn_edges_hslc(4) = max(abs(results_hslc.lpo(:,end)));

results_lshc.conf = zeros(1,results_lshc.params.trials);
results_hslc.conf = zeros(1,results_hslc.params.trials);

results_lshc.conf(find(abs(results_lshc.lpo(:,end))>=bn_edges_lshc(1) & abs(results_lshc.lpo(:,end))<=bn_edges_lshc(2))') = 1;
results_lshc.conf(find(abs(results_lshc.lpo(:,end))>=bn_edges_lshc(2) & abs(results_lshc.lpo(:,end))<=bn_edges_lshc(3))') = 2;
results_lshc.conf(find(abs(results_lshc.lpo(:,end))>=bn_edges_lshc(3) & abs(results_lshc.lpo(:,end))<=bn_edges_lshc(4))') = 3;
results_hslc.conf(find(abs(results_hslc.lpo(:,end))>=bn_edges_hslc(1) & abs(results_hslc.lpo(:,end))<=bn_edges_hslc(2))') = 1;
results_hslc.conf(find(abs(results_hslc.lpo(:,end))>=bn_edges_hslc(2) & abs(results_hslc.lpo(:,end))<=bn_edges_hslc(3))') = 2;
results_hslc.conf(find(abs(results_hslc.lpo(:,end))>=bn_edges_hslc(3) & abs(results_hslc.lpo(:,end))<=bn_edges_hslc(4))') = 3;

boots = 500;% number of bootstraps to get PK
hpr_ridge = logspace(-1, 5, 7);
hpr_ar1 = 0.0;
hpr_curvature = logspace(-1, 5, 7);
expt_type_noise = 2;
expt_type_ratio = 1;
standardize = 0; % z-score data or not
folds = 10; % how many folds of cross validation
num_frames = 10; % number of stimulus frames
cases = 2; % ratio case and noise case, so 2 cases. phase 2 for noise and 1 ratio

%%
model_lshc_mean_lpo = mean(abs(results_lshc.lpo(:,end)))
model_hslc_mean_lpo = mean(abs(results_hslc.lpo(:,end)))
subj_lshc_mean_conf = subj_info.across_sub_conf(1,1)+subj_info.across_sub_conf(1,2)*2+subj_info.across_sub_conf(1,3)*3
subj_hslc_mean_conf = subj_info.across_sub_conf(2,1)+subj_info.across_sub_conf(2,2)*2+subj_info.across_sub_conf(2,3)*3
%%
for j=1:cases
    phase = 3-j;
    if phase==expt_type_noise
        disp('Starting analysis for NOISE data of model');
        cs_txt = 'noise';
        results = results_lshc;
    elseif phase==expt_type_ratio
        disp('Starting analysis for RATIO data of model');
        cs_txt = 'ratio';
        results = results_hslc;
    end
    [conf_analysis{j}] = run_confidence_regression_model(results,boots,phase,hpr_ridge,hpr_ar1,hpr_curvature,standardize,folds);
    for cs=1:6
        temporal_kernel(cs,j,:) = prctile(conf_analysis{j}{cs}.params_boot(:, 1:num_frames), 50);
        norm_temporal_kernel(cs,j,:) = temporal_kernel(cs,j,:)/mean(temporal_kernel(cs,j,:));
    end
    disp(['All analysis complete for case ' cs_txt]);
    disp('-----------------------------------------------------------------------------------------------------');
end
%%
f = figure();
set(f,'defaultLegendAutoUpdate','off');
for cs=1:6
    subplot(2,3,cs);
    plot(1:num_frames,squeeze(norm_temporal_kernel(cs,1,1:num_frames)),'-or','LineWidth',2);
    cs1_sub = squeeze(norm_temporal_kernel(cs,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(cs,2,1:num_frames)),'-ob','LineWidth',2);
    cs2_sub = squeeze(norm_temporal_kernel(cs,2,1:num_frames));
    legend({['LSHC'],['HSLC']},'box','off');
    hold('on');
    axis('tight');
    yline(0.0,'k','linewidth',2);
    if cs<=3
        ylim([-2 4]);
    else
        ylim([-8 10]);
    end
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title([conf_analysis{1}{cs}.case ' comparison'],'fontsize',20);
end
%%
for j=1:cases
    phase = 3-j;
    if phase==expt_type_noise
        disp('Starting analysis for NOISE data of model');
        cs_txt = 'noise';
        results = results_lshc;
    elseif phase==expt_type_ratio
        disp('Starting analysis for RATIO data of model');
        cs_txt = 'ratio';
        results = results_hslc;
    end
    for cs=1:6
        lo_temporal_kernel(cs,j,:) = (prctile(conf_analysis{j}{cs}.params_boot(:, 1:num_frames), 50) - prctile(conf_analysis{j}{cs}.params_boot(:, 1:num_frames), 16))/(mean(temporal_kernel(cs,j,:)));
        hi_temporal_kernel(cs,j,:) = (prctile(conf_analysis{j}{cs}.params_boot(:, 1:num_frames), 84) - prctile(conf_analysis{j}{cs}.params_boot(:, 1:num_frames), 50))/(mean(temporal_kernel(cs,j,:)));
        norm_temporal_kernel(cs,j,:) = temporal_kernel(cs,j,:)/mean(temporal_kernel(cs,j,:));
    end
    disp(['All analysis complete for case ' cs_txt]);
    disp('-----------------------------------------------------------------------------------------------------');
end
f4 = figure();
set(f4,'defaultLegendAutoUpdate','off');
for cs=1:6
    subplot(2,3,cs);
    errorbar(1:num_frames,squeeze(norm_temporal_kernel(cs,1,1:num_frames)),...
        squeeze(lo_temporal_kernel(cs,1,1:num_frames)),squeeze(hi_temporal_kernel(cs,1,1:num_frames)),'-or','LineWidth',2);
    cs1_sub = squeeze(norm_temporal_kernel(cs,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    errorbar(1:num_frames,squeeze(norm_temporal_kernel(cs,2,1:num_frames)),...
        squeeze(lo_temporal_kernel(cs,2,1:num_frames)),squeeze(hi_temporal_kernel(cs,2,1:num_frames)),'-ob','LineWidth',2);
    cs2_sub = squeeze(norm_temporal_kernel(cs,2,1:num_frames));
    legend({['LSHC'],['HSLC']},'box','off');
    hold('on');
    axis('tight');
    yline(0.0,'k','linewidth',2);
    if cs<=3
        ylim([-2 4]);
    else
        ylim([-8 10]);
    end
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title([conf_analysis{1}{cs}.case ' comparison'],'fontsize',20);
end
%%
f1 = figure();
set(f1,'defaultLegendAutoUpdate','off');

[f_lshc,x_lshc] = ksdensity(abs(results_lshc.lpo(:,end)),'support','positive');
[f_hslc,x_hslc] = ksdensity(abs(results_hslc.lpo(:,end)),'support','positive');
LH1(1) = plot(x_lshc,f_lshc,'r','linewidth',2);
L1{1} = 'LSHC';
hold on;
xline(bn_edges_lshc(2),'--r','linewidth',1);xline(bn_edges_lshc(3),'--r','linewidth',1);
LH1(2) = plot(x_hslc,f_hslc,'b','linewidth',2);
L1{2} = 'HSLC';
hold on;
xline(bn_edges_hslc(2),'--b','linewidth',1);xline(bn_edges_hslc(3),'--b','linewidth',1);
legend(LH1,L1,'box','off','fontsize',15);
xlim([0,20]);
ax = gca;
set(ax, 'box','off');
ax.LineWidth=2;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xlabel('absolute LPO');
ylabel('density (ksdensity)');
title('Conf in model w.r.t avg subject','fontsize',20);

%%
bin_num = 21;
bins_sig_lshc = linspace(min(mean(results_lshc.data,2)),max(mean(results_lshc.data,2)),bin_num);
bins_sig_hslc = linspace(min(mean(results_hslc.data,2)),max(mean(results_hslc.data,2)),bin_num);

for i=1:bin_num-1
    
    bin_mid_lshc(i) = (bins_sig_lshc(i) + bins_sig_lshc(i+1))/2;
    bin_mid_hslc(i) = (bins_sig_hslc(i) + bins_sig_hslc(i+1))/2;
    indx_lshc{i} = find(mean(results_lshc.data,2)>=bins_sig_lshc(i) & mean(results_lshc.data,2)<=bins_sig_lshc(i+1));
    indx_hslc{i} = find(mean(results_hslc.data,2)>=bins_sig_hslc(i) & mean(results_hslc.data,2)<=bins_sig_hslc(i+1));
    
    conf_binned_lshc(1,i) = sum(results_lshc.conf(indx_lshc{i})==1)/length(indx_lshc{i});
    err_conf_binned_lshc(1,i) = sqrt((conf_binned_lshc(1,i) * (1-conf_binned_lshc(1,i)))/(length(indx_lshc{i})));
    conf_binned_lshc(2,i) = sum(results_lshc.conf(indx_lshc{i})==2)/length(indx_lshc{i});
    err_conf_binned_lshc(2,i) = sqrt((conf_binned_lshc(2,i) * (1-conf_binned_lshc(2,i)))/(length(indx_lshc{i})));
    conf_binned_lshc(3,i) = sum(results_lshc.conf(indx_lshc{i})==3)/length(indx_lshc{i});
    err_conf_binned_lshc(3,i) = sqrt((conf_binned_lshc(3,i) * (1-conf_binned_lshc(3,i)))/(length(indx_lshc{i})));
    
    conf_binned_hslc(1,i) = sum(results_hslc.conf(indx_hslc{i})==1)/length(indx_hslc{i});
    err_conf_binned_hslc(1,i) = sqrt((conf_binned_hslc(1,i) * (1-conf_binned_hslc(1,i)))/(length(indx_hslc{i})));
    conf_binned_hslc(2,i) = sum(results_hslc.conf(indx_hslc{i})==2)/length(indx_hslc{i});
    err_conf_binned_hslc(2,i) = sqrt((conf_binned_hslc(2,i) * (1-conf_binned_hslc(2,i)))/(length(indx_hslc{i})));
    conf_binned_hslc(3,i) = sum(results_hslc.conf(indx_hslc{i})==3)/length(indx_hslc{i});
    err_conf_binned_hslc(3,i) = sqrt((conf_binned_hslc(3,i) * (1-conf_binned_hslc(3,i)))/(length(indx_hslc{i})));
    
end

clr = {'r', 'g', 'b'};
% clr = {'m', 'c', 'g'};
figure();
subplot(1,2,1)
for cc=1:3
    %     LH2(cc) = plot(bin_mid_lshc,conf_binned_lshc(cc,:),['-o' clr{cc}],'linewidth',2);
    LH2(cc) = errorbar(bin_mid_lshc,conf_binned_lshc(cc,:),err_conf_binned_lshc(cc,:),['-o' clr{cc}],'linewidth',2);
    L2{cc} = ['conf=' num2str(cc)];
    hold on;
end
legend(LH2,L2,'box','off','fontsize',15);
ylabel('Fraction reported conf','fontsize',20);
xlabel('Signal','fontsize',20);
title('LSHC','fontsize',20)
ax = gca;
set(ax, 'box','off');
ax.LineWidth=2;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,2,2)
for cc=1:3
%     LH3(cc) = plot(bin_mid_hslc,conf_binned_hslc(cc,:),['-o' clr{cc}],'linewidth',2);
    LH3(cc) = errorbar(bin_mid_hslc,conf_binned_hslc(cc,:),err_conf_binned_hslc(cc,:),['-o' clr{cc}],'linewidth',2);
    L3{cc} = ['conf=' num2str(cc)];
    hold on;
end
legend(LH3,L3,'box','off','fontsize',15);
ylabel('Fraction reported conf','fontsize',20);
xlabel('Signal','fontsize',20);
title('HSLC','fontsize',20)
ax = gca;
set(ax, 'box','off');
ax.LineWidth=2;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

%%
bin_num_cic = 11;
bins_sig_lshc_cic = linspace(min(mean(results_lshc.data,2)),max(mean(results_lshc.data,2)),bin_num_cic);
bins_sig_hslc_cic = linspace(min(mean(results_hslc.data,2)),max(mean(results_hslc.data,2)),bin_num_cic);
for i=1:bin_num_cic-1
    bin_mid_lshc_cic(i) = (bins_sig_lshc_cic(i) + bins_sig_lshc_cic(i+1))/2;
    bin_mid_hslc_cic(i) = (bins_sig_hslc_cic(i) + bins_sig_hslc_cic(i+1))/2;
    indx_lshc_cic{i} = find(mean(results_lshc.data,2)>=bins_sig_lshc_cic(i) & mean(results_lshc.data,2)<=bins_sig_lshc_cic(i+1));
    indx_hslc_cic{i} = find(mean(results_hslc.data,2)>=bins_sig_hslc_cic(i) & mean(results_hslc.data,2)<=bins_sig_hslc_cic(i+1));
    
    indx_lshc_correct_cic{i} = indx_lshc_cic{i}(results_lshc.choices(indx_lshc_cic{i})==results_lshc.ideal_choices(indx_lshc_cic{i}));
    indx_hslc_correct_cic{i} = indx_hslc_cic{i}(results_hslc.choices(indx_hslc_cic{i})==results_hslc.ideal_choices(indx_hslc_cic{i}));
    
    indx_lshc_incorrect_cic{i} = setdiff(indx_lshc_cic{i},indx_lshc_correct_cic{i});
    indx_hslc_incorrect_cic{i} = setdiff(indx_hslc_cic{i},indx_hslc_correct_cic{i});
    
    
    conf_binned_lshc_correct(i) = mean(results_lshc.conf(indx_lshc_correct_cic{i}));
    err_conf_binned_lshc_correct(i) = sqrt(var(results_lshc.conf(indx_lshc_correct_cic{i}))/length(indx_lshc_correct_cic{i}));
    conf_binned_lshc_incorrect(i) = mean(results_lshc.conf(indx_lshc_incorrect_cic{i}));
    err_conf_binned_lshc_incorrect(i) = sqrt(var(results_lshc.conf(indx_lshc_incorrect_cic{i}))/length(indx_lshc_incorrect_cic{i}));
    
    conf_binned_hslc_correct(i) = mean(results_hslc.conf(indx_hslc_correct_cic{i}));
    err_conf_binned_hslc_correct(i) = sqrt(var(results_hslc.conf(indx_hslc_correct_cic{i}))/length(indx_hslc_correct_cic{i}));
    conf_binned_hslc_incorrect(i) = mean(results_hslc.conf(indx_hslc_incorrect_cic{i}));
    err_conf_binned_hslc_incorrect(i) = sqrt(var(results_hslc.conf(indx_hslc_incorrect_cic{i}))/length(indx_hslc_incorrect_cic{i}));

end
figure();
subplot(1,2,1)
% LH4(1) = plot(bin_mid_lshc_cic,conf_binned_lshc_correct,'-or','linewidth',2);
LH4(1) = errorbar(bin_mid_lshc_cic,conf_binned_lshc_correct,err_conf_binned_lshc_correct,'-or','linewidth',2);
L4{1} = 'Correct';
hold on;
% LH4(2) = plot(bin_mid_lshc_cic,conf_binned_lshc_incorrect,'--or','linewidth',2);
LH4(2) = errorbar(bin_mid_lshc_cic,conf_binned_lshc_incorrect,err_conf_binned_lshc_incorrect,'--or','linewidth',2);
L4{2} = 'Incorrect';
legend(LH4,L4,'box','off','fontsize',15);
ylabel('Mean conf','fontsize',20);
xlabel('Signal','fontsize',20);
title('LSHC','fontsize',20)
ylim([0.0 3]);
ax = gca;
set(ax, 'box','off');
ax.LineWidth=2;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,2,2)
% LH5(1) = plot(bin_mid_hslc_cic,conf_binned_hslc_correct,'-ob','linewidth',2);
LH5(1) = errorbar(bin_mid_hslc_cic,conf_binned_hslc_correct,err_conf_binned_hslc_correct,'-ob','linewidth',2);
L5{1} = 'Correct';
hold on;
% LH5(2) = plot(bin_mid_hslc_cic,conf_binned_hslc_incorrect,'--ob','linewidth',2);
LH5(2) = errorbar(bin_mid_hslc_cic,conf_binned_hslc_incorrect,err_conf_binned_hslc_incorrect,'--ob','linewidth',2);
L5{2} = 'Incorrect';
legend(LH5,L5,'box','off','fontsize',15);
ylabel('Mean conf','fontsize',20);
xlabel('Signal','fontsize',20);
title('HSLC','fontsize',20);
ylim([0.0 3]);
ax = gca;
set(ax, 'box','off');
ax.LineWidth=2;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;