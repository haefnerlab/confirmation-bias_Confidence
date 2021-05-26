clear all; close all; clc;
boots = 100;% number of bootstraps to get PK

hpr_ridge = logspace(-3, 3, 7);
hpr_ar1 = 0.0;
hpr_curvature = logspace(-3, 3, 7);
expt_type_noise = 2;
expt_type_ratio = 1;
standardize = 0; % z-score data or not
folds = 10; % how many folds of cross validation
num_frames = 10; % number of stimulus frames
cases = 2; % ratio case and noise case, so 2 cases. phase 2 for noise and 1 ratio
dir = 'RawDataNewConfidence';
subjects = {...
    'bpgConfidencetest-subject01';
    'bpgConfidencetest-subject02';
    'bpgConfidencetest-subject03';
    'bpgConfidencetest-subject05';
    'bpgConfidencetest-subject06';
    'bpgConfidencetest-subject07';
    'bpgConfidencetest-subject08';
    'bpgConfidencetest-subject09';
    'bpgConfidencetest-subject10';
    'bpgConfidencetest-subject11'
    };

subjects_id = {'subj 1';'subj 2';'subj 3';'subj 4'; 'subj 5'; 'subj 6';'subj 7';'subj 8';'subj 9';'subj 10'};
[num_sub,~] = size(subjects);

disp('Starting to find best hyperparameters for NOISE condition data across subjects....');
% [best_hprs_noise] = CustomRegression.combined_hprs_search(subjects,expt_type_noise,hpr_ridge,hpr_ar1,hpr_curvature,standardize,folds,dir);
best_hprs_noise = [1 0 443.6687];
disp('Starting to find best hyperparameters for RATIO condition data across subjects....');
% [best_hprs_ratio] = CustomRegression.combined_hprs_search(subjects,expt_type_ratio,hpr_ridge,hpr_ar1,hpr_curvature,standardize,folds,dir);
best_hprs_ratio = [17.1907 0 5.08022];
%%
for i=1:(num_sub)
    tic;
    figure();
    for j=1:cases
        phase = 3-j;
        if phase==expt_type_noise
            best_hprs = best_hprs_noise;
            disp(['Starting analysis for NOISE data of Subject ' num2str(i) ' ...']);
        elseif phase==expt_type_ratio
            best_hprs = best_hprs_ratio;
            disp(['Starting analysis for RATIO data of Subject ' num2str(i) ' ...']);
        end
        [params_boot,sobl,abbl_exp,best_hprs,trials,bin_centers,~,~,conf,conf_th,conf_tr,...
            subj_resp,subj_resp_err,confidence_resp,confidence_resp_err,ntrial_subj,conf_hist,conf_rt, ch_rt, direct_acc_conf, data, acc_m_l, acc_h_m,...
            log_bernoulli{i,j}]...
            = run_analysis_both_ratio_and_noise(subjects{i},phase,boots,best_hprs,standardize,dir);
        accuracy_m_l(i,j, :) = acc_m_l;
        accuracy_h_m(i,j, :) = acc_h_m;
        confidence_hist(i,j,:) = conf_hist;
        confidence_rt(i,j,:,:) = conf_rt;
        acc_conf(i,j,:,:) = direct_acc_conf;
        choice_rt(i,j,:,:) = ch_rt;
        confidence(i,j,:) = conf;
        confidence_threshold(i,j) = conf_th;
        confidence_trials_under_threshold(i,j) = conf_tr;
        subj_accuracy_mean(i,j,:) = subj_resp;
        subj_accuracy_err(i,j,:) = subj_resp_err;
        confidence_wrt_accuracy(i,j,:) = confidence_resp;
        confidence_wrt_accuracy_err(i,j,:) = confidence_resp_err;
        trials_accuracy_confidence_compare(i,j,:) = ntrial_subj;
        alpha(i,j,:) = [prctile(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)), 50) std(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)))/sqrt(size(params_boot,1))];
%         alpha(i,j,:) = [prctile(params_boot(:,end).^2, 50) std(params_boot(:,end).^2)/sqrt(size(params_boot,1))];
        bias(i,j) = prctile(params_boot(:, end-1), 50);
        temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 50);
        norm_temporal_kernel(i,j,:) = temporal_kernel(i,j,:)/mean(temporal_kernel(i,j,:));
        lo_temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
        hi_temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
        num_trials(i,j) = trials;
        all_exp(i,j,:,:) = abbl_exp;
        beta_all(i,j,:) = abbl_exp(:,2);
        beta(i,j) = prctile(squeeze(all_exp(i,j,:,2)),50);
        all_linear(i,j,:,:) = sobl;
        norm_all_linear(i,j,:,:) = [sobl(:,1)/mean(temporal_kernel(i,j,:)) sobl(:,2)/mean(temporal_kernel(i,j,:)) sobl(:,3) sobl(:,4)];
        slope(i,j) = prctile(squeeze(all_linear(i,j,:,1)),50);
        slope_all(i,j,:) = sobl(:,1);
        norm_slope_all(i,j,:) = norm_all_linear(i,j,:,1);
        norm_slope(i,j) = prctile(squeeze(norm_all_linear(i,j,:,1)),50);
        hprs_used(i,j,:) = best_hprs;
        data_sub{i,j}  = data;
        
        subplot(cases,3,3+3*(j-1))
        errorbar(1:num_frames,squeeze(temporal_kernel(i,j,1:num_frames)),squeeze(lo_temporal_kernel(i,j,:)),squeeze(hi_temporal_kernel(i,j,:)),'k','LineWidth',1)
        xlabel('Frames');
        ylabel('Weights');
        axis('tight');
        
        if phase==2
            subplot(cases,3,1+3*(j-1))
            plot((1:length(data.noise)), data.noise);
            hold on;
            plot(linspace(1,length(data.noise),100),confidence_threshold(i,j) * ones(1,100),'r','LineWidth',1);
            xlabel('Trials');
            ylabel('Noise Level');
            axis('tight');
            
            subplot(cases,3,2+3*(j-1))
            bins = 10;
            subject_pm_curve =[];
            uniq_vals = linspace(-0.8,0.8,bins);
            tr_kappa = data.sign_noise;
            noise_signal = data.ideal_frame_signals;
            for tt=1:(length(uniq_vals)-1)
                subj_resp(i,j,tt) = mean(data.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
                ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
            end
            noise_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
            errorbar(noise_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
            subject_pm_curve = alpha(i,j,1) * 0.5 + (1 - alpha(i,j,1)) * sigmoid(bias(i,j) + noise_signal*squeeze(temporal_kernel(i,j,:)));
            for tt=1:(length(uniq_vals)-1)
                subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))));
                ntrial_subj_pred(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
            end
            hold on;
            plot(noise_vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
            yline(0.5,'--k');
            xline(0.0,'--k');
            xlabel('Signed Kappa');
            ylabel('Percent chose left');
            xlim([-0.8 0.8])
            ylim([0.0 1.0])
        else
            subplot(cases,3,1+3*(j-1))
            plot((1:length(data.ratio)), data.ratio);
            hold on;
            plot(linspace(1,length(data.ratio),100),confidence_threshold(i,j) * ones(1,100),'r','LineWidth',1);
            xlabel('Trials');
            ylabel('Ratio Level');
            axis('tight');
            
            subplot(cases,3,2+3*(j-1))
            bins = 10;
            subject_pm_curve = [];
            uniq_vals = linspace(0,1,bins);
            tr_ratio = data.true_ratio;
            noise_signal = sign(data.ideal_frame_signals);
            for tt=1:(length(uniq_vals)-1)
                subj_resp(i,j,tt) = mean(data.choice(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
                ntrial_subj(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
            end
            ratio_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
            errorbar(ratio_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
            subject_pm_curve = alpha(i,j,1) * 0.5 + (1 - alpha(i,j,1)) * sigmoid(bias(i,j) + noise_signal*squeeze(temporal_kernel(i,j,:)));
            for tt=1:(length(uniq_vals)-1)
                subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)))));
                ntrial_subj_pred(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
            end
            hold on;
            plot(ratio_vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
            yline(0.5,'--k');
            xline(0.5,'--k');
            xlabel('Left Ratio');
            ylabel('Percent chose left');
            xlim([0.0 1.0])
            ylim([0.0 1.0])
        end
    end
    sgtitle(['Top Row: Noise (' num2str(confidence_trials_under_threshold(i,1)) ' trials)' ' and Bottom Row: Ratio (' num2str(confidence_trials_under_threshold(i,2)) ' trials)'  ' for Subject ' num2str(i)]);
    disp(['All analysis complete for Subject ' num2str(i) ' !!!!']);
    toc;
    disp('-----------------------------------------------------------------------------------------------------');
end
hold off;
%% Subjectwise plots
for i=1:(num_sub)
    figure();
    for j=1:cases
        phase = 3-j;
        data = data_sub{i,j};
        subplot(cases,3,3+3*(j-1))
        errorbar(1:num_frames,squeeze(temporal_kernel(i,j,1:num_frames)),squeeze(lo_temporal_kernel(i,j,:)),squeeze(hi_temporal_kernel(i,j,:)),'k','LineWidth',1)
        xlabel('Frames');
        ylabel('Weights');
        axis('tight');
        
        if phase==2
            subplot(cases,3,1+3*(j-1))
            plot((1:length(data.noise)), data.noise);
            hold on;
            plot(linspace(1,length(data.noise),100),confidence_threshold(i,j) * ones(1,100),'r','LineWidth',1);
            xlabel('Trials');
            ylabel('Noise Level');
            axis('tight');
            
            subplot(cases,3,2+3*(j-1))
            bins = 10;
            subject_pm_curve =[];
            uniq_vals = linspace(-0.8,0.8,bins);
            tr_kappa = data.sign_noise;
            noise_signal = data.ideal_frame_signals;
            for tt=1:(length(uniq_vals)-1)
                subj_resp(i,j,tt) = mean(data.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
                ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
            end
            noise_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
            errorbar(noise_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
            subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
            for tt=1:(length(uniq_vals)-1)
                subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))));
                ntrial_subj_pred(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
            end
            hold on;
            plot(noise_vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
            yline(0.5,'--k');
            xline(0.0,'--k');
            xlabel('Signed Kappa');
            ylabel('Percent chose left');
            xlim([-0.8 0.8])
            ylim([0.0 1.0])
        else
            subplot(cases,3,1+3*(j-1))
            plot((1:length(data.ratio)), data.ratio);
            hold on;
            plot(linspace(1,length(data.ratio),100),confidence_threshold(i,j) * ones(1,100),'r','LineWidth',1);
            xlabel('Trials');
            ylabel('Ratio Level');
            axis('tight');
            
            subplot(cases,3,2+3*(j-1))
            bins = 10;
            subject_pm_curve = [];
            uniq_vals = linspace(0,1,bins);
            tr_ratio = data.true_ratio;
            noise_signal = sign(data.ideal_frame_signals);
            for tt=1:(length(uniq_vals)-1)
                subj_resp(i,j,tt) = mean(data.choice(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
                ntrial_subj(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
            end
            ratio_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
            errorbar(ratio_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
            subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
            for tt=1:(length(uniq_vals)-1)
                subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)))));
                ntrial_subj_pred(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
            end
            hold on;
            plot(ratio_vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
            yline(0.5,'--k');
            xline(0.5,'--k');
            xlabel('Left Ratio');
            ylabel('Percent chose left');
            xlim([0.0 1.0])
            ylim([0.0 1.0])
        end
    end
    sgtitle(['Top Row: Noise (' num2str(confidence_trials_under_threshold(i,1)) ' trials)' ' and Bottom Row: Ratio (' num2str(confidence_trials_under_threshold(i,2)) ' trials)'  ' for Subject ' num2str(i)]);
end
%%
%%
bin_num = 15;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
for sub=1:num_sub
    subplot(2,5,sub)
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(log_bernoulli{sub,1});
    choice_used = data_sub{sub,1}.choice(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Log likelihood odds','Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Check how good weights predict noise stimulus behavior','fontsize',30);

bin_num = 15;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
for sub=1:num_sub
    subplot(2,5,sub)
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(log_bernoulli{sub,2});
    choice_used =  data_sub{sub,2}.choice(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Log likelihood odds','Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Check how good weights predict ratio stimulus behavior','fontsize',30);


%%

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    errorbar(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),squeeze(lo_temporal_kernel(i,1,:)),squeeze(hi_temporal_kernel(i,1,:)),'r','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    errorbar(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),squeeze(lo_temporal_kernel(i,2,:)),squeeze(hi_temporal_kernel(i,2,:)),'b','LineWidth',2);
    legend({['noise ' '(' num2str(num_trials(i,1)) ' trials)'],['ratio ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Temporal weights for all subjects','fontsize',30);

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'r','LineWidth',2);
    xlabel('Frames');
    ylabel('Norm Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'b','LineWidth',2);
    legend({['noise ' '(' num2str(num_trials(i,1)) ' trials)'],['ratio ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Normalized temporal weights for all subjects','fontsize',30);

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1)),'-or','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1)),'-ob','LineWidth',2);
    legend({['noise ' '(' num2str(num_trials(i,1)) ' trials)'],['ratio ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Exponential temporal weights for all subjects','fontsize',30);


figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50),'-or','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50),'-ob','LineWidth',2);
    legend({['noise ' '(' num2str(num_trials(i,1)) ' trials)'],['ratio ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Linear temporal weights for all subjects','fontsize',30);

%%
f = figure();
set(f,'defaultLegendAutoUpdate','off');
subplot(2,2,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),'r','LineWidth',0.2);
    cs1_sub(i,:) = squeeze(temporal_kernel(i,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'b','LineWidth',0.2);
    cs2_sub(i,:) = squeeze(temporal_kernel(i,2,1:num_frames));
    legend({['noise trials'],['ratio trials']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-or','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-ob','LineWidth',4);
title('Temporal Weights');

subplot(2,2,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'r','LineWidth',0.2);
    cs1_sub(i,:) = squeeze(norm_temporal_kernel(i,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'b','LineWidth',0.2);
    cs2_sub(i,:) = squeeze(norm_temporal_kernel(i,2,1:num_frames));
    legend({['noise trials'],['ratio trials']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-or','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-ob','LineWidth',4);
title('Normalized Temporal Weights');

subplot(2,2,4)
for i=1:(num_sub)
    plot(1:num_frames,prctile(squeeze(norm_all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(norm_all_linear(i,1,:,1)),50),'r','LineWidth',0.2);
    cs1_sub(i,:) = prctile(squeeze(norm_all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(norm_all_linear(i,1,:,1)),50);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(norm_all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(norm_all_linear(i,2,:,1)),50),'b','LineWidth',0.2);
    cs2_sub(i,:) = prctile(squeeze(norm_all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(norm_all_linear(i,2,:,1)),50);
    legend({['noise trials'],['ratio trials']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-or','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-ob','LineWidth',4);
title('Normalized Linear Weights');

subplot(2,2,3);
for i=1:(num_sub)
    plot(1:num_frames,prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1)),'r','LineWidth',0.2);
    cs1_sub(i,:) = prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1)),'b','LineWidth',0.2);
    cs2_sub(i,:) = prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1));
    legend({['noise trials'],['ratio trials']});
    title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-or','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-ob','LineWidth',4);
title('Exponential Weights');

sgtitle('All types of temporal weights fit across subjects','fontsize',30);

%% Precise Summary Figure
figure();
ax1 = subplot(2,6,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b');
    hold on;
end
hold on;
axis('tight')
xlabel('Frames');
ylabel('Norm Weights for Noise');
hold on;
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-ob','LineWidth',2);

ax2 = subplot(2,6,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'r');
    hold('on');
end
hold on;
xlabel('Frames');
ylabel('Norm Weights for Ratio');
hold on;
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-or','LineWidth',2);
hold on;
axis('tight')

linkaxes([ax1,ax2],'y')

subplot(2,6,3);
scatter(norm_slope(:,1),norm_slope(:,2),'k','filled');
err1 = std(squeeze(norm_slope_all(:,1,:))',1)./sqrt(boots);
err2 = std(squeeze(norm_slope_all(:,2,:))',1)./sqrt(boots);
v1 = var(squeeze(norm_slope(:,1))',1);
v2 = var(squeeze(norm_slope(:,2))',1);
hold on;
errorbar(norm_slope(:,1),norm_slope(:,2),err1, 'horizontal', 'LineStyle', 'none','color','b');
hold on;
errorbar(norm_slope(:,1),norm_slope(:,2),err2, 'vertical', 'LineStyle', 'none','color','r');
hold on;
scatter(norm_slope(:,1),norm_slope(:,2),'k','filled');
mn_s = -1.0;
mx_s = 0.0;
hold on;
plot(linspace(-2,2,10),linspace(-2,2,10),'k','LineWidth',2);
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),'g','filled');
xlabel('Norm Slopes for Noise ');
ylabel('Norm Slopes for Ratio');
hold('on');
axis('tight')

subplot(2,6,4);
scatter(beta(:,1),beta(:,2),'k','filled');
err1 = std(squeeze(beta_all(:,1,:))',1)./sqrt(boots);
err2 = std(squeeze(beta_all(:,2,:))',1)./sqrt(boots);
v1 = var(squeeze(beta(:,1))',1);
v2 = var(squeeze(beta(:,2))',1);
hold on;
errorbar(beta(:,1),beta(:,2),err1, 'horizontal', 'LineStyle', 'none','color','b');
hold on;
errorbar(beta(:,1),beta(:,2),err2, 'vertical', 'LineStyle', 'none','color','r');
hold on;
scatter(beta(:,1),beta(:,2),'k','filled');
mn_b = -1.0;
mx_b = 0;
hold on;
plot(linspace(-2,2,10),linspace(-2,2,10),'k','LineWidth',2);
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* beta(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* beta(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),'g','filled');
xlabel('Beta for Noise ');
ylabel('Beta for Ratio');
hold('on');
hold('on');
axis('tight')

subplot(2,6,5);
for i=1:(num_sub)
    plot(noise_vals,squeeze(subj_resp_pred(i,1,:)),'b','Linewidth',0.5);
    hold('on');
end
plot(noise_vals,squeeze(mean(subj_resp_pred(:,1,:),1)),'-ob','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1);
hold on;
xline(0.0,'k','Linewidth',1);
hold on;
ylim([0.0 1.0]);
xlim([-0.75 0.75]);
xlabel('Signed Kappa');
ylabel('Percent chose left for noise stimulus');

subplot(2,6,6);
for i=1:(num_sub)
    plot(ratio_vals,squeeze(subj_resp_pred(i,2,:)),'r','Linewidth',0.5);
    hold('on');
end
plot(ratio_vals,squeeze(mean(subj_resp_pred(:,2,:),1)),'-or','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1);
hold on;
xline(0.5,'k','Linewidth',1);
hold on;
ylim([0.0 1.0]);
xlim([0.05 0.95]);
xlabel('Ratio');
ylabel('Percent chose left for ratio stimulus');

subplot(2,6,7)
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'b');
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'r');
hold on;
scatter(mn1,mn2,'k','filled')
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
scatter(sum(mn_slp1,2),sum(mn_slp2,2),'g','filled');
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',2);
ylim([1.4 2.25]);
xlabel('Confidence in Noise Condition');
ylabel('Confidence in Ratio Condition');

subplot(2,6,8)
for i=1:(num_sub)
    vals(i,:) = [mean(squeeze(confidence(i,1,:))) mean(squeeze(confidence(i,2,:)))];
end
x = linspace(1,num_sub,num_sub);
bar(x,vals)
hold('on');
xlabel('Subject Number');
ylabel('Confidence');
legend({'Noise condition';'Ratio Condition'})

ax9 = subplot(2,6,9);
for i=1:(num_sub)
    [t,indt] = sort(squeeze(subj_accuracy_mean(i,1,:)));
    plot(t,squeeze(confidence_wrt_accuracy(i,1,indt)),'-ob','LineWidth',1.5);
    hold on;
    errorbar(squeeze(subj_accuracy_mean(i,1,:)),squeeze(confidence_wrt_accuracy(i,1,:)),squeeze(subj_accuracy_err(i,j,:)), 'horizontal', 'LineStyle', 'none','Color','k');
    hold on;
    errorbar(squeeze(subj_accuracy_mean(i,1,:)),squeeze(confidence_wrt_accuracy(i,1,:)),squeeze(confidence_wrt_accuracy_err(i,j,:)), 'vertical', 'LineStyle', 'none','Color','k');
    hold on;    
end
xlabel('Accuracy');
ylabel('Confidence');

ax10 = subplot(2,6,10);
for i=1:(num_sub)
    [t,indt] = sort(squeeze(subj_accuracy_mean(i,2,:)));
    plot(t,squeeze(confidence_wrt_accuracy(i,2,indt)),'-or','LineWidth',1.5);
    hold on;
    errorbar(squeeze(subj_accuracy_mean(i,2,:)),squeeze(confidence_wrt_accuracy(i,2,:)),squeeze(subj_accuracy_err(i,j,:)), 'horizontal', 'LineStyle', 'none','Color','k');
    hold on;
    errorbar(squeeze(subj_accuracy_mean(i,2,:)),squeeze(confidence_wrt_accuracy(i,2,:)),squeeze(confidence_wrt_accuracy_err(i,j,:)), 'vertical', 'LineStyle', 'none','Color','k');
    hold on;
end
xlabel('Accuracy');
ylabel('Confidence');
linkaxes([ax9,ax10],'y')

ax11 = subplot(2,6,11);
for i=1:(num_sub)
    errorbar([1 2 3],squeeze(confidence_rt(i,1,1,:)),squeeze(confidence_rt(i,1,2,:)), 'vertical', 'LineStyle', 'none','Color','k');
    hold on;
    plot([1 2 3],squeeze(confidence_rt(i,1,1,:)),'-ob','LineWidth',1.5);
    hold on;
end
xlabel('Confidence');
ylabel('Reaction time to report Confidence');

ax12 = subplot(2,6,12);
for i=1:(num_sub)
    errorbar([1 2 3],squeeze(confidence_rt(i,2,1,:)),squeeze(confidence_rt(i,2,2,:)), 'vertical', 'LineStyle', 'none','Color','k');
    hold on;
    plot([1 2 3],squeeze(confidence_rt(i,2,1,:)),'-or','LineWidth',1.5);
    hold on;
end
xlabel('Confidence');
ylabel('Reaction time to report Confidence');
linkaxes([ax11,ax12],'y')


%% Subjectwise Figure
figure();
ax1 = subplot(2,4,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'-o','LineWidth',1);
    hold on;
end
hold on;
xlim([1 num_frames]);
hold on;
yline(0.0,'k','Linewidth',1);
xlabel('Frames');
ylabel('Norm Weights for Noise');
hold on;
legend(subjects_id);

ax2 = subplot(2,4,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'-o','LineWidth',1);
    hold('on');
end
hold on;
yline(0.0,'k','Linewidth',1);
xlabel('Frames');
ylabel('Norm Weights for Ratio');
hold on;
legend(subjects_id);
xlim([1 num_frames]);

linkaxes([ax1,ax2],'y')

subplot(2,4,3);
for i=1:(num_sub)
    plot(noise_vals,squeeze(subj_resp_pred(i,1,:)),'-o','Linewidth',1);
    hold('on');
end
yline(0.5,'k','Linewidth',1);
hold on;
xline(0.0,'k','Linewidth',1);
hold on;
ylim([0.0 1.0]);
xlim([-0.75 0.75]);
legend(subjects_id);
xlabel('Signed Kappa');
ylabel('Percent chose left for noise stimulus');

subplot(2,4,4);
for i=1:(num_sub)
    plot(ratio_vals,squeeze(subj_resp_pred(i,2,:)),'-o','Linewidth',1);
    hold('on');
end
yline(0.5,'k','Linewidth',1);
hold on;
xline(0.5,'k','Linewidth',1);
hold on;
ylim([0.0 1.0]);
xlim([0.05 0.95]);
legend(subjects_id);
xlabel('Ratio');
ylabel('Percent chose left for ratio stimulus');

subplot(2,4,5)
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end
hold on;
for i=1:num_sub
    scatter(mn1(i),mn2(i),'filled')
    hold on;
end
legend(subjects_id,'AutoUpdate','off');
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',2);
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'k');
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'k');
xlabel('Confidence in Noise Condition');
ylabel('Confidence in Ratio Condition');

ax9 = subplot(2,4,6);
for i=1:(num_sub)
    [t,indt] = sort(squeeze(subj_accuracy_mean(i,1,:)));
    plot(t,squeeze(confidence_wrt_accuracy(i,1,indt)),'-o','LineWidth',1.5);
    hold on;
end
legend(subjects_id,'AutoUpdate','off');
for i=1:num_sub
    errorbar(squeeze(subj_accuracy_mean(i,1,:)),squeeze(confidence_wrt_accuracy(i,1,:)),squeeze(subj_accuracy_err(i,j,:)), 'horizontal', 'LineStyle', 'none','Color','k');
    hold on;
    errorbar(squeeze(subj_accuracy_mean(i,1,:)),squeeze(confidence_wrt_accuracy(i,1,:)),squeeze(confidence_wrt_accuracy_err(i,j,:)), 'vertical', 'LineStyle', 'none','Color','k');
    hold on;
end
xlim([0.5 1.0]);
xlabel('Accuracy');
ylabel('Confidence');

ax10 = subplot(2,4,7);
for i=1:(num_sub)
    [t,indt] = sort(squeeze(subj_accuracy_mean(i,2,:)));
    plot(t,squeeze(confidence_wrt_accuracy(i,2,indt)),'-o','LineWidth',1.5);
    hold on;
end
legend(subjects_id,'AutoUpdate','off');
for i=1:num_sub
    errorbar(squeeze(subj_accuracy_mean(i,2,:)),squeeze(confidence_wrt_accuracy(i,2,:)),squeeze(subj_accuracy_err(i,j,:)), 'horizontal', 'LineStyle', 'none','Color','k');
    hold on;
    errorbar(squeeze(subj_accuracy_mean(i,2,:)),squeeze(confidence_wrt_accuracy(i,2,:)),squeeze(confidence_wrt_accuracy_err(i,j,:)), 'vertical', 'LineStyle', 'none','Color','k');
    hold on;
end
xlim([0.5 1.0]);
xlabel('Accuracy');
ylabel('Confidence');
linkaxes([ax9,ax10],'y')

subplot(2,4,8)
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(accuracy_m_l(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(accuracy_h_m(sb,1,:)),0.95);%0.67);
    v1(sb) = var(squeeze(accuracy_m_l(sb,1,:)));
    v2(sb) = var(squeeze(accuracy_h_m(sb,1,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
    
    [mn11(sb),l11(sb),h11(sb),~] = meanci(squeeze(accuracy_m_l(sb,2,:)),0.95);%0.67);
    [mn22(sb),l22(sb),h22(sb),~] = meanci(squeeze(accuracy_h_m(sb,2,:)),0.95);%0.67);
    v11(sb) = var(squeeze(accuracy_m_l(sb,2,:)));
    v22(sb) = var(squeeze(accuracy_h_m(sb,2,:)));
    l11(sb) = mn11(sb) - l11(sb);
    h11(sb) = h11(sb) - mn11(sb);
    l22(sb) = mn22(sb) - l22(sb);
    h22(sb) = h22(sb) - mn22(sb);
end
hold on;
for i=1:num_sub
    scatter([mn1(i) mn11(i)],[mn2(i) mn22(i)],200,'filled');
    hold on;
end
legend(subjects_id,'AutoUpdate','off');
for i=1:num_sub
    errorbar(mn1(i),mn2(i),l1(i),h1(i), 'horizontal', 'LineStyle', 'none','color', 'k');
    hold on;
    errorbar(mn1(i),mn2(i),l11(i),h11(i), 'vertical', 'LineStyle', 'none','color', 'k');
    hold on;
    errorbar(mn11(i),mn22(i),l22(i),h22(i), 'vertical', 'LineStyle', 'none','color', 'k');
    hold on;
    errorbar(mn11(i),mn22(i),l2(i),h2(i), 'horizontal', 'LineStyle', 'none','color', 'k');
    hold on;
end
for i=1:num_sub
    d = double([mn11(i)-mn1(i) mn22(i)-mn2(i)]);
    q = quiver(mn1(i),mn2(i),d(1),d(2),'AutoScale','off');
    q.AutoScaleFactor=0.9;
    q.ShowArrowHead='on';
    q.AlignVertexCenters='on';
    q.LineWidth=1.5;
    q.MaxHeadSize=0.4;
    q.Head.LineStyle = 'solid';
    q.Color='black';
end
plot(linspace(-0.3,0.4,10),linspace(-0.3,0.4,10),'--k','LineWidth',1);
hold on;
xlabel('Accuracy Medium -  Accuracy Low');
ylabel('Accuracy High -  Accuracy Medium');

figure();
for i=1:num_sub
    subplot(2,5,i)
    plot([1 2 3],squeeze(confidence_hist(i,1,:)),'-ob','LineWidth',2)
    hold on;
    plot([1 2 3],squeeze(confidence_hist(i,2,:)),'-or','LineWidth',2)
    ylabel('Prob reported Confidence');
    xlabel('Confidence');
    title( subjects_id{i})
    legend({'Noise condition';'Ratio Condition'})
end

%%
figure()
subplot(2,4,3);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),'r');
    hold on;
end
hold on;
axis('tight')
xlabel({'Number of';' stimulus frames'},'FontSize',20);
ylabel('Temporal weights','FontSize',20);
hold on;
LH(1) = plot(1:num_frames,mean(squeeze((temporal_kernel(:,1,1:num_frames))),1),'-or','LineWidth',2);
L{1} = 'LSHC';
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'b');
    hold('on');
end
hold on;
xlabel({'Number of';' stimulus frames'},'FontSize',20);
ylabel('Temporal weights','FontSize',20);
hold on;
LH(2) = plot(1:num_frames,mean(squeeze((temporal_kernel(:,2,1:num_frames))),1),'-ob','LineWidth',2);
L{2} = 'HSLC';
hold on;
axis('tight')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylim([-0.1 0.25]);
legend(LH,L, 'Fontsize',20, 'Box','off');

subplot(2,4,1);
uniq_vals = linspace(-0.8,0.8,bins);
noise_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
for i=1:(num_sub)
    plot(noise_vals,squeeze(subj_resp_pred(i,1,:)),'r','Linewidth',0.5);
    hold('on');
end
plot(noise_vals,squeeze(mean(subj_resp_pred(:,1,:),1)),'-or','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1.5);
hold on;
xline(0.0,'k','Linewidth',1.5);
hold on;
ylim([0.0 1.0]);
xlim([-0.725 0.725]);
xlabel({'Signed Kappa';'(energy in favor of +45)'},'FontSize',20);
ylabel({'Prob chose +45 degrees';' for HSLC stimuli'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(2,4,2);
for i=1:(num_sub)
    plot(ratio_vals,squeeze(subj_resp_pred(i,2,:)),'b','Linewidth',0.5);
    hold('on');
end
plot(ratio_vals,squeeze(mean(subj_resp_pred(:,2,:),1)),'-ob','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1.5);
hold on;
xline(0.5,'k','Linewidth',1.5);
hold on;
ylim([0.0 1.0]);
xlim([0.05 0.95]);
xlabel({'Fraction of stimuli';' oriented +45 degrees'},'FontSize',20);
ylabel({'Prob chose +45 degrees';' for HSLC stimuli'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(2,4,4);
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'r','linewidth',1.5);
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'b','linewidth',1.5);
hold on;
scatter(mn1,mn2,'k','filled')
hold on;
ss = size(norm_slope(:,1),1);
mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
scatter(sum(mn_slp1,2),sum(mn_slp2,2),'g','filled');
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',1.5);
ylim([1.0 2.5]);
xlabel({'Mean confidence';' in LSHC trials'},'FontSize',20);
ylabel({'Mean confidence';' in HSLC trials'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;


%%
figure();
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end

subplot(1,2,1)
x_vals = mn1-mn2;
y_vals1 = beta(:,1)-beta(:,2);
scatter(x_vals,y_vals1,200,'k','filled')
temp = corrcoef(x_vals,y_vals1);
corr_bias_conf1 = temp(1,2);
xlabel({'Confidence difference (LSHC - HSLC)'},'FontSize',20);
ylabel({'Beta of weights difference (LSHC - HSLC)'},'FontSize',20);
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf1)],'FontSize',20,'FontWeight','bold');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,2,2)
x_vals = mn1-mn2;
y_vals2 = norm_slope(:,1) - norm_slope(:,2);
scatter(x_vals,y_vals2,200,'k','filled')
temp = corrcoef(x_vals,y_vals2);
corr_bias_conf2 = temp(1,2);
xlabel({'Confidence difference (LSHC - HSLC)'},'FontSize',20);
ylabel({'Norm slope of weights difference (LSHC - HSLC)'},'FontSize',20);
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf2)],'FontSize',20,'FontWeight','bold');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;


%%
figure();
ax9 = subplot(2,3,1);
for i=1:(num_sub)
    plot([1 2 3],squeeze(acc_conf(i,1,1,:)),'-or','LineWidth',1.5);
    hold on;
    errorbar([1 2 3],squeeze(acc_conf(i,1,1,:)),squeeze(acc_conf(i,1,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for LSHC', 'Fontsize',20);
ylabel({'Accuracy in LSHC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;

ax10 = subplot(2,3,4);
for i=1:(num_sub)
    plot([1 2 3],squeeze(acc_conf(i,2,1,:)),'-ob','LineWidth',1.5);
    hold on;
    errorbar([1 2 3],squeeze(acc_conf(i,2,1,:)),squeeze(acc_conf(i,2,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for HSLC', 'Fontsize',20);
ylabel({'Accuracy in HSLC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
linkaxes([ax9,ax10],'y')

ax11 = subplot(2,3,2);
for i=1:(num_sub)
    plot([1 2 3],squeeze(confidence_rt(i,1,1,:)),'-or','LineWidth',1.5);
    hold on;
    errorbar([1 2 3],squeeze(confidence_rt(i,1,1,:)),squeeze(confidence_rt(i,1,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for LSHC', 'Fontsize',20);
ylabel({'RT in ms to report';' Confidence in LSHC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;

ax12 = subplot(2,3,5);
for i=1:(num_sub)
    plot([1 2 3],squeeze(confidence_rt(i,2,1,:)),'-ob','LineWidth',1.5);
    hold on;
    errorbar([1 2 3],squeeze(confidence_rt(i,2,1,:)),squeeze(confidence_rt(i,2,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for HSLC', 'Fontsize',20);
ylabel({'RT in ms to report';' Confidence in HSLC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
linkaxes([ax11,ax12],'y')

ax13 = subplot(2,3,3);
for i=1:(num_sub)
    plot([1 2 3],squeeze(choice_rt(i,1,1,:)),'-or','LineWidth',1.5);
    hold on;
    errorbar([1 2 3],squeeze(choice_rt(i,1,1,:)),squeeze(choice_rt(i,1,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for LSHC', 'Fontsize',20);
ylabel({'RT in ms to report';' choice in LSHC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;

ax14 = subplot(2,3,6);
for i=1:(num_sub)
    plot([1 2 3],squeeze(choice_rt(i,2,1,:)),'-ob','LineWidth',1.5);
    hold on;
    errorbar([1 2 3],squeeze(choice_rt(i,2,1,:)),squeeze(choice_rt(i,2,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for HSLC', 'Fontsize',20');
ylabel({'RT in ms to report';' choice in HSLC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
linkaxes([ax13,ax14],'y')


%%
figure()
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'r','linewidth',1.5);
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'b','linewidth',1.5);
hold on;
scatter(mn1,mn2,'k','filled')
hold on;
ss = size(norm_slope(:,1),1);
mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
scatter(sum(mn_slp1,2),sum(mn_slp2,2),'g','filled');
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',1.5);
ylim([1.0 2.5]);
xlabel({'Mean confidence';' in LSHC trials'},'FontSize',20);
ylabel({'Mean confidence';' in HSLC trials'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,4,3)
for i=1:(num_sub)
    plot([1 2 3],squeeze(acc_conf(i,1,1,:)),'-or','LineWidth',2);
    hold on;
    errorbar([1 2 3],squeeze(acc_conf(i,1,1,:)),squeeze(acc_conf(i,1,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for LSHC', 'Fontsize',20);
ylabel({'Accuracy in LSHC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;

ax10 = subplot(1,4,4);
for i=1:(num_sub)
    plot([1 2 3],squeeze(acc_conf(i,2,1,:)),'-ob','LineWidth',2);
    hold on;
    errorbar([1 2 3],squeeze(acc_conf(i,2,1,:)),squeeze(acc_conf(i,2,2,:)), 'vertical', 'LineStyle', 'none','Color','k','linewidth',2);
    hold on;
end
xlabel('Confidence for HSLC', 'Fontsize',20);
ylabel({'Accuracy in HSLC'}, 'Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;

subplot(1,4,1);
uniq_vals = linspace(-0.8,0.8,bins);
noise_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
for i=1:(num_sub)
    plot(noise_vals,squeeze(subj_resp_pred(i,1,:)),'r','Linewidth',0.5);
    hold('on');
end
plot(noise_vals,squeeze(mean(subj_resp_pred(:,1,:),1)),'-or','Linewidth',4);
hold('on');
yline(0.5,'k','Linewidth',1.5);
hold on;
xline(0.0,'k','Linewidth',1.5);
hold on;
ylim([0.0 1.0]);
xlim([-0.725 0.725]);
xlabel({'Signed Kappa';'(energy in favor of +45)'},'FontSize',20);
ylabel({'Prob chose +45 degrees';' for LSHC stimuli'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,4,2);
for i=1:(num_sub)
    plot(ratio_vals,squeeze(subj_resp_pred(i,2,:)),'b','Linewidth',0.5);
    hold('on');
end
plot(ratio_vals,squeeze(mean(subj_resp_pred(:,2,:),1)),'-ob','Linewidth',4);
hold('on');
yline(0.5,'k','Linewidth',1.5);
hold on;
xline(0.5,'k','Linewidth',1.5);
hold on;
ylim([0.0 1.0]);
xlim([0.05 0.95]);
xlabel({'Fraction of stimuli';' oriented +45 degrees'},'FontSize',20);
ylabel({'Prob chose +45 degrees';' for HSLC stimuli'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;


%%
figure()
subplot(1,3,1)
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'r');
    hold on;
end
hold on;
axis('tight')
xlabel({'Number of';' stimulus frames'},'FontSize',20);
ylabel('Temporal weights','FontSize',20);
hold on;
LH(1) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-or','LineWidth',4);
L{1} = 'LSHC';
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'b');
    hold('on');
end
hold on;
xlabel({'Number of';' stimulus frames'},'FontSize',20);
ylabel({'Normalized';' temporal weights'},'FontSize',20);
hold on;
LH(2) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-ob','LineWidth',4);
L{2} = 'HSLC';
hold on;
axis('tight')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,3,2)
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'r','linewidth',3);
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'b','linewidth',3);
hold on;
scatter(mn1,mn2,100,'k','filled')
hold on;
ss = size(norm_slope(:,1),1);
mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
% scatter(sum(mn_slp1,1),sum(mn_slp2,1),'g','filled');
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',1.5);
ylim([1.0 2.5]);
xlabel({'Mean confidence';' in LSHC trials'},'FontSize',20);
ylabel({'Mean confidence';' in HSLC trials'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end

subplot(1,3,3)
x_vals = mn1-mn2;
y_vals1 = beta(:,1)-beta(:,2);
scatter(x_vals,y_vals1,200,'k','filled')
[temp,p] = corrcoef(x_vals,y_vals1);
corr_bias_conf1 = temp(1,2);
xlabel({'Confidence difference';' (LSHC - HSLC)'},'FontSize',20);
ylabel({'Difference of exponential';' slope parameter of weights';' (LSHC - HSLC)'},'FontSize',20);
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf1)],'FontSize',20,'FontWeight','bold');
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf1) ',p=' num2str(p(1,2))],'FontSize',20,'FontWeight','bold');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

%%
figure()
subplot(1,3,1)
x_vals = mn1;
y_vals1 = beta(:,1);
scatter(x_vals,y_vals1,200,'k','filled')
[temp,p] = corrcoef(x_vals,y_vals1);
corr_bias_conf1 = temp(1,2);
xlabel({'Confidence in LSHC'},'FontSize',20);
ylabel({'Beta of weights in LSHC'},'FontSize',20);
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf1) ',p=' num2str(p(1,2))],'FontSize',20,'FontWeight','bold');
title(['corr = ' num2str(corr_bias_conf1) ', p=' num2str(p(1,2))],'Fontsize',15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,3,2)
x_vals = mn2;
y_vals1 = beta(:,2);
scatter(x_vals,y_vals1,200,'k','filled')
[temp,p] = corrcoef(x_vals,y_vals1);
corr_bias_conf1 = temp(1,2);
xlabel({'Confidence in HSLC'},'FontSize',20);
ylabel({'Beta of weights in HSLC'},'FontSize',20);
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf1) ',p=' num2str(p(1,2))],'FontSize',20,'FontWeight','bold');
title(['corr = ' num2str(corr_bias_conf1) ', p=' num2str(p(1,2))],'Fontsize',15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(1,3,3)
x_vals = mn1-mn2;
y_vals1 = beta(:,1)-beta(:,2);
scatter(x_vals,y_vals1,200,'k','filled')
[temp,p] = corrcoef(x_vals,y_vals1);
corr_bias_conf1 = temp(1,2);
xlabel({'Confidence difference (LSHC - HSLC)'},'FontSize',20);
ylabel({'Beta of weights difference (LSHC - HSLC)'},'FontSize',20);
text(0.0,-0.7,['corr = ' num2str(corr_bias_conf1) ',p=' num2str(p(1,2))],'FontSize',20,'FontWeight','bold');
title(['corr = ' num2str(corr_bias_conf1) ', p=' num2str(p(1,2))],'Fontsize',15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;


%%
chosen_sub = [5 9];
figure();
ax1=subplot(1,4,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'r');
    hold on;
end
hold on;
axis('tight')
hold on;
LH(1) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-or','LineWidth',4);
L{1} = 'LSHC';
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'b');
    hold('on');
end
hold on;
xlabel({'Number of stimulus frames'},'FontSize',20);
ylabel({'Temporal weights'},'FontSize',20);
hold on;
LH(2) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-ob','LineWidth',4);
L{2} = 'HSLC';
plot(1:10, zeros(1,10),'-k','linewidth',2);
hold on;
axis('tight')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylim([-1 3.5]);
xlim([1 10])
xticks([1 5 10])
legend(LH,L, 'Fontsize',20, 'Box','off');

subplot(1,4,1)
for i=1:(num_sub)
    if i==5||i==9
        plot([1 2],beta(i,:),'--ko','LineWidth',2.0);
    else
        plot([1 2],beta(i,:),'-ko','LineWidth',2.0);
    end
    hold on;
end
hold on;
scatter(ones(1,num_sub),beta(:,1),200,'r','filled');
hold on;
scatter(ones(1,num_sub)*2,beta(:,2),200,'b','filled');
hold on;
ylabel({'Exponential slope parameter of weights'},'Fontsize',20)
xlim([0.9 2.1]);
ylim([-0.75 0.25]);
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel',{'LSHC' 'HSLC'});
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
legend(LH,L, 'Fontsize',20, 'Box','off');

subplot(1,2,2)
for sb=1:num_sub
    [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
    [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
    l1(sb) = mn1(sb) - l1(sb);
    h1(sb) = h1(sb) - mn1(sb);
    l2(sb) = mn2(sb) - l2(sb);
    h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'r','linewidth',3);
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'b','linewidth',3);
hold on;
scatter(mn1,mn2,100,'k','filled')
hold on;
ss = size(norm_slope(:,1),1);
mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
% scatter(sum(mn_slp1,1),sum(mn_slp2,1),'g','filled');
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',2);
ylim([1.0 2.5]);
xlabel({'Mean confidence in LSHC trials'},'FontSize',20);
ylabel({'Mean confidence in HSLC trials'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

%%
chosen_sub = [5 9];
figure()
ax1=subplot(1,4,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'r');
    hold on;
end
hold on;
axis('tight')
% xlabel({'Number of';' stimulus frames'},'FontSize',20);
hold on;
LH(1) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-or','LineWidth',4);
L{1} = 'LSHC';
% title('Noise Norm Weights')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
% ax2=subplot(2,5,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'b');
    hold('on');
end
hold on;
% xlabel({'Number of';' stimulus frames'},'FontSize',20);
xlabel({'Number of stimulus frames'},'FontSize',20);
ylabel({'Temporal weights'},'FontSize',20);
hold on;
LH(2) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-ob','LineWidth',4);
L{2} = 'HSLC';
% title('Ratio Norm Weights')
plot(1:10, zeros(1,10),'-k','linewidth',2);
hold on;
axis('tight')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylim([-1 3.5]);
xlim([1 10])
xticks([1 5 10])
legend(LH,L, 'Fontsize',20, 'Box','off');

subplot(1,4,1)
for i=1:(num_sub)
    if i==5||i==9
        plot([1 2],beta(i,:),'--ko','LineWidth',2.0);
        
    else
        plot([1 2],beta(i,:),'-ko','LineWidth',2.0);
    end
    hold on;
end
hold on;
scatter(ones(1,num_sub),beta(:,1),200,'r','filled');
hold on;
scatter(ones(1,num_sub)*2,beta(:,2),200,'b','filled');
hold on;
% xlabel('Accumulated evidence','Fontsize',20)
ylabel({'Exponential slope parameter of weights'},'Fontsize',20)
xlim([0.9 2.1]);
ylim([-0.75 0.25]);
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel',{'LSHC' 'HSLC'});
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
legend(LH,L, 'Fontsize',20, 'Box','off');


subplot(1,2,2)
for sb=1:num_sub
%     l1(sb) = prctile(squeeze(confidence(sb,1,:)), 50) - prctile(squeeze(confidence(sb,1,:)), 16);
%     h1(sb) = prctile(squeeze(confidence(sb,1,:)), 84) - prctile(squeeze(confidence(sb,1,:)), 50);
%     l2(sb) = prctile(squeeze(confidence(sb,2,:)), 50) - prctile(squeeze(confidence(sb,2,:)), 16);
%     h2(sb) = prctile(squeeze(confidence(sb,2,:)), 84) - prctile(squeeze(confidence(sb,2,:)), 50);
    l1(sb) = prctile(squeeze(confidence(sb,1,:))', 50) - prctile(squeeze(confidence(sb,1,:))', 2.5);
    h1(sb) = prctile(squeeze(confidence(sb,1,:))', 97.5) - prctile(squeeze(confidence(sb,1,:))', 50);
    l2(sb) = prctile(squeeze(confidence(sb,2,:))', 50) - prctile(squeeze(confidence(sb,2,:))', 2.5);
    h2(sb) = prctile(squeeze(confidence(sb,2,:))', 97.5) - prctile(squeeze(confidence(sb,2,:))', 50);
%     [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
%     [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
%     l1(sb) = mn1(sb) - l1(sb);
%     h1(sb) = h1(sb) - mn1(sb);
%     l2(sb) = mn2(sb) - l2(sb);
%     h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'r','linewidth',3);
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'b','linewidth',3);
hold on;
scatter(mn1,mn2,100,'k','filled')
hold on;
ss = size(norm_slope(:,1),1);
% mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
% mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
% scatter(sum(mn_slp1,2),sum(mn_slp2,2),'g','filled');
scatter(mean(mn1),mean(mn2),100,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2);
% hold on;
% errorbar(mean(mn1),mean(mn2),std(mn1')/sqrt(num_sub), 'horizontal', 'LineStyle', 'none','color', 'g','linewidth',3);
% hold on;
% errorbar(mean(mn1),mean(mn2),std(mn2')/sqrt(num_sub), 'vertical', 'LineStyle', 'none','color', 'g','linewidth',3);
% hold on;
% l_mn1 = prctile(squeeze(mn1), 50) - prctile(squeeze(mn1), 2.5);
% h_mn1 = prctile(squeeze(mn1), 97.5) - prctile(squeeze(mn1), 50);
% l_mn2 = prctile(squeeze(mn2), 50) - prctile(squeeze(mn2), 2.5);
% h_mn2 = prctile(squeeze(mn2), 97.5) - prctile(squeeze(mn2), 50);
% l_mn1 = prctile(squeeze(mean(confidence(:,1,:),1))', 50) - prctile(squeeze(mean(confidence(:,1,:),1))', 2.5);
% h_mn1 = prctile(squeeze(mean(confidence(:,1,:),1))', 97.5) - prctile(squeeze(mean(confidence(:,1,:),1))', 50);
% l_mn2 = prctile(squeeze(mean(confidence(:,2,:),1))', 50) - prctile(squeeze(mean(confidence(:,2,:),1))', 2.5);
% h_mn2 = prctile(squeeze(mean(confidence(:,2,:),1))', 97.5) - prctile(squeeze(mean(confidence(:,2,:),1))', 50);
hold on;
errorbar(mean(mn1),mean(mn2),2*std(mn1)/sqrt(num_sub), 'horizontal', 'LineStyle', 'none','color', 'g','linewidth',1.5);
hold on;
errorbar(mean(mn1),mean(mn2),2*std(mn2)/sqrt(num_sub), 'vertical', 'LineStyle', 'none','color', 'g','linewidth',1.5);
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',2);
ylim([1.0 2.5]);
xlim([1 3]);
xlabel({'Mean confidence in LSHC trials'},'FontSize',20);
ylabel({'Mean confidence in HSLC trials'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
%%
fx = @(mu) normcdf(abs(mean(mu)-0.0)./std(mu),'upper');
for sub=1:num_sub
    confidence_diff(sub,:) = confidence(sub,1,:) - confidence(sub,2,:);
    beta_diff(sub,:) = beta_all(sub,1,:) - beta_all(sub,2,:);
    norm_slope_diff(sub,:) = norm_slope_all(sub,1,:) - norm_slope_all(sub,2,:);
    p_values_conf(sub) = fx(squeeze(confidence_diff(sub,:)));
    p_values_beta(sub) = fx(squeeze(beta_diff(sub,:)));
    p_values_slope(sub) = fx(squeeze(norm_slope_diff(sub,:)));
    if (p_values_conf(sub))<0.001
        stars_conf{sub} = '***';
    elseif (p_values_conf(sub))<0.01
        stars_conf{sub} = '**';
    elseif (p_values_conf(sub))<0.05
        stars_conf{sub} = '*';
    else
        stars_conf{sub} = ' ';
    end
    if (p_values_slope(sub))<0.001
        stars_slope{sub} = '***';
    elseif (p_values_slope(sub))<0.01
        stars_slope{sub} = '**';
    elseif (p_values_slope(sub))<0.05
        stars_slope{sub} = '*';
    else
        stars_slope{sub} = ' ';
    end
    if (p_values_beta(sub))<0.001
        stars_beta{sub} = '***';
    elseif (p_values_beta(sub))<0.01
        stars_beta{sub} = '**';
    elseif (p_values_beta(sub))<0.05
        stars_beta{sub} = '*';
    else
        stars_beta{sub} = ' ';
    end
    
    conf_diff_mn(sub) = mean(squeeze(confidence_diff(sub,:)));
    beta_diff_mn(sub) = mean(squeeze(beta_diff(sub,:)));
    slope_diff_mn(sub) = mean(squeeze(norm_slope_diff(sub,:)));
    l_conf_diff(sub) = prctile(squeeze(confidence_diff(sub,:)), 50) - prctile(squeeze(confidence_diff(sub,:)), 2.5);
    h_conf_diff(sub) = prctile(squeeze(confidence_diff(sub,:)), 97.5) - prctile(squeeze(confidence_diff(sub,:)), 50);
    l_beta_diff(sub) = prctile(squeeze(beta_diff(sub,:)), 50) - prctile(squeeze(beta_diff(sub,:)), 2.5);
    h_beta_diff(sub) = prctile(squeeze(beta_diff(sub,:)), 97.5) - prctile(squeeze(beta_diff(sub,:)), 50);
    l_slope_diff(sub) = prctile(squeeze(norm_slope_diff(sub,:)), 50) - prctile(squeeze(norm_slope_diff(sub,:)), 2.5);
    h_slope_diff(sub) = prctile(squeeze(norm_slope_diff(sub,:)), 97.5) - prctile(squeeze(norm_slope_diff(sub,:)), 50);   
end


figure();
subplot(3,1,1)
hold on;
for sub=1:num_sub
    if sub==5 || sub==9
        bar(sub,mean(confidence_diff(sub,:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        hold on;
    else
        bar(sub,mean(confidence_diff(sub,:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
        hold on;
    end
end
hold on;
errorbar(1:num_sub,mean(confidence_diff,2),squeeze(l_conf_diff),squeeze(h_conf_diff),'ok','LineWidth', 1.5)
hold on;
text([1:num_sub]-0.1,squeeze(mean(confidence_diff,2))+0.2,stars_conf,'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Confidence in LSHC ';' - Confidence in HSLC'},'Fontsize',20);
yline(0.0,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
ylim([-0.25 1]);
xlim([0.5 10.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = [1 2 3 4 5 6 7 8 9 10];
hold on;

% subplot(3,1,2)
% hold on;
% for sub=1:num_sub
%     if sub==5 || sub==9
%         bar(sub,mean(beta_diff(sub,:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
%         hold on;
%     else
%         bar(sub,mean(beta_diff(sub,:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
%         hold on;
%     end
% end
% bar(1:num_sub,mean(beta_diff,2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
% hold on;
% errorbar(1:num_sub,mean(beta_diff,2),squeeze(l_beta_diff),squeeze(h_beta_diff),'ok','LineWidth', 1.5)
% hold on;
% text([1:num_sub]-0.1,ones(1,num_sub)*0.5,stars_beta,'FontSize',15,'FontWeight','bold');
% xlabel('Subject Number','Fontsize',20)
% ylabel({'Beta in LSHC ';' - Beta in HSLC'},'Fontsize',20);
% yline(0.0,'k','LineWidth',1.5);
% ax = gca;
% ax.LineWidth=2;
% ylim([-1 1]);
% xlim([0.5 10.5]);
% ax.XAxis.FontSize = 20;
% ax.YAxis.FontSize = 20;
% ax.XTick = [1 2 3 4 5 6 7 8 9 10];
% hold on;

% subplot(3,1,3)
% hold on;
% bar(1:num_sub,mean(norm_slope_diff,2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
% hold on;
% % bar([5 9],mean(confidence_diff([5 9],:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
% hold on;
% errorbar(1:num_sub,mean(norm_slope_diff,2),squeeze(l_slope_diff),squeeze(h_slope_diff),'ok','LineWidth', 1.5)
% hold on;
% text([1:num_sub]-0.1,ones(1,num_sub)*0.5,stars_slope,'FontSize',15,'FontWeight','bold');
% xlabel('Subject Number','Fontsize',20)
% ylabel({'Slope in LSHC ';' - Slope in HSLC'},'Fontsize',20);
% yline(0.0,'k','LineWidth',1.5);
% ax = gca;
% ax.LineWidth=2;
% ylim([-1 1]);
% xlim([0.5 10.5]);
% ax.XAxis.FontSize = 20;
% ax.YAxis.FontSize = 20;
% ax.XTick = [1 2 3 4 5 6 7 8 9 10];
% hold on;

%%
fx_mn = @(mu,nt) normcdf(abs(mean(mu)-0.0)./(std(mu)/sqrt(nt)),'upper');
confidence_diff_all = mean(confidence_diff,2);
confidence_diff_all = confidence_diff_all(:);
beta_diff_all = mean(beta_diff,2);
beta_diff_all = beta_diff_all(:);
norm_slope_diff_all = mean(norm_slope_diff,2);
norm_slope_diff_all = norm_slope_diff_all(:);
p_values_conf_all = fx_mn(squeeze(confidence_diff_all),num_sub);
p_values_beta_all = fx_mn(squeeze(beta_diff_all),num_sub);
p_values_slope_all = fx_mn(squeeze(norm_slope_diff_all),num_sub);
if (p_values_conf_all)<(0.001)
    stars_conf_all = '***';
elseif (p_values_conf_all)<(0.01)
    stars_conf_all = '**';
elseif (p_values_conf_all)<(0.05)
    stars_conf_all = '*';
else
    stars_conf_all = ' ';
end
if (p_values_slope_all)<(0.001)
    stars_slope_all = '***';
elseif (p_values_slope_all)<(0.01)
    stars_slope_all = '**';
elseif (p_values_slope_all)<(0.05)
    stars_slope_all = '*';
else
    stars_slope_all = ' ';
end
if (p_values_beta_all)<(0.001)
    stars_beta_all = '***';
elseif (p_values_beta_all)<(0.01)
    stars_beta_all = '**';
elseif (p_values_beta_all)<(0.05)
    stars_beta_all = '*';
else
    stars_beta_all = ' ';
end

conf_diff_mn_all = mean(squeeze(confidence_diff_all));
beta_diff_mn_all = mean(squeeze(beta_diff_all));
slope_diff_mn_all = mean(squeeze(norm_slope_diff_all));
% l_conf_diff_all = prctile(squeeze(confidence_diff_all), 50) - prctile(squeeze(confidence_diff_all), 2.5);
% h_conf_diff_all = prctile(squeeze(confidence_diff_all), 97.5) - prctile(squeeze(confidence_diff_all), 50);
% l_beta_diff_all = prctile(squeeze(beta_diff_all), 50) - prctile(squeeze(beta_diff_all), 2.5);
% h_beta_diff_all = prctile(squeeze(beta_diff_all), 97.5) - prctile(squeeze(beta_diff_all), 50);
% l_slope_diff_all = prctile(squeeze(norm_slope_diff_all), 50) - prctile(squeeze(norm_slope_diff_all), 2.5);
% h_slope_diff_all = prctile(squeeze(norm_slope_diff_all), 97.5) - prctile(squeeze(norm_slope_diff_all), 50);


figure();
hold on;
bar([1 2 3],[conf_diff_mn_all beta_diff_mn_all slope_diff_mn_all], 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
errorbar([1 2 3],[conf_diff_mn_all beta_diff_mn_all slope_diff_mn_all],2*[std(confidence_diff_all)/sqrt(num_sub) std(beta_diff_all)/sqrt(num_sub) std(norm_slope_diff_all)/sqrt(num_sub)],...
    2*[std(confidence_diff_all)/sqrt(num_sub) std(beta_diff_all)/sqrt(num_sub) std(norm_slope_diff_all)/sqrt(num_sub)],'ok','LineWidth', 1.5)
hold on;
text([1 2 3]-0.1,ones(1,3)*0.5,{stars_conf_all; stars_beta_all; stars_slope_all},'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Difference between LSHC';'and HSLC conditions'},'Fontsize',20);
yline(0.0,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
ylim([-2 1]);
xlim([0.5 3.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = [1 2 3];
ax.XTickLabel = {'Conf'; 'Beta'; 'Slope'};
hold on;


%%
chosen_sub = [5 9];
figure()
ax1=subplot(2,2,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'r');
    hold on;
end
hold on;
axis('tight')
% xlabel({'Number of';' stimulus frames'},'FontSize',20);
hold on;
LH(1) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-or','LineWidth',4);
L{1} = 'LSHC';
% title('Noise Norm Weights')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
% ax2=subplot(2,5,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'b');
    hold('on');
end
hold on;
% xlabel({'Number of';' stimulus frames'},'FontSize',20);
xlabel({'Number of stimulus frames'},'FontSize',20);
ylabel({'Temporal weights'},'FontSize',20);
hold on;
LH(2) = plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-ob','LineWidth',4);
L{2} = 'HSLC';
% title('Ratio Norm Weights')
plot(1:10, zeros(1,10),'-k','linewidth',2);
hold on;
axis('tight')
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylim([-1 3.5]);
xlim([1 10])
xticks([1 5 10])
legend(LH,L, 'Fontsize',20, 'Box','off');

subplot(2,2,2)
for i=1:(num_sub)
    if i==5||i==9
        plot([1 2],beta(i,:),'--ko','LineWidth',2.0);
        
    else
        plot([1 2],beta(i,:),'-ko','LineWidth',2.0);
    end
    hold on;
end
hold on;
scatter(ones(1,num_sub),beta(:,1),200,'r','filled');
hold on;
scatter(ones(1,num_sub)*2,beta(:,2),200,'b','filled');
hold on;
% xlabel('Accumulated evidence','Fontsize',20)
ylabel({'Exponential slope parameter'; 'of weights'},'Fontsize',20)
xlim([0.9 2.1]);
ylim([-0.75 0.25]);
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel',{'LSHC' 'HSLC'});
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
legend(LH,L, 'Fontsize',20, 'Box','off');


subplot(2,2,3)
for sb=1:num_sub
%     l1(sb) = prctile(squeeze(confidence(sb,1,:)), 50) - prctile(squeeze(confidence(sb,1,:)), 16);
%     h1(sb) = prctile(squeeze(confidence(sb,1,:)), 84) - prctile(squeeze(confidence(sb,1,:)), 50);
%     l2(sb) = prctile(squeeze(confidence(sb,2,:)), 50) - prctile(squeeze(confidence(sb,2,:)), 16);
%     h2(sb) = prctile(squeeze(confidence(sb,2,:)), 84) - prctile(squeeze(confidence(sb,2,:)), 50);
    l1(sb) = prctile(squeeze(confidence(sb,1,:))', 50) - prctile(squeeze(confidence(sb,1,:))', 2.5);
    h1(sb) = prctile(squeeze(confidence(sb,1,:))', 97.5) - prctile(squeeze(confidence(sb,1,:))', 50);
    l2(sb) = prctile(squeeze(confidence(sb,2,:))', 50) - prctile(squeeze(confidence(sb,2,:))', 2.5);
    h2(sb) = prctile(squeeze(confidence(sb,2,:))', 97.5) - prctile(squeeze(confidence(sb,2,:))', 50);
%     [mn1(sb),l1(sb),h1(sb),~] = meanci(squeeze(confidence(sb,1,:)),0.95);%0.67);
%     [mn2(sb),l2(sb),h2(sb),~] = meanci(squeeze(confidence(sb,2,:)),0.95);%0.67);
    v1(sb) = var(squeeze(confidence(sb,1,:)));
    v2(sb) = var(squeeze(confidence(sb,2,:)));
%     l1(sb) = mn1(sb) - l1(sb);
%     h1(sb) = h1(sb) - mn1(sb);
%     l2(sb) = mn2(sb) - l2(sb);
%     h2(sb) = h2(sb) - mn2(sb);
end
hold on;
errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none','color', 'r','linewidth',3);
hold on;
errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none','color', 'b','linewidth',3);
% set(eb(1), 'color', 'b')
% set(eb(2), 'color', 'r')
hold on;
scatter(mn1,mn2,100,'k','filled')
hold on;
ss = size(norm_slope(:,1),1);
mn_slp1 = ((1./v1)./sum(1./v1)) .* mn1;
mn_slp2 = ((1./v2)./sum(1./v2)) .* mn2;
scatter(sum(mn_slp1,2),sum(mn_slp2,2),100,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2);
% scatter(mean(mn1),mean(mn2),100,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2);
% hold on;
% errorbar(mean(mn1),mean(mn2),std(mn1')/sqrt(num_sub), 'horizontal', 'LineStyle', 'none','color', 'g','linewidth',3);
% hold on;
% errorbar(mean(mn1),mean(mn2),std(mn2')/sqrt(num_sub), 'vertical', 'LineStyle', 'none','color', 'g','linewidth',3);
% hold on;
% l_mn1 = prctile(squeeze(mn1), 50) - prctile(squeeze(mn1), 2.5);
% h_mn1 = prctile(squeeze(mn1), 97.5) - prctile(squeeze(mn1), 50);
% l_mn2 = prctile(squeeze(mn2), 50) - prctile(squeeze(mn2), 2.5);
% h_mn2 = prctile(squeeze(mn2), 97.5) - prctile(squeeze(mn2), 50);
% l_mn1 = prctile(squeeze(mean(confidence(:,1,:),1))', 50) - prctile(squeeze(mean(confidence(:,1,:),1))', 2.5);
% h_mn1 = prctile(squeeze(mean(confidence(:,1,:),1))', 97.5) - prctile(squeeze(mean(confidence(:,1,:),1))', 50);
% l_mn2 = prctile(squeeze(mean(confidence(:,2,:),1))', 50) - prctile(squeeze(mean(confidence(:,2,:),1))', 2.5);
% h_mn2 = prctile(squeeze(mean(confidence(:,2,:),1))', 97.5) - prctile(squeeze(mean(confidence(:,2,:),1))', 50);
hold on;
errorbar(sum(mn_slp1,2),sum(mn_slp2,2),2*std(mn1)/sqrt(num_sub), 'horizontal', 'LineStyle', 'none','color', 'g','linewidth',1.5);
hold on;
errorbar(sum(mn_slp1,2),sum(mn_slp2,2),2*std(mn2)/sqrt(num_sub), 'vertical', 'LineStyle', 'none','color', 'g','linewidth',1.5);
plot(linspace(0,3,10),linspace(0,3,10),'k','LineWidth',2);
ylim([1.0 2.5]);
xlabel({'Mean confidence'; 'in LSHC trials'},'FontSize',20);
ylabel({'Mean confidence'; 'in HSLC trials'},'FontSize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(2,2,4)
hold on;
for sub=1:num_sub
    if sub==5 || sub==9
        bar(sub,mean(confidence_diff(sub,:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        hold on;
    else
        bar(sub,mean(confidence_diff(sub,:),2), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
        hold on;
    end
end
hold on;
errorbar(1:num_sub,mean(confidence_diff,2),squeeze(l_conf_diff),squeeze(h_conf_diff),'ok','LineWidth', 1.5)
hold on;
text([1:num_sub]-0.25,squeeze(mean(confidence_diff,2))+0.2,stars_conf,'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Confidence in LSHC ';' - Confidence in HSLC'},'Fontsize',20);
yline(0.0,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
ylim([-0.25 1]);
xlim([0.5 10.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = [1 2 3 4 5 6 7 8 9 10];
hold on;