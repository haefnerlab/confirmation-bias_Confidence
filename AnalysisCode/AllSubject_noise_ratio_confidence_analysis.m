clear all; close all; clc;
boots = 500;% number of bootstraps to get PK
thresh_trials_only = 0;
hpr_ridge = logspace(-1, 5, 7);
hpr_ar1 = 0.0;
hpr_curvature = logspace(-1, 5, 7);
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
perf_lo = 0.3;
perf_hi = 0.7;

% disp('Starting to find best hyperparameters for NOISE condition data across subjects....');
% [best_hprs_noise] = [1 0 100];
% disp('Starting to find best hyperparameters for RATIO condition data across subjects....');
% [best_hprs_ratio] = [10 0 10];
%%
for i=1:(num_sub)
    tic;
    for j=1:cases
        phase = 3-j;
        if phase==expt_type_noise
%             best_hprs = best_hprs_noise;
            disp(['Starting analysis for NOISE data of Subject ' num2str(i) ' ...']);
        elseif phase==expt_type_ratio
%             best_hprs = best_hprs_ratio;
            disp(['Starting analysis for RATIO data of Subject ' num2str(i) ' ...']);
        end
        [conf_analysis{i}{j},thresh(i,j,:),count_trials_thresh(i,j)] = run_confidence_regression_both_ratio_and_noise(subjects{i},phase,perf_lo,perf_hi,boots,hpr_ridge,hpr_ar1,hpr_curvature,standardize,folds,dir,thresh_trials_only);
        for cs=1:6
            temporal_kernel(cs,i,j,:) = prctile(conf_analysis{i}{j}{cs}.params_boot(:, 1:num_frames), 50);
            norm_temporal_kernel(cs,i,j,:) = temporal_kernel(cs,i,j,:)/mean(temporal_kernel(cs,i,j,:));
        end
    end
    disp(['All analysis complete for Subject ' num2str(i) ' !!!!']);
    toc;
    disp('-----------------------------------------------------------------------------------------------------');
end


%%
f = figure();
set(f,'defaultLegendAutoUpdate','off');
for cs=1:6
    subplot(2,3,cs);
    for i=1:(num_sub)
        plot(1:num_frames,squeeze(norm_temporal_kernel(cs,i,1,1:num_frames)),'-or','LineWidth',0.2);
        cs1_sub(i,:) = squeeze(norm_temporal_kernel(cs,i,1,1:num_frames));
        xlabel('Frames');
        ylabel('Weights');
        hold('on');
        plot(1:num_frames,squeeze(norm_temporal_kernel(cs,i,2,1:num_frames)),'-ob','LineWidth',0.2);
        cs2_sub(i,:) = squeeze(norm_temporal_kernel(cs,i,2,1:num_frames));
        legend({['noise trials'],['ratio trials']});
        title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
        hold('on');
        axis('tight');
    end
    yline(0.0,'k','linewidth',2);
%         plot(1:num_frames,mean(cs1_sub,1),'-or','LineWidth',2);
%         plot(1:num_frames,mean(cs2_sub,1),'-ob','LineWidth',2);
%     ylim([-10 10]);
    title(['Weights for ' conf_analysis{1}{1}{cs}.case ' comparison']);
end
% sgtitle('All types of temporal weights fit across subjects','fontsize',30);

%%

for cs=1:6
    f = figure(cs+1);
    set(f,'defaultLegendAutoUpdate','off');
    
    for i=1:(num_sub)
        subplot(2,5,i);
        plot(1:num_frames,squeeze(norm_temporal_kernel(cs,i,1,1:num_frames)),'-or','LineWidth',2);
        cs1_sub(i,:) = squeeze(norm_temporal_kernel(cs,i,1,1:num_frames));
        xlabel('Frames');
        ylabel('Weights');
        hold('on');
        plot(1:num_frames,squeeze(norm_temporal_kernel(cs,i,2,1:num_frames)),'-ob','LineWidth',2);
        cs2_sub(i,:) = squeeze(norm_temporal_kernel(cs,i,2,1:num_frames));
        legend({['noise trials'],['ratio trials']});
        title(['Noise Vs Ratio Stimulus for Subject ' num2str(i)])
        hold('on');
        axis('tight');
        ylim([-2 5]);
        yline(0.0,'k','linewidth',2);
    end
    
    %     plot(1:num_frames,mean(cs1_sub,1),'-or','LineWidth',4);
    %     plot(1:num_frames,mean(cs2_sub,1),'-ob','LineWidth',4);
    
    sgtitle(['Weights for ' conf_analysis{1}{1}{cs}.case ' comparison']);
end
