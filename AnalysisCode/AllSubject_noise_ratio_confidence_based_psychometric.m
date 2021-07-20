clear all; close all; clc;
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
clr = {'r', 'g', 'b'};
for i=1:(num_sub)
    tic;
    f = figure();
    set(f,'defaultLegendAutoUpdate','off');
    for j=1:cases
        phase = 3-j;
        if phase==expt_type_noise
            disp(['Starting analysis for NOISE data of Subject ' num2str(i) ' ...']);
        elseif phase==expt_type_ratio
            disp(['Starting analysis for RATIO data of Subject ' num2str(i) ' ...']);
        end
        [pm_fit{i}{j}, uniq_vals_ps{i}{j}, yvals{i}{j}, stderrs{i}{j},...
            data1{i}{j}, data2{i}{j}, data3{i}{j},...
            num_trials1{i}{j}, num_trials2{i}{j}, num_trials3{i}{j}]...
            = compute_psychometric(subjects{i},phase,dir);
        subplot(1,2,j)
        for k =1:3
            if k==1
                data = data1{i}{j};
                txt{1} = ['conf=1 trials=' num2str(num_trials1{i}{j})];
            elseif k==2
                data = data2{i}{j};
                txt{2} = ['conf=2 trials=' num2str(num_trials2{i}{j})];
            else
                data = data3{i}{j};
                txt{3} = ['conf=3 trials=' num2str(num_trials3{i}{j})];
            end
            if j==1
                bins = 10;
                subj_resp = [];
                uniq_vals = [];
                ntrial_subj = [];
                uniq_vals = linspace(-0.8,0.8,bins);
                tr_kappa = data.sign_noise;
                for tt=1:(length(uniq_vals)-1)
                    subj_resp(i,j,tt) = mean(data.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
                    ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
                end
                noise_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
                errorbar(noise_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'ok','MarkerFaceColor',clr{k},'linewidth',1);
                hold on;
            else
                bins = 10;
                subj_resp = [];
                uniq_vals = [];
                ntrial_subj = [];
                subject_pm_curve = [];
                uniq_vals = linspace(0,1,bins);
                tr_ratio = data.true_ratio;
                for tt=1:(length(uniq_vals)-1)
                    subj_resp(i,j,tt) = mean(data.choice(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
                    ntrial_subj(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
                end
                ratio_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
                errorbar(ratio_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'ok','MarkerFaceColor',clr{k},'linewidth',1);
                hold on;
            end
            hold on;
            subject_pm_curve_psig{i}{j}{k} = (1-pm_fit{i}{j}{k}.Fit(3)-pm_fit{i}{j}{k}.Fit(4))*arrayfun(@(x) pm_fit{i}{j}{k}.options.sigmoidHandle(x,pm_fit{i}{j}{k}.Fit(1),pm_fit{i}{j}{k}.Fit(2)), uniq_vals_ps{i}{j}{k})+pm_fit{i}{j}{k}.Fit(4);
            LH(k)=plot(uniq_vals_ps{i}{j}{k}, subject_pm_curve_psig{i}{j}{k},clr{k},'linewidth',2);
            L{k} = txt{k};
            hold on;
        end
        yline(0.5,'--k', 'linewidth',2);
        legend(LH,L,'box','off','fontsize',15,'location','northwest');
        ylabel('Percent chose +45', 'fontsize',20);
        if j==1
            xline(0.0,'--k', 'linewidth',2);
            xlabel('Energy favoring +45', 'fontsize',20);
            title('LSHC condition', 'fontsize',20);
            xlim([-0.8 0.8]);
        else
            xline(0.5,'--k', 'linewidth',2);
            title('HSLC condition', 'fontsize',20);
            xlabel('Ratio of frames favoring +45', 'fontsize',20);
            xlim([0.0 1.0]);
        end
        ylim([0 1]);
        ax = gca;
        set(ax, 'box','off');
        ax.LineWidth=2;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        hold on;
    end
%     pause;
end
