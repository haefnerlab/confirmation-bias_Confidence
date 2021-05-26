function [best_hprs] = combined_hprs_search(subjectID, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir)
% initialize
datadir = fullfile(pwd, dir);
[num_sub,~] = size(subjectID);

for sub =1:num_sub
    data = LoadAllSubjectData(subjectID{sub},expt_type,datadir);
    signal_raw = [];
    choice_raw = [];
    for k = 1:size(data.ideal_frame_signals, 1)
        signal_raw = [signal_raw; data.ideal_frame_signals(k, :)];
        choice_raw = [choice_raw data.choice(k)];
    end
    disp(['Searching hyperparameters for Subject ' num2str(sub)]);
    [~, log_likelihoods(sub,:,:,:,:)] = CustomRegression.xValidatePKwithlapseSabya(signal_raw, choice_raw, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds);
    log_likelihood_summed(sub,:,:,:) = mean(log_likelihoods(sub,:,:,:,:),5);
    disp(['Done searching for ' num2str(sub) '/' num2str(num_sub) ' subjects...']);
end
sz = size(log_likelihoods);
avg_ll =  mean(log_likelihood_summed,1);
[~, imax] = max(avg_ll(:));
[iRidge, iAR1, iCurve] = ind2sub(sz(2:4), imax);
% Err on the side of less regularization by choosing smoothing that is one order of magnitude less than the best.
% iRidge = max(iRidge-1, 1);
% iAR1 = max(iAR1-1, 1);
% iCurve = max(iCurve-1, 1);
best_hprs = [hpr_ridge(iRidge), hpr_ar1(iAR1), hpr_curvature(iCurve)];

disp('Searching of hyperparameters complete!!');
disp (['Best hyperparameters across subjects is: ' num2str(best_hprs)]);

end