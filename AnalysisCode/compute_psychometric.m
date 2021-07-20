function [pm_fit, uniq_vals_ps, yvals, stderrs,...
    data1, data2, data3, num_trials1, num_trials2, num_trials3] = compute_psychometric(subjectID,expt_type,dir)
datadir = fullfile(pwd, dir);
data = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded for the subject! ');
data.conf = data.conf + 1; % because values are 0, 1 and 2 for pressing keys 1, 2 and 3. Convert confidence values to 1, 2 and 3.

num_trials1 = 0;
num_trials2 = 0;
num_trials3 = 0;
k1 = 1;
k2 = 1;
k3 = 1;
for k = 1:size(data.ideal_frame_signals, 1)
    if data.conf(k)==1
        data1.sign_noise(k1) = data.sign_noise(k);
        data1.choice(k1) = data.choice(k);
        data1.true_ratio(k1) = data.true_ratio(k);
        data1.accuracy(k1) = data.accuracy(k);
        data1.ideal_answer(k1) = data.ideal_answer(k);
        data1.correct_answer(k1) = data.correct_answer(k);
        data1.ratio(k1) = data.ratio(k);
        data1.ideal_frame_signals(k1,:) = data.ideal_frame_signals(k,:);
        num_trials1 = num_trials1 + 1;
        k1 = k1 + 1;
    elseif data.conf(k)==2
        data2.sign_noise(k2) = data.sign_noise(k);
        data2.choice(k2) = data.choice(k);
        data2.true_ratio(k2) = data.true_ratio(k);
        data2.accuracy(k2) = data.accuracy(k);
        data2.ideal_answer(k2) = data.ideal_answer(k);
        data2.correct_answer(k2) = data.correct_answer(k);
        data2.ratio(k2) = data.ratio(k);
        data2.ideal_frame_signals(k2,:) = data.ideal_frame_signals(k,:);
        num_trials2 = num_trials2 + 1;
        k2 = k2 + 1;
    elseif data.conf(k)==3
        data3.sign_noise(k3) = data.sign_noise(k);
        data3.choice(k3) = data.choice(k);
        data3.true_ratio(k3) = data.true_ratio(k);
        data3.accuracy(k3) = data.accuracy(k);
        data3.ideal_answer(k3) = data.ideal_answer(k);
        data3.correct_answer(k3) = data.correct_answer(k);
        data3.ratio(k3) = data.ratio(k);
        data3.ideal_frame_signals(k3,:) = data.ideal_frame_signals(k,:);
        num_trials3 = num_trials3 + 1;
        k3 = k3 + 1;
    end
end

if expt_type==1
    ep = 1;
    
elseif expt_type==2
    ep = -2;
end
for possible_conf=1:3
    if possible_conf==1
        data_use = data1;
    elseif possible_conf==2
        data_use = data2;
    elseif  possible_conf==3
        data_use = data3;
    end
    
    [pm_fit{possible_conf}, uniq_vals_ps{possible_conf}, yvals{possible_conf}, stderrs{possible_conf}] = GaborPsychometric(data_use, ep);
end

end