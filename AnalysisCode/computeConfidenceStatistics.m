function [subj_resp,subj_resp_err,confidence,confidence_err,ntrial_subj, confidence_hist,confidence_rt,choice_rt,direct_acc_conf] = computeConfidenceStatistics(data,exp)
c = data.conf(:);
c1 = sum(c==1);
r1 = mean(data.conf_reaction(data.conf==1));
r1_err = std(data.conf_reaction(data.conf==1))/sqrt(c1);
r11 = mean(data.reaction_time(data.conf==1));
r11_err = std(data.reaction_time(data.conf==1))/sqrt(c1);
r111 = sum(data.accuracy(data.conf==1))/c1;
r111_err = sqrt(r111 * (1 - r111))/sqrt(c1);
c2 = sum(c==2);
r2 = mean(data.conf_reaction(data.conf==2));
r2_err = std(data.conf_reaction(data.conf==2))/sqrt(c2);
r22 = mean(data.reaction_time(data.conf==2));
r22_err = std(data.reaction_time(data.conf==2))/sqrt(c2);
r222 = sum(data.accuracy(data.conf==2))/c2;
r222_err = sqrt(r222 * (1 - r222))/sqrt(c2);
c3 = sum(c==3);
r3 = mean(data.conf_reaction(data.conf==3));
r3_err = std(data.conf_reaction(data.conf==3))/sqrt(c3);
r33 = mean(data.reaction_time(data.conf==3));
r33_err = std(data.reaction_time(data.conf==3))/sqrt(c3);
r333 = sum(data.accuracy(data.conf==3))/c3;
r333_err = sqrt(r333 * (1 - r333))/sqrt(c3);
confidence_hist = [c1 c2 c3]./(c1 + c2 + c3);
confidence_rt = [r1 r2 r3; r1_err r2_err r3_err];
choice_rt = [r11 r22 r33; r11_err r22_err r33_err];
direct_acc_conf = [r111 r222 r333; r111_err r222_err r333_err];
bins = 4;
if exp==2
    uniq_vals = linspace(-0.8,0.8,bins);
    tr_kappa = data.sign_noise;
    for tt=1:(length(uniq_vals)-1)
        subj_resp(tt) = mean(data.accuracy(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
        confidence(tt) = mean(data.conf(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
        ntrial_subj(tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
        confidence_err(tt) = std(data.conf(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))/sqrt(ntrial_subj(tt));
    end
    subj_resp_err = squeeze((subj_resp(:)).*(1-subj_resp(:))./sqrt(ntrial_subj(:)));
    
else
    uniq_vals = linspace(0,1,bins);
    tr_ratio = data.true_ratio;
    for tt=1:(length(uniq_vals)-1)
        subj_resp(tt) = mean(data.accuracy(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
        confidence(tt) = mean(data.conf(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
        ntrial_subj(tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
        confidence_err(tt) = std(data.conf(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)))/sqrt(ntrial_subj(tt));    
    end
    subj_resp_err = squeeze((subj_resp(:)).*(1-subj_resp(:))./sqrt(ntrial_subj(:)));
    
    
end