function [cf] = fit_logistic(Data,type)
if type==2
    uniq_vals = unique(Data.sign_noise);
    init_val = [0 1 0 0.2 0];
    yvals = arrayfun(@(u) mean(Data.choice(Data.sign_noise == u) == +1), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(Data.sign_noise == u), uniq_vals);
    stderrs = arrayfun(@(u) std(Data.choice(Data.sign_noise == u) == +1), uniq_vals) ./ sqrt(num_trials_at_vals);
else
    init_val = [0 1 0.5 0.2 1];
    uniq_vals = unique(Data.true_ratio);
    yvals = arrayfun(@(u) mean(Data.choice(Data.true_ratio == u)), uniq_vals);
    num_trials_at_vals = arrayfun(@(u) sum(Data.true_ratio == u), uniq_vals);
    stderrs = arrayfun(@(u) std(Data.choice(Data.true_ratio == u)), uniq_vals) ./ sqrt(num_trials_at_vals);
end
    function lse_use = compute_lse(params)
        lse_use = lse_val(uniq_vals,yvals,params);
    end
options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e5,'FiniteDifferenceStepSize', 1e-3);
try
[cf,~] = fminunc(@compute_lse, rand(5,1), options);
catch
    [cf,~] = fminunc(@compute_lse, rand(5,1), options);
end

    function lse = lse_val(x,y,params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        E = params(5);
        lse = sum((D+(A-D)/((1+(x/C)^B)^E) - y).^2);
    end

end