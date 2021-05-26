function name = getModelStringID(params, drop_cs_terms)
switch lower(params.model)
    case 'is_temporal'
        name = sprintf('is_temporal%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_ns%d_noise%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, params.samples, ...
            params.noise, params.lapse);
        if ~params.importance_norm
            name = [name '_unnorm'];
        end
        if ~params.exact_category_info
            name = [name '_pnotmatched'];
        end
    case 'is_spatial'
        name = sprintf('is_spatial_cog%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_ns%d_noise%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, params.samples, ...
            params.noise, params.lapse);
        if ~params.importance_norm
            name = [name '_unnorm'];
        end
        if ~params.exact_category_info
            name = [name '_pnotmatched'];
        end
    case 'vb'
        name = sprintf('vb_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_step%.2e_noise%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, ...
            params.step_size, params.noise, params.lapse);
    case 'vb-czx'
        name = sprintf('vb_cxz_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f_vx%.2f_gam%.2f_nu%d_step%.2e_noise%.2e_lapse%.2e', ...
            params.trials, params.frames, params.category_info, params.sensory_info, params.var_s, ...
            params.p_match, params.prior_C, params.var_x, params.gamma, params.updates, ...
            params.step_size, params.noise, params.lapse);
    case 'ideal'
        name = sprintf('ideal_%dx%d_cinfo%.3f_sinfo%.3f_vs%.2f_pm%.2f_pC%.2f', params.trials, params.frames, ...
            params.category_info, params.sensory_info, params.var_s, params.p_match, params.prior_C);
    otherwise
        error('Unrecognized model type: %s', params.model);
end

if nargin >= 2 && drop_cs_terms
    name = regexprep(name, '_cinfo[0-9.]+_sinfo[0-9.]+_vs[0-9.]+_pm[0-9.]+', '');
end

end