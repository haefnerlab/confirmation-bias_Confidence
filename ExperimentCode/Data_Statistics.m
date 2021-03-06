ps = .51:.02:.99;
params = Model.newModelParams();
[ss, cc] = meshgrid(ps);
meanev = zeros(size(ss));
varev = zeros(size(ss));
meandat = zeros(size(ss));
vardat = zeros(size(ss));
%%
for ii=1:numel(ss)
    params.sensory_info = ss(ii);
    params.category_info = cc(ii);
    params.var_s = Model.getEvidenceVariance(ss(ii));
    params.p_match = cc(ii);
    [data, categories] = Model.genDataWithParams(params);
    data = data .* categories;
    llo = Model.logLikelihoodOdds(params, data);
    meanev(ii) = mean(llo(:));
    varev(ii) = var(llo(:));
    meandat(ii) = mean(data(:));
    vardat(ii) = var(data(:));
end
%%
subplot(2,3,1);
imagesc('XData', ps, 'YData', ps, 'CData', meanev);
set(gca, 'YDir', 'normal');
axis image; colorbar;
title('mean LLO');
subplot(2,3,2);
imagesc('XData', ps, 'YData', ps, 'CData', sqrt(varev));
set(gca, 'YDir', 'normal');
axis image; colorbar;
title('\sigma of LLO');
subplot(2,3,3);
imagesc('XData', ps, 'YData', ps, 'CData', meanev ./ sqrt(varev));
set(gca, 'YDir', 'normal');
axis image; colorbar;
title('\mu/\sigma of LLO');
subplot(2,3,4);
imagesc('XData', ps, 'YData', ps, 'CData', meandat);
set(gca, 'YDir', 'normal');
axis image; colorbar;
title('mean s');
subplot(2,3,5);
imagesc('XData', ps, 'YData', ps, 'CData', sqrt(vardat));
set(gca, 'YDir', 'normal');
axis image; colorbar;
title('\sigma of s');
subplot(2,3,6);
imagesc('XData', ps, 'YData', ps, 'CData', meandat ./ sqrt(vardat));
set(gca, 'YDir', 'normal');
axis image; colorbar;
title('\mu/\sigma of s');