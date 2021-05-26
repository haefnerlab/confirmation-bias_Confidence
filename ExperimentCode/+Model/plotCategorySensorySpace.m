% function [correct, fig, fig1, fig2, opt_fig] = plotCategorySensorySpace(category_infos, sensory_infos, params, optimize, optim_grid_size)
function [correct, pc_fig, conf_fig, opt_fig] = plotCategorySensorySpace(category_infos, sensory_infos, params, optimize, optim_grid_size)
%PLOTCATEGORYSENSORYSPACE make category_info vs sensory_info plots for the
%given params.
switch lower(params.model)
    case 'is_temporal'
        model_text = 'temporal';
    case 'is_spatial'
        model_text = 'spatial';
    case 'vb'
        model_text = 'variational';
    case 'vb-czx'
        model_text = 'variational-czx';
    case 'ideal'
        model_text = 'ideal';
end

savedir = fullfile('+Model', 'figures');
if ~exist(savedir, 'dir'), mkdir(savedir); end

% savedir1 = fullfile('+Model', 'saved results');
% if ~exist(savedir1, 'dir'), mkdir(savedir1); end

if nargin < 4, optimize = {}; end
if nargin < 5, optim_grid_size = 11; end

if exist('optimize', 'var') && ~isempty(optimize)
    optim_prefix = Model.getOptimPrefix(optimize, optim_grid_size);
    results_uid = Model.getModelStringID(params, true);
    optim_results_uid = ['optim_CS_' optim_prefix '_' results_uid];
    [optim_params, ~, ~] = LoadOrRun(@Model.optimizeParamsOverCS, ...
        {category_infos, sensory_infos, params, optimize, optim_grid_size}, ...
        fullfile(params.save_dir, optim_results_uid));
else
    optimize = {};
    optim_params = [];
end

if strcmpi(params.model, 'ideal') && ~isempty(optimize)
    error('Nothing to optimize for the ideal observer');
end

[ss, cc] = meshgrid(sensory_infos, category_infos);

% Preallocate return values.
correct = nan(size(ss));
confidence = nan(size(ss));
samples_var = nan(size(ss));

for i=1:numel(ss)
    params_copy = params;
    % Set data-generating parameters.
    params_copy.sensory_info = ss(i);
    params_copy.category_info = cc(i);
    % Set variances for this pair of category- & sensory-info values. (That is, assume that the
    % model 'knows' the environment statistics)
    params_copy.var_s = Model.getEvidenceVariance(ss(i));
    %     params_copy.p_match = cc(i);
    if params.exact_category_info
        params_copy.p_match = cc(i);
    else
        temp = random('normal', cc(i), 0.05, 1, 1);%
        if temp<0.51
            temp = 0.51;
        elseif temp>0.99
            temp = 0.99;
        end
        params_copy.p_match = temp;
    end
    
    % TODO - smarter setting of seed?
    params_copy.seed = randi(1000000000);
    
    if isempty(optimize)
        % Run the model on the given params.
        results_uid = Model.getModelStringID(params_copy);
        results = LoadOrRun(@Model.runVectorized, {params_copy}, ...
            fullfile(params.save_dir, results_uid));
    else
        % Run the model using the best params at this value of category and sensory info.
        for iVar=1:length(optimize)
            params_copy.(optimize{iVar}) = optim_params(i).(optimize{iVar});
        end
        results_uid = Model.getModelStringID(params_copy);
        results = LoadOrRun(@Model.runVectorized, {params_copy}, ...
            fullfile(params.save_dir, results_uid), '-verbose');
    end
    
    [~, correct_categories] = Model.genDataWithParams(results.params);
    
    correct(i) = mean(results.choices == correct_categories);
    confidence(i) = meanci(abs(results.lpo(:, end)), 0.67);
    temp_samples = reshape(results.samples,params.trials,params.samples*params.frames*params.updates);
    samples_var(i) = meanci(var(temp_samples,1,2),0.67);
    if mod(i,10)==1, disp(i); end
end
% confidence_norm = (confidence-min(confidence(:)))./(max(confidence(:))-min(confidence(:)));

% Plot percent correct
pc_fig = figure();
imagesc(correct, [0.5 1.0]); axis image; temp1 = colorbar('FontSize',20); temp1.Label.String = 'Percent correct';
% Add contour line at threshold
hold on; contour(imfilter(correct, smoothkernel(9), 'replicate'), [0.7 0.7], '-w', 'LineWidth', 2);
category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 3)));%5)));
sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 3)));%5)));
set(gca, 'YTick', category_tick_indices);
set(gca, 'XTick', sensory_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('Sensory Information');
ylabel('Category Information');
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
title(['Performance of ' model_text ' model'],'FontSize', 20);

% figname = ['CSSpace_' Model.getModelStringID(params, true) '.fig'];
% saveas(gcf, fullfile(savedir, figname));
% filename = ['CSSpace_' Model.getModelStringID(params, true) '.mat'];
% saveas(fullfile(savedir1, filename));


% Plot confidence
conf_fig = figure();
% imagesc(confidence,[0.0 max(confidence(:))]); axis image; colorbar;
imagesc(sqrt(confidence),[0.0 max(sqrt(confidence(:)))]); axis image; temp = colorbar('FontSize',20); temp.Label.String = 'sqrt(Mean of absolute LPO)';%set(gca,'ColorScale','log');
% Add contour line at threshold
hold on; contour(imfilter(correct, smoothkernel(9), 'replicate'), [0.7 0.7], '-w', 'LineWidth', 2);
category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 3)));%5)));
sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 3)));%5)));
set(gca, 'YTick', category_tick_indices);
set(gca, 'XTick', sensory_tick_indices);
set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
set(gca, 'YDir', 'Normal');
xlabel('Sensory Information');
ylabel('Category Information');
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
title(['Confidence as mean of absolute LPO in ' model_text ' model'],'FontSize', 20);

% figname = ['CSSpaceLPO_' Model.getModelStringID(params, true) '.fig'];
% saveas(gcf, fullfile(savedir, figname));

% Plot sample variance
% fig2 = figure();
% % imagesc(confidence,[0.0 max(confidence(:))]); axis image; colorbar;
% imagesc(1.0./samples_var,[(1.0./max(samples_var(:))) (1.0./min(samples_var(:)))]); axis image; colorbar;
% % Add contour line at threshold
% hold on; contour(imfilter(correct, smoothkernel(9), 'replicate'), [0.7 0.7], '-w', 'LineWidth', 2);
% category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
% sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
% set(gca, 'YTick', category_tick_indices);
% set(gca, 'XTick', sensory_tick_indices);
% set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
% set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
% set(gca, 'YDir', 'Normal');
% xlabel('SI');
% ylabel('CI');
% title('Certainty');
% 
% figname = ['CSSpaceSamplesVar_' Model.getModelStringID(params, true) '.fig'];
% saveas(gcf, fullfile(savedir, figname));

% filename = ['CSSpaceLPO_' Model.getModelStringID(params, true) '.mat'];
% saveas(fullfile(savedir1, filename));


% fig2 = figure();
% scatter(correct(:),confidence(:));
% xlabel('Accuracy');
% ylabel('Confidence');
% title('Accuracy vs Confidence');
%
% figname = ['CSSpaceLPOvsCorrect_' Model.getModelStringID(params, true) '.fig'];
%
% saveas(gcf, fullfile(savedir, figname));
%

% Plot value of optimized parameters.
for i=length(optimize):-1:1
    opt_fig(i) = figure();
    
    % Unravel optimal param values
    opt_param_values = reshape([optim_params.(optimize{i})]', size(ss));
    
    imagesc(opt_param_values); axis image; colorbar; colormap cool; 
    category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
    sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
    set(gca, 'YTick', category_tick_indices);
    set(gca, 'XTick', sensory_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('Sensory Information');
    ylabel('Category Information');
    title(['Optimized value of ' optimize{i}]);
    figname = ['CSSpace_optim_' optimize{i} '_' Model.getModelStringID(params, true) '.fig'];
    saveas(gcf, fullfile(savedir, figname));
end
end

function kernel = smoothkernel(n)
x = linspace(-2,2,n);
[xx,yy] = meshgrid(x,x);
kernel = exp(-(xx.^2 + yy.^2));
kernel = kernel / sum(kernel(:));
end