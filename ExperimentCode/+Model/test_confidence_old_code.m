
diff = 0.02;
category_infos = 0.51:diff:0.99;
sensory_infos = 0.51:diff:0.99;
params = Model.newModelParams();
samples = [5];
ns = length(samples);
category_info_thresh_all = zeros(length(samples));
sensory_info_thresh_all = zeros(length(samples));

for i=1:length(samples)
    disp(i);
    params.samples = samples(i);
    [correct(i,:,:), fig(i), cont{i}, tmp_smp{i}, lpo{i}] = Model.plotCategorySensorySpace(category_infos, sensory_infos, params);%, { 'gamma'}, 10);
end

%%
for i=1:ns
    [sen, cat, tr, fr, s, up] = size(tmp_smp{i});
    Samples_obtained{i} = squeeze(reshape(squeeze(tmp_smp{i}),sen, cat, tr, fr * s * up));
    confidence(i,:,:) = mean(abs(lpo{i}),3);
%     certainty(i,:,:) = mean(-sqrt(var(Samples_obtained{i},0,4)/(fr * s * up)),3);
end
confidence = (confidence-min(confidence(:)))./(max(confidence(:))-min(confidence(:)));
% certainty = (certainty-min(certainty(:)))./(max(certainty(:))-min(certainty(:)));
disp('Confidence and Certainty Computations Done!!');
%% Plots for Confidence and Accuracy for all ns and collapsing CI and SI
disp('Preparing figure!!');
figure();
hold on;
for i=1:ns
    subplot(2,ns,i)
    hold on
    temp1 = confidence(i,:,:);
    temp2 = correct(i,:,:);
    temp = corrcoef(temp1(:),temp2(:));
    R(i) = temp(1,2);
    plot(temp1(:),temp2(:),'bx');
    xlabel('Confidence across all CI and SI');
    ylabel('Accuracy');
%     text(0.1,0.45,['Corr = ', num2str(round(R(i),2))]);
%     xlim([0 1]);
%     ylim([0.4 1]);
    title({['Confidence for samples = ' num2str(samples(i))],['Corr = ', num2str(round(R(i),2))]})
end
hold on;
for i=1:ns
    hold on;
    subplot(2,ns,ns+i)
    temp1 = certainty(i,:,:);
    temp2 = correct(i,:,:);
    temp = corrcoef(temp1(:),temp2(:));
    R2(i) = temp(1,2);
    plot(temp1(:),temp2(:), 'ro');
    xlabel('Confidence across all CI and SI');
    ylabel('Accuracy');
%     xlim([0 1]);
%     ylim([0.4 1]);
%     text(0.4,0.3,['Corr = ', num2str(round(R2(i),2))]);
    title({['Certainty for samples = ' num2str(samples(i))],['Corr = ', num2str(round(R2(i),2))]})
end
suptitle('Confidence and Certainty Measures Versus Accuracy')
%% Plots for Confidence and Accuracy for all CI and SI and ns
disp('Preparing figure!!');
figure()
hold on;
k = 1;
for i=1:ns
%     subaxis(ns,2, k, 'sh', 0.0, 'sv', 0.05, 'padding', 0.0, 'margin', 0.1);
    subplot(ns,3,k)
    imagesc(squeeze(confidence(i,:,:))); caxis manual; caxis([0 1]);colorbar; axis image
    category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
    sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
    set(gca, 'YTick', category_tick_indices);
    set(gca, 'XTick', sensory_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('SI');
    ylabel('CI');
    title(['Confidence for x samples = ' num2str(samples(i))])
    
    k = k + 1;
%     subaxis(ns,2, k, 'sh', 0.0, 'sv', 0.05, 'padding', 0.0, 'margin', 0.1);
    subplot(ns,3,k)
    imagesc(squeeze(correct(i,:,:)));caxis manual; caxis([0 1]);colorbar; axis image;
    category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
    sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
    set(gca, 'YTick', category_tick_indices);
    set(gca, 'XTick', sensory_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('SI');
    ylabel('CI');
    title(['Accuracy for x samples = ' num2str(samples(i))])
    
    k = k + 1;
   
    subplot(ns,3,k)
    imagesc(squeeze(certainty(i,:,:))); caxis manual; caxis([0 1]);colorbar; axis image
    category_tick_indices = round(linspace(1, length(category_infos), min(length(category_infos), 5)));
    sensory_tick_indices = round(linspace(1, length(sensory_infos), min(length(sensory_infos), 5)));
    set(gca, 'YTick', category_tick_indices);
    set(gca, 'XTick', sensory_tick_indices);
    set(gca, 'YTickLabel', arrayfun(@num2str, category_infos(category_tick_indices), 'UniformOutput', false));
    set(gca, 'XTickLabel', arrayfun(@num2str, sensory_infos(sensory_tick_indices), 'UniformOutput', false));
    set(gca, 'YDir', 'Normal');
    xlabel('SI');
    ylabel('CI');
    title(['Certainty for x samples = ' num2str(samples(i))])
    k = k + 1;
    hold on;
end
suptitle('Confidence and Certainty Measures Versus Accuracy')



%% Plots for Confidence and ns for fixed accuracy of 70%

figure();
l = 0;
for i=1:ns
    ind = cont{i};
    ind_x = round(ind(1,2:end-1));
    ind_y = round(ind(2,2:end-1));
    confidence_fixed_accuracy{i} = squeeze(confidence(i,ind_x,ind_y));
    confidence_fixed_accuracy_mat(1,l+1:l+length(ind_x)*length(ind_y)) = (reshape(confidence(i,ind_x,ind_y),1,length(ind_x)*length(ind_y)));
    ns_mat(1,l+1:l+length(ind_x)*length(ind_y)) = samples(i);
    l = length(confidence_fixed_accuracy_mat);
    
    confidence_noise_condition(i) = confidence(i,round(ind_x(1)),round(ind_y(1)));
    confidence_ratio_condition(i) = confidence(i,round(ind_x(end)),round(ind_y(end)));
    
    
end
subplot(1,2,1)
plot(samples,confidence_noise_condition,'-o','LineWidth',2);
ylabel('Confidence for 70% performance');
xlabel('Number of x Samples per update (drawn in parallel)');
title('Confidence as LPO of C (normalized)')
hold on;
plot(samples,confidence_ratio_condition,'-o','LineWidth',2);
legend({'noise condition'; 'ratio condition'})
hold on;


for i=1:ns
    ind = cont{i};
    ind_x = round(ind(1,2:end-1));
    ind_y = round(ind(2,2:end-1));
    certainty_fixed_accuracy{i} = squeeze(certainty(i,ind_x,ind_y));
    certainty_fixed_accuracy_mat(1,l+1:l+length(ind_x)*length(ind_y)) = (reshape(certainty(i,ind_x,ind_y),1,length(ind_x)*length(ind_y)));
    ns_mat2(1,l+1:l+length(ind_x)*length(ind_y)) = samples(i);
    l = length(certainty_fixed_accuracy_mat);
    
    certainty_noise_condition(i) = certainty(i,round(ind_x(1)),round(ind_y(1)));
    certainty_ratio_condition(i) = certainty(i,round(ind_x(end)),round(ind_y(end)));
    
    
end

subplot(1,2,2)
plot(samples,certainty_noise_condition,'-o','LineWidth',2);
ylabel('Certainty for 70% performance');
xlabel('Number of x Samples per update (drawn in parallel)');
title('Certainty as -std of samples of x (normalized)')
hold on;
plot(samples,certainty_ratio_condition,'-o','LineWidth',2);
legend({'noise condition'; 'ratio condition'})
% suptitle('Confidence as a proxy for number of samples');