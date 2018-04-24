%% plot performance of estimators of L1 distance from uniform distribution, known support set size
% before running this script, run get_data_L1_distance_from_uniform_test_v3.m to generate
% data_L1_distance_from_uniform_performance_test_v3.mat

load('data_L1_distance_from_uniform_performance_test_v3')

rmse_vec = squeeze(sqrt(mean((est_mat-L1_distance_from_uniform_list).^2)));
stderr_vec = squeeze(std((est_mat-L1_distance_from_uniform_list).^2))./sqrt(num_trials);

L1_distance_from_uniform_estimators = {'MLE','valiant','PML_approx'};
distribution_labels = {'Uniform','Mix 2 Uniforms','Zipf(-1)','Zipf(-0.6)'};

legend_labels = {'MLE','VV','APML'};
marker_shapes = 's^o*t';
color_list = [0 0 1 ; 0 0.75 0 ; 1 0 0 ; 1 0.5 0.5 ; 0 0 0];
line_width = 1;

% figure
d = 0;
for distr_ix = 1:length(distribution_labels)
    d = d + 1;
    subplot(2, length(distribution_labels)/2, d);
    hold on
    est_num = 0;
    for i = 1:length(L1_distance_from_uniform_estimators)
        est_num = est_num + 1;
        y = rmse_vec(est_num,:,distr_ix);
        y_error_bars = stderr_vec(est_num,:,distr_ix);
        line_obj = errorbar(n_list, y, y_error_bars, 'ko-','color',color_list(est_num,:),...
            'linewidth',line_width,'marker',marker_shapes(i));
    end
    if d == 4
        legend_obj = legend(legend_labels);
        set(legend_obj,'location','northeast');
        set(legend_obj,'box','off');
        set(legend_obj,'FontWeight','bold');
        set(legend_obj,'fontsize',16)
        set(legend_obj,'interpreter','latex');
        set(legend_obj,'position',[0.7370    0.3241    0.1683    0.1257])
    end
    if d <= Inf
        title(distribution_labels{distr_ix},'fontweight','bold','fontsize',12)
    end
    if d >= 0
        xlabel('sample size','fontsize',12)
    end
    set(gca,'xlim',[10,1e8])
    set(gca,'xscale','log')
    set(gca,'fontsize',8)
    set(gca,'yscale','log')
    box on
    grid on
    set(gca,'xminorgrid','off')
    set(gca,'yminorgrid','off')
end