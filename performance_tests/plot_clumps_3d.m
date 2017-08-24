% plots 3d clumps

close all
figure

filename_list = {...
    '/home/dmitrip/w/inference_without_permutation/labeled_counts.csv',...
    '/home/dmitrip/w/inference_without_permutation/labeled_counts.csv',...
    '/home/dmitrip/w/inference_without_permutation/labeled_counts.csv'};
%     'labeled_counts_1.csv',...
filename_list = {...
    'labeled_counts_2.csv',...
    'labeled_counts_2.csv',...
    'labeled_counts_2.csv'}
colors_ = colormap('lines')';
rand_perm = randperm(size(colors_,2));

rand_perm(1:5) = [64; 26; 55; 35; 3];

for file_num = 1:length(filename_list)
    subplot(1,length(filename_list),file_num);
    hold on
    
    data = csvread(filename_list{file_num});
    points = data(:,2:4);
    mults = data(:,5);
    
    %
    hold on
    % plot3(data(:,2), data(:,3), data(:,4), 'b.')
    u = unique(data(:,1));
    colors = colors_(:,rand_perm(1:length(u)));
    %     colors = rand(3,length(u));
    %     for i = 1:size(colors,2)
    %         [~,r] = max(colors(:,i));
    %         colors(r,i) = 1
    %     end
    
    alpha_ = 0.2;
    for i = 1:length(u)
        i_u = find(data(:,1) == u(i));
        plot3(data(i_u,2), data(i_u,3), data(i_u,4), '.', 'markersize', 20, 'color',colors(:,i));
        for j = 1:length(i_u)
            c = i_u(j);
            plot3(data(c,2), data(c,3), data(c,4), '.',...
                'markersize', 20*mults(c),...
                'color',colors(:,i));
        end
        
        if length(i_u) >= 4
            C = convhull(points(i_u,1), points(i_u,2), points(i_u,3));
            
            for j = 1:size(C,1)
                c = i_u(C(j,:));
                patch(points(c,1), points(c,2), points(c,3), [1,0,0], ...
                    'facecolor', colors(:,i), 'edgecolor', alphacolor(colors(:,i),0.5));
            end
        elseif length(i_u) == 3
            c = i_u;
            patch(points(c,1), points(c,2), points(c,3), [1,0,0], ...
                'facecolor', colors(:,i), 'edgecolor', colors(:,i));
        elseif length(i_u) <= 2
            c = i_u;
            plot3(points(c,1), points(c,2), points(c,3), 'ko-', 'color', colors(:,i),...
                'markersize',1);
            %         elseif length(i_u) == 1
            c = i_u;
            plot3(points(c,1), points(c,2), points(c,3), 'ko', ...
                'markersize', 10, 'color', colors(:,i));
%             if length(i_u) == 1
%                 plot3(points(c,1), points(c,2), points(c,3), 'k', ...
%                     'markersize', 10, 'color', colors(:,i),...
%                     'markerfacecolor', alphacolor(colors(:,i),0.5));
%                 plot3(data(i_u,2), data(i_u,3), data(i_u,4), '.', 'markersize', 20, 'color',colors(:,i));
%             end
        end
    end
    
    xlim_ = get(gca,'xlim');
    ylim_ = get(gca,'ylim');
    zlim_ = get(gca,'zlim');
    
    set(gca,'xlim',[0 xlim_(2)])
    set(gca,'ylim',[0 ylim_(2)])
    set(gca,'zlim',[0 zlim_(2)])
    
    grid minor
    
    camproj('perspective')
    axis vis3d
    grid on
    alpha(alpha_)
    %     set(gca,'color',[1,1,1].*0.75)
    %     box on
    xlabel('$n_1$','fontsize',16,'interpreter','latex')
    ylabel('$n_2$','fontsize',16,'interpreter','latex')
    zlabel('$n_3$','fontsize',16,'interpreter','latex')
    view([40 10])
    zoom(0.9)
    
    if file_num == 1
        view([20 10])
    elseif file_num == 2
        view([110 10])
    end
    
    
    if file_num <= 1
        view([0 0]+5)
    elseif file_num == 2
        view([90 0]+5)
        %         set(gca,'xdir','reverse')
    elseif file_num == 3
        view([0 90])
        set(gca,'ylim',[0 15])
    end
    box on
end

set(gcf,'position',[675         524        1103         450])

set(gcf,'position',[319         538        1469         423])

%%
% from https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(gcf,'units','inches')
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'filename','-dpdf','-r0')

% view([    0.7071    0.7071    0.0000   -0.7071 ; ...
%    -0.2649    0.2649    0.9272   -0.4636 ; ...
%    -0.6556    0.6556   -0.3746    8.8476 ; ...
%          0         0         0    1.0000])