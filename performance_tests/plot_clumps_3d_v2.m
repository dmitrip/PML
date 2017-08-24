% plots 3d clumps

close all
figure

filename_list = {...
    '/home/dmitrip/w/inference_without_permutation/labeled_counts.csv',...
    '/home/dmitrip/w/inference_without_permutation/labeled_counts.csv',...
    '/home/dmitrip/w/inference_without_permutation/labeled_counts.csv'};

colors_ = colormap('lines')';
rand_perm = randperm(size(colors_,2));

for file_num = 1:length(filename_list)
    subplot(2,2,file_num);
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
    xlabel('$n_1$','fontsize',18,'interpreter','latex')
    ylabel('$n_2$','fontsize',18,'interpreter','latex')
    zlabel('$n_3$','fontsize',18,'interpreter','latex')
    view([40 10])
%     zoom(0.85)
    
    if file_num == 1
        view([0 0] + 5)
    elseif file_num == 2
        view([90 0] + 5)
    elseif file_num == 3
        view([0 90])
    end
    
    
%     if file_num <= 1
%         view([0 0]+5)
%     elseif file_num == 2
%         view([90 0]+5)
%     elseif file_num == 3
%         view([0 90]+5)
%     end
end

subplot(2,2,4)
hold on

S = 100;
p_1 = ones(S,1)./S;
p_2 = ones(S,1)./S;
p_3 = [(5/(2*S)).*ones(S/5,1) ; (5/(8*S)).*ones(4*S/5,1)];

plot(-1,-1,'b.','markersize',20)
plot(-2,-2,'r.','markersize',20)
plot(-3,-3,'g.','color',[0 1 0].*0.5,'markersize',20)

plot(p_1, 'b-','linewidth',4)
plot(p_2, 'r-','linewidth',2)
plot(p_3, 'g.','color',[0 1 0].*0.5)
set(gca,'xlim',[1 S])
set(gca,'ylim',[0 max(p_3)*1.2])
xlabel('$x$','interpreter','latex','fontsize',18)
legend_obj = legend(...
    '$p^{(1)}_x$ (Uniform)',...
    '$p^{(2)}_x$ (Uniform)',...
    '$p^{(3)}_x$ (Mix 2 Unif.)');
set(legend_obj, 'interpreter', 'latex', 'fontsize',16)

set(gcf,'position',[675         524        1103         450])
set(gcf,'position',[961     1   960   973])
set(gcf,'position',[956   -32   970   948])

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