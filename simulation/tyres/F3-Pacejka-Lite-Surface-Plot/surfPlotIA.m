% SURFPLOTIA plots FY vs. FZ vs. SA plot for multiple cambers (IAs). Requires 
% TireDataAnalysis.mlx to have been run. 

clear Xplot Yplot Z

for n=1:3
    [Xplot,Yplot] = meshgrid(S_M{1,n,1}, FZ_binvalues);
    
    for m=1:5
        Z(m,:) = F_M{1,n,m};
    end
    
    Xplot=Xplot(:,2:8:end);
    Yplot=Yplot(:,2:8:end);
    Z = Z(:,2:8:end);
    
    surf(Xplot,Yplot,Z,'EdgeColor',[n==1 n==2 n==3],'FaceColor',[n==1 n==2 n==3],'FaceAlpha',0.4)
    hold on
    clear Z
end

xlabel('Slip angle ({\circ})','Rotation',20)
ylabel('FZ (N)','Rotation',-30)
zlabel('F_y (N)')
title('F_y vs. FZ vs. SA for multiple IAs')
legend('0 deg IA', '2 deg IA', '4 deg IA')
hold off