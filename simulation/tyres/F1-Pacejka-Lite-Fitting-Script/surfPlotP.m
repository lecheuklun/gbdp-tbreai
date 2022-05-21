% SURFPLOTP plots FY vs. FZ vs. SA plot for multiple pressures (Ps). Requires 
% TireDataAnalysis.mlx to have been run. 

clear Xplot Yplot Z

n=1;
for i=1:3
    [Xplot,Yplot] = meshgrid(S_M{i,n,1}, FZ_binvalues);
    
    for m=1:5
        Z(m,:) = F_M{i,n,m};
    end
    
    Xplot=Xplot(:,2:8:end);
    Yplot=Yplot(:,2:8:end);
    Z = Z(:,2:8:end);
    
    surf(Xplot,Yplot,Z,'EdgeColor',[i==1 i==2 i==3],'FaceColor',[i==1 i==2 i==3],'FaceAlpha',0.4)
    hold on
    clear Z
end

xlabel('Slip angle ({\circ})','Rotation',20)
ylabel('FZ (N)','Rotation',-30)
zlabel('F_y (N)')
title('F_y vs. FZ vs. SA for multiple pressures')
legend('10 psi', '12 psi', '14 psi')
hold off