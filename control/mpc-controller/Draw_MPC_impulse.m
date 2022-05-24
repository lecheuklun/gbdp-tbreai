%   ______________________________________________________________________
%   **********************************************************************
%   Function uses code proposed in:
%   MMehrez,MPC-and-MHE-implementation-in-MATLAB-using-Casadi,(2021),GitHub
%   https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
%   ______________________________________________________________________
%   **********************************************************************

function Draw_MPC_impulse (t,xx,xx1,u_cl,x_ref,N,rob_diam, conesx,conesy,cones_diam,cone_horizon)
% Displays the results of the simulation 


line_width = 1.5;
fontsize_labels = 14;

x_r_1 = [];
y_r_1 = [];


figure(1)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);


for k = 1:size(xx,2)
  
    
    % plot exhibited trajectory
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1); vx= xx(4,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on 
    

    
    hold off
    ylabel('$y$-position (m)','interpreter','latex')
    xlabel('$x$-position (m)','interpreter','latex')
    title('Impulse Input Analysis of MPC Controller', 'interpreter', 'latex' );
    axis([0 27 -0.5 1])
    box on;
    grid on 
    drawnow
    
end
hold on
stairs([0,2,2,30],[0,100,0,0], 'b--')



% Plot the Control Actions and Velocity
figure(2)
subplot(3,1,1)
stairs(t,u_cl(:,1),'k','linewidth',1.5); axis([0 t(end) -1000 1000])
xlabel('time (seconds)', 'interpreter', 'latex')
ylabel('Frx (N)','interpreter', 'latex')
grid on
subplot(3,1,2)
stairs(t,u_cl(:,2),'r','linewidth',1.5); axis([0 t(end) -0.85 0.85])
xlabel('time (seconds)','interpreter', 'latex')
ylabel('$\sigma$ (rad)','interpreter', 'latex')
subplot(3,1,3)
stairs(t,xx(4,1:end-1,1),'b','linewidth',1.5); axis([0 t(end) 0 5])
xlabel('time (seconds)','interpreter', 'latex')
ylabel('Velocity (m/s)','interpreter', 'latex')
grid on