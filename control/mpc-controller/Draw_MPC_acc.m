%   ______________________________________________________________________
%   **********************************************************************
%   Function uses code proposed in:
%   MMehrez,MPC-and-MHE-implementation-in-MATLAB-using-Casadi,(2021),GitHub
%   https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi
%   ______________________________________________________________________
%   **********************************************************************

function Draw_MPC_acc(t,xx,xx1,u_cl,x_ref,N,rob_diam,conesx,conesy,cones_diam,cone_horizon)
% Displays the results of the simulation for the acceleration based events

line_width = 1.5;
fontsize_labels = 14;

x_r_1 = [];
y_r_1 = [];

% car circle
r = rob_diam/2; 
ang=0:0.005:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

% cone circle
r = cones_diam/2;  
xp_obs=r*cos(ang);
yp_obs=r*sin(ang);

figure(1)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

a=cone_horizon/2;
cone_sample=1;
xs = x_ref(:,a);
cones_x = conesx(1:cone_horizon);
cones_y = conesy(1:cone_horizon);


for k = 1:size(xx,2)
  
    h_t = 0.14; w_t=0.09; %triangle parameters

    % finds line which car needs to cross to update cones and reference
    grad = (cones_y(3)-cones_y(4))/(cones_x(3)-cones_x(4));
    if (abs(grad) == inf)
        c = cones_x(3);
    elseif (grad == 0)
        c = cones_y(3); 
    else
        c = cones_y(3)-grad*cones_x(3); 
    end
    
    % updates cones and reference if line is crossed
    if grad == 0
        if (cones_x(4)>cones_x(3) && xx(2,k,1)<c)
            a=a+1;
            cone_sample = cone_sample+2;
        elseif (cones_x(4)<cones_x(3) && xx(2,k,1)>c)
            a=a+1;
            cone_sample = cone_sample+2;
        end    
    elseif (grad == inf && xx(1,k,1)>c)
        a=a+1;
        cone_sample = cone_sample+2;
    elseif (grad == -inf && xx(1,k,1)<c)
        a=a+1;
        cone_sample = cone_sample+2;     
    elseif (abs(grad) ~= inf && grad ~= 0)
        if(cones_y(3)>cones_y(4))
            if (grad>0 && xx(2,k,1)<grad*xx(1,k,1)+c)
                a=a+1;
                cone_sample = cone_sample+2;
            elseif (grad<0 && xx(2,k,1)>grad*xx(1,k,1)+c)
                a=a+1;
                cone_sample = cone_sample+2;
            end
        else
            if (grad>0 && xx(2,k,1)>grad*xx(1,k,1)+c)
                a=a+1;
                cone_sample = cone_sample+2;
            elseif (grad<0 && xx(2,k,1)<grad*xx(1,k,1)+c)
                a=a+1;
                cone_sample = cone_sample+2;
            end
        end
    end
    
    % Set reference position
    xs = x_ref(:,a);
    
    % Set cones positions
    cones_x = conesx(cone_sample:cone_horizon+cone_sample-1);
    cones_y = conesy(cone_sample:cone_horizon+cone_sample-1); 
    
    % plot reference state
    x1 = xs(1); y1 = xs(2); th1 = xs(3); 
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];
    fill(x1_tri, y1_tri, 'g'); 
    hold on;
    
    % plot exhibited trajectory
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1); vx= xx(4,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];
    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on 
    
    % plot prediction
    if k < size(xx,2) 
        plot(xx1(1:N,1,k),xx1(1:N,2,k),'r--*')
        for j = 2:N+1
            plot(xx1(j,1,k)+xp,xx1(j,2,k)+yp,'--r'); 
        end
    end
    
    fill(x1_tri, y1_tri, 'r'); % plot car current position  
    plot(x1+xp,y1+yp,'--r'); % plot car circle
    
    % plot cone circles    
    for n = 1:length(cones_x)
        plot(cones_x(n)+xp_obs,cones_y(n)+yp_obs,'--b'); 
    end
    
    hold off
    ylabel('$y$-position (m)','interpreter','latex')
    xlabel('$x$-position (m)','interpreter','latex')
    title(['Velocity of Car: ', num2str(vx)], 'interpreter', 'latex' );
    axis([0 100 -50 50])
    box on;
    grid on 
    drawnow
    
end


% Plot the 
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
stairs(t,xx(4,1:end-1,1),'b','linewidth',1.5); axis([0 t(end) 0 12])
xlabel('time (seconds)','interpreter', 'latex')
ylabel('Velocity (m/s)','interpreter', 'latex')
grid on