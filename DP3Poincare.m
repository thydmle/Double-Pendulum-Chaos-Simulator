% Thy Doan Mai Le
% Double Pendulum with l1 > l2 ; 4/5/2018
% Please pray that my code works

%initialize initial angles in radians
init = [pi/8 0 pi/4 0];

%interval of integration
l_1 = 4;               %length of first pendulum (m)
l_2 = 1;               %length of second pendulum (m)
g = 9.8;

tfinal = 60000;
period1 = 2*pi*sqrt(l_1/g);
period2 = 2*pi*sqrt(l_2/g);

num_periods1 = floor(tfinal./period1);
interp_t1 = [1:num_periods1].*period1;

num_periods2 = floor(tfinal./period2);
interp_t2 = [1:num_periods2].*period2;

[t, y] = ode45(@DoublePendulum, [0 tfinal], init);

% taking out the stuff that I need
theta1 = y(:,1);        % the position of theta_1
theta2 = y(:,3);        % position of theta_2
theta1_dot = y(:,2);     % velocity of theta_1
theta2_dot = y(:,4);   % velocity of theta_2

theta_1interp = interp1(t, theta1, interp_t);
theta_1dotinterp = interp1(t, theta1_dot, interp_t);
theta_2interp = interp1(t, theta2, interp_t2);
theta_2dotinterp = interp1(t, theta2_dot, interp_t2);

figure(2)
plot( theta1, theta1_dot, 'b.',theta_1interp, theta_1dotinterp, 'r.'); 
xlabel('$$\theta_1$$', 'interpreter', 'latex');
ylabel('$$\dot{{\theta_1}}$$', 'interpreter', 'latex');
title('Poincare Map for Mass 1');

figure(3)
plot(theta2, theta2_dot, 'b.',theta_2interp, theta_2dotinterp, 'r.');
xlabel('$$\theta_2$$', 'interpreter', 'latex');
ylabel('$$\dot{\theta_2}$$', 'interpreter', 'latex');
title('Poincare Map for Mass 2');



function yprime = DoublePendulum(t, y)
    g = 9.8;
    l_1 = 4;
    l_2 = 1;
    C = cos(y(1) - y(3));
    S = sin(y(1)- y(3));
    Q = 1;
    omega_1 = sqrt(g/l_1);
    omega_2 = sqrt(g/l_2);
    yprime =  [y(4); -2.*omega_1.^2.*sin(y(3)) + 2.*Q.*S.*(y(2).^2) + 2.*C.*omega_2.^2.*sin(y(1)) + C.*S.*(y(4).^2); ...
        y(2); ((-2./Q).*omega_2.^2.*sin(y(1)) - S./Q.*(y(4).^2) + C.*omega_1.^2.*sin(y(3)) - C.*S.*(y(2).^2))];
    
end

