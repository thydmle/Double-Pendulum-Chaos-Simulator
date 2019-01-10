q% Thy Doan Mai Le
% Double Pendulum with l1 = l2 ; 4/5/2018
% Please pray that my code works

%initialize initial angles in radians
init = [pi/8 0 pi/4 0];
g = 9.8;
%interval of integration
l_1 = 1;               %length of first pendulum (m)
l_2 = 1;               %length of second pendulum (m)

tfinal = 50000;
period1 = 2*pi*sqrt(l_1/g);
period2 = 2*pi*sqrt(l_2/g);

[t, y] = ode45(@DoublePendulum, [0 tfinal], init);

theta1_pos = y(:,1);
theta1_vel = y(:,2);
theta2_pos = y(:,3);
theta2_vel = y(:,4);

num_periods1 = floor(tfinal./(period1));
interp_t = [1:num_periods1].*period1;

num_periods2 = floor(tfinal./(period2));
interp_t2 = [1:num_periods2].*period2;

theta_1interp = interp1(t, theta1_pos, interp_t);
theta_1dotinterp = interp1(t, theta1_vel, interp_t);
theta_2interp = interp1(t, theta2_pos, interp_t2);
theta_2dotinterp = interp1(t, theta2_vel, interp_t2);

figure
plot(theta1_pos, theta1_vel, 'b.',theta_1interp, theta_1dotinterp, 'r.' )
xlabel('Angular Position');
ylabel('Angular Velocity');
title('Poincare map of mass 1');

figure
plot(theta2_pos, theta2_vel, 'b.', theta_2interp, theta_2dotinterp, 'r.')
xlabel('Angular Position');
ylabel('Angular Velocity');
title('Poincare map of mass 2');

function yprime = DoublePendulum(t, y)
    g = 9.8;
    l_1 = 1;
    l_2 = 1;
    C = cos(y(1) - y(3));
    S = sin(y(1) - y(3));
    Q = 1;
    omega_1 = sqrt(g/l_1);
    omega_2 = sqrt(g/l_2);
    yprime = [y(4); -2.*omega_1.^2.*sin(y(3)) + 2.*Q.*S.*(y(2).^2) + 2.*C.*omega_2.^2.*sin(y(1)) + C.*S.*(y(4).^2); ...
        y(2); ((-2./Q).*omega_2.^2.*sin(y(1)) - S./Q.*(y(4).^2) + C.*omega_1.^2.*sin(y(3)) - C.*S.*(y(2).^2))];
        
end


