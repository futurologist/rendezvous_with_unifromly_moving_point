function main()
    h = 0.3;
    a0 = 0.1;
    u_max = .8; % velocity threshold
    xB0 = [0; 0; 0];
    vB0 = [-1; 2; 0];
    xA0 = [ 7; 12; 0] + xB0;
    vA0 = [1; 5; 0]/7;
    %vA0 =  [2; -1; 0];
    L_unit = cross(xA0 - xB0, vA0 - vB0);
    L_unit =  L_unit / norm(L_unit);
    t = 0;
    xB = xB0;
    x = xA0;
    v = vA0;
    hold on
    grid on
    %axis([-200 20 -100 350])
    plot(0, 0, 'bo')

    % STEP 0 (the motion as described in the text above):
    n = floor(sqrt(2)*norm(vA0 - vB0)/(a0*h));
    for i=1:n
        [t, x, v, xB] = E(t, x, v, xB, vB0, a0, L_unit, h);
        plot(x(1), x(2), 'ro');
        plot(xB(1), xB(2), 'bo');
        pause(0.1)
    end
    u = v - vB0;
    norm_u = norm(u);
    % short additional decceleration so that A attains velocity v = vB0 
    t0 = t + norm_u/a0;    
    n = floor((t0 - t)/h); 
    a = - a0 * u / norm_u;    
    for i=1:n
        [t, x, v, xB] = ET(t, x, v, xB, vB0, a, h);
        plot(x(1), x(2), 'ro');
        plot(xB(1), xB(2), 'bo');
        pause(0.1)
    end
    [t, x, v, xB] = ET(t, x, v, xB, vB0, a, t0-t);
    plot(x(1), x(2), 'ro');
    plot(xB(1), xB(2), 'bo');
    pause(0.1)
    
    % STEP 1 (uniform acceleration of magnitude a0):
    v = vB0; 
    a = x-xB;
    norm_y0 = norm(a);
    a = - a0 * a / norm_y0;
    %t2 = t1 + sqrt( norm_y/a0 ); 
    accel_time = min( u_max/a0, sqrt( norm_y0/a0 ) );
    t1 = t0 + accel_time;
    n = floor((t1 - t0)/h);     
    for i=1:n
       [t, x, v, xB] = ET(t, x, v, xB, vB0, a, h);
       plot(x(1), x(2), 'bo');
       plot(xB(1), xB(2), 'ro');
       pause(0.1)
    end 
    [t, x, v, xB] = ET(t, x, v, xB, vB0, a, t1-t);
    plot(x(1), x(2), 'bo');
    plot(xB(1), xB(2), 'ro');
    pause(0.1)
    
    % STEP 2 (uniform straight-line motion): 
    norm_y1 = norm(x-xB);
    norm_y12 = max(0, norm_y0 - 2*(norm_y0 - norm_y1));
    t12 = norm_y12 / norm(v-vB0)
    t = t + t12
    n12 = floor(t12/h)
    for i=1:n12
        x = x + h*v; 
        xB = xB + h*vB0;
        plot(x(1), x(2), 'ro');
        plot(xB(1), xB(2), 'bo');
        pause(0.1)
    end
    x = x + (t12-n12*h)*v; 
    xB = xB + (t12-n12*h)*vB0;
    plot(x(1), x(2), 'ro');
    plot(xB(1), xB(2), 'bo');
    pause(0.1)

    % STEP 3 (uniform deceleration of magnitude a0, symmetric to STEP 1): 
    a = -a;
    for i=1:n % t2 + (t2-t1)
       [t, x, v, xB] = ET(t, x, v, xB, vB0, a, h);
       plot(x(1), x(2), 'bo');
       plot(xB(1), xB(2), 'ro');
       pause(0.1)
    end
    [t, x, v, xB] = ET(t, x, v, xB, vB0, a, t0+t12+2*accel_time-t);
    plot(x(1), x(2), 'bo');
    plot(xB(1), xB(2), 'ro');
    pause(0.1)
    norm(x-xB)
    norm(v-vB0)
end


% Here are the additional functions that are used in the main code above:

% change of coordinates from world coordinates x, v 
% to coordinates y, u from spaship B's point of view:
function [y, u] = change(x, v, xB, vB0)
    y = x - xB;
    u = v - vB0;
end

% inverse chage of coordinates from y, u to x, v
function [x, v] = inv_change(y, u, xB, vB0)
    x = y + xB;
    v = u + vB0;
end

% solution to the second system of differential equations for a step h:
function [y_out, u_out] = R(y, u, a0, L_unit, h)
   omega = a0 / (sqrt(2) * norm(u));
   L_x_u = cross(L_unit, u);
   cos_omega_h = cos(omega*h);
   sin_omega_h = sin(omega*h);
   omega = 2*omega;
   y_out = y + (L_x_u ...
       + sin_omega_h * u  -  cos_omega_h * L_x_u) / omega;      
   u_out = cos_omega_h * u  +  sin_omega_h * L_x_u;
end

% solution to the first system of differential equations for a step h:
function [y_out, u_out] = T(y, u, a0, h)
    sqrt_2 = sqrt(2);
    u_unit = u / norm(u);  
    y_out = y  +  h * u/2  -  h^2 * a0 * u_unit/ (4*sqrt_2);
    u_out = u - h * a0 * u_unit / sqrt_2;
end

% approximate solution of the original system of differential equations for step h
% i.e. the sum of furst and second systems of differential equations:
function [t_out, x_out, v_out, xB_out] = E(t, x, v, xB, vB0, a0, L_unit, h)
   t_out = t + h;
   [y, u] = change(x, v, xB, vB0);
   [y, u] = R(y, u, a0, L_unit, h/2);
   [y, u] = T(y, u, a0, h);
   [y, u] = R(y, u, a0, L_unit, h/2);
   xB_out = xB + h*vB0;
   [x_out, v_out] = inv_change(y, u, xB_out, vB0);  
end

% straight-line motion with constant acceleration: 
function [t_out, x_out, v_out, xB_out] = ET(t, x, v, xB, vB0, a, h)
    t_out = t + h;
    [y, u] = change(x, v, xB, vB0);
    y = y  +  h * u  +  h^2 * a / 2;
    u = u + h * a;
    xB_out = xB + h*vB0;
    [x_out, v_out] = inv_change(y, u, xB_out, vB0);
end



