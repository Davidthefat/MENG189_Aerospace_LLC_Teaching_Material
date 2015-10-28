function DropSim(n)
%DROPSIM Simulates a trajectory of an unpowered rocket projectile.
%   n: maximum number of discrete steps in simulation

%Initialize all variables used in calculations
x0 = [0, 1750];  %initial positions (x,y) (m)
v0 = [60, 60];  %initial velocities (x,y) (m/s)
g = 9.81;       %gravitational acceleration (m/s^2)
dt = 0.005;     %time step (sec)
m = 8;          %mass of projectile (kg)
length = 20;    %length of projectile (m)
cg = 2;         %distance of center of gravity from tip (m)
cp = 9;         %distance of center of pressure from tip (m)
cD = 0.8;       %coefficient of drag of cylindrical body
cDF = 1.28;     %coefficient of drag of fin
ref_A = 0.025;  %reference area of drag (m^2)
fin_A = 0.0775; %area of fins for lift calculations (m^2)
rho = 1.225;    %density of air (kg/m^3)

%Initialize 2 x n matrices of zeros to hold values of position, velocity and
%   acceleration, respectively
x = zeros(2,n);
v = zeros(2,n);
a = zeros(2,n);

%Initialize n sized arrays of zeroes to hold values of angle of projectile,
%   angular velocity and dirction of velocity of vehicle
B = zeros(1,n);
w = zeros(1,n);
T = zeros(1,n);

%Set the first columns of position and velocity to initial values 
x(1,1) = x0(1);
x(2,1) = x0(2);
v(1,1) = v0(1);
v(2,1) = v0(2);

%Set the first columns of angles to reflect initial values 
B(1) = -atan2(v(2,1),v(1,1));
T(1) = B(1); 

%Set the end step number to maximum
e = n;

%Begin iterations for calculations
for i = 2:n
    %Calculate direction of flow
    T(i) = -atan2(v(2,i-1),v(1,i-1));
    %Calculate net velocity of projectile
    v_net = sqrt(v(1,i-1)^2 + v(2,i-1)^2);
    %Calculate the angle of attack of projectile
    A = T(i)-B(i-1);
    %Calculate coefficient of lift based on flat plate approximation
    cL = 2*pi*sin(A);
    %Calculate lift of fins
    L = cL * rho * v_net^2 * fin_A / 2;
    %Calculate drag of fins as the body rotates
    D_fin = cDF * rho * ((cp-cg)*w(i-1))^2 * fin_A / 2;
    %Calculate angular velocity
    w(i) = w(i-1)+(L-D_fin)*(cp-cg)*dt/((length-cg)/length*m);
    %Calculate angle of projectile respective to environment
    B(i) = B(i-1)+w(i)*dt;
    %Calculate drag of projectile assuming flow is perfectly perpendicular
    %   to projectile
    D = cD * rho * v_net^2 * ref_A / 2;
    %Horizontal component of drag is only force in the horizontal direction
    a(1,i) = -D*cos(T(i))/m;
    %Vertical component of acceleration includes gravity and drag
    a(2,i) = D*sin(T(i))/m-g;
    %Velocity and acceleration are calculated using Euler method
    v(1,i) = v(1,i-1)+dt*a(1,i);
    v(2,i) = v(2,i-1)+dt*a(2,i);
    x(1,i) = x(1,i-1)+dt*v(1,i);
    x(2,i) = x(2,i-1)+dt*v(2,i);
    %Once projectile reaches ground, exit calculation loop
    if x(2,i) < 0
        %e is the step number at the end of calculations
        e = i;
        break;
    end
end
%Plot the important values
figure;
plot(x(1,1:e),x(2,1:e));
title('Flight Path');
xlabel('X (m)');
ylabel('Y (m)');
figure;
plot((1:e)*dt, B(1,1:e).*180/pi);
title('Angle of Projectile VS Time');
ylabel('Degrees');
xlabel('Sec');
figure;
plot((1:e)*dt, (T(1,1:e)-B(1,1:e)).*180/pi);
title('Angle of Attack VS Time');
ylabel('Degrees');
xlabel('Sec');
figure;
plot((1:e)*dt, v(2,1:e));
title('Vertical Velocity of Projectile VS Time');
ylabel('m/s');
xlabel('Sec');

%Print out results
fprintf('Time to Impact: %3.4f Sec.\n',e*dt);
fprintf('Distance Downrange: %3.4f m.\n',x(1,e));
end

