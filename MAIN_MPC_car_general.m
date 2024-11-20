close all
clear all
clc
warning off;

%% Create an object for the support functions.
constants=initial_constants();

%% Load the constant values needed in the main file
Ts=constants('Ts');
outputs=constants('outputs');
hz = constants('hz');
inputs=constants('inputs');
trajectory=constants('trajectory');

%% Create the time array
t = 0:Ts:constants('time_length');

%% Import trajectory generation values
[x_dot_ref,y_dot_ref,psi_ref,X_ref,Y_ref]=trajectory_generator(t);
sim_length=length(t);

%% Generate the reference signal array
refSignals=zeros(1,length(X_ref)*outputs);
k=1;
for i = 1:outputs:length(refSignals)
   refSignals(i)=x_dot_ref(k);
   refSignals(i+1)=psi_ref(k);
   refSignals(i+2)=X_ref(k);
   refSignals(i+3)=Y_ref(k);
   k=k+1;
end
clear i k

%% Load the initial states
x_dot=x_dot_ref(1);
y_dot=y_dot_ref(1);
psi=psi_ref(1);
psi_dot=0;
X=X_ref(1);
Y=Y_ref(1);

%% Create state arrays
states=[x_dot,y_dot,psi,psi_dot,X,Y+30];
statesTotal=zeros(length(t),length(states));
statesTotal(1,:)=states;

%% Accelerations
x_dot_dot=0;
y_dot_dot=0;
psi_dot_dot=0;

accelerations=[x_dot_dot,y_dot_dot,psi_dot_dot];
accelerations_total=zeros(length(t),length(accelerations));

%% Initiate the controller - simulation loops
U1=0;
U2=0;
UTotal=zeros(length(t),2);
UTotal(1,1)=U1;
UTotal(1,2)=U2;

du=zeros(inputs*hz,1);


du=zeros(inputs*hz,1);
times = zeros(sim_length-1, 1);
error_x = 0;
error_y = 0;
error_psi = 0;
error_x_dot = 0;
%% Start with the loop
k=1; % for reading reference signals
for i =1:sim_length-1
    
    %% Generate discrete LPV Ad, Bd, Cd, Dd matrices
    [Ad, Bd, Cd, Dd]=state_space(states,U1,U2);
    
    %% Generating the current state and the reference vector
    x_aug_t=[states';U1;U2];
    
    k=k+outputs;
    if k+outputs*hz-1 <= length(refSignals)
        r=refSignals(k:k+outputs*hz-1);
    else
        r=refSignals(k:length(refSignals));
        hz=hz-1;
    end
    
    
    %% Generate simplification matrices for the cost function
    [Hdb,Fdbt,Cdb,Adc,G,ht] = mpc_simplification(Ad,Bd,Cd,Dd,hz,x_aug_t,du);
    ft=[x_aug_t',r]*Fdbt;
    
    
    tic;
    %% Calling the optimizer (quadprog)
    
    % Cost function in quadprog: min(du)*1/2*du'Hdb*du+f'du
    % Hdb must be positive definite for the problem to have finite minimum.
    options = optimoptions('quadprog','Display', 'off','LinearSolver','dense');
    
    [du,fval]=quadprog(Hdb,ft,G,ht,[],[],[],[],du,options);
    %[du,fval]=quadprog(Hdb,ft,[],[],[],[],[],[],[],options);
    times(i)=toc;
    
    if length(du)==0
        'The solver could not find the solution'
    end
    
    U1=U1+du(1);
    U2=U2+du(2);
    
    UTotal(i+1,1)=U1;
    UTotal(i+1,2)=U2;
    
    % Simulate the new states
    time_interval=(Ts)/30;
    T = (Ts)*(i-1):time_interval:Ts*(i-1)+(Ts);
    [T,x]=ode45(@(t,x) open_loop_new_states(t,x,[U1,U2]),T,states);
    
    states=x(end,:);
    statesTotal(i+1,:)=states;
    if i >sim_length/2
        error_x = error_x +abs(states(5)-refSignals(4*i+3));
        error_y = error_y +abs(states(6)-refSignals(4*i+4));
        error_psi = error_psi +abs(states(3)-refSignals(4*i+2));
        error_x_dot = error_x_dot +abs(states(1)-refSignals(4*i+1));
    end
    
    % Accelerations
    x_dot_dot=(x(end,1)-x(end-1,1))/time_interval;
    y_dot_dot=(x(end,2)-x(end-1,2))/time_interval;
    psi_dot_dot=(x(end,4)-x(end-1,4))/time_interval;
    
    accelerations=[x_dot_dot,y_dot_dot,psi_dot_dot];
    accelerations_total(i+1,:)=accelerations;
    
    if mod(i,500)==0
        'Progress (%) '
        i/sim_length*100
    end
    
end



avg_time = mean(times);
fprintf('Thời gian chạy trung bình: %.6f giây\n', avg_time);
fprintf('Thời gian chạy lâu nhất trong một vòng lặp: %.6f giây\n', max(times));
fprintf('error_x: %.6f \n', error_x*2/sim_length);
fprintf('error_y: %.6f \n', error_y*2/sim_length);
fprintf('error_psi: %.6f \n', error_psi*2/sim_length);
fprintf('error_x_dot: %.6f \n', error_x_dot*2/sim_length);
%% Plot the results

% Plot the trajectory
figure;
plot(X_ref,Y_ref,'--b','LineWidth',2)
hold on
plot(statesTotal(:,5),statesTotal(:,6),'r','LineWidth',1)
grid on;
xlabel('x-position [m]','FontSize',15)
ylabel('y-position [m]','FontSize',15)
%plotObstacles(Obstacles);
legend({'position-ref','position'},'Location','southeast','FontSize',15)

% Plot the inputs
figure;
subplot(2,1,1)
plot(t,UTotal(:,1),'b','LineWidth',2)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('steering wheel angle (delta) [rad]','FontSize',15)
legend({'delta'},'Location','southeast','FontSize',15)

hold on
subplot(2,1,2)
plot(t,UTotal(:,2),'b','LineWidth',2)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('applied acceleration (a) [m/s^2]','FontSize',15)
legend({'a'},'Location','southeast','FontSize',15)

% Plot Psi, X, Y
figure;
subplot(3,1,1)
plot(t,psi_ref,'--b','LineWidth',2)
hold on
plot(t,statesTotal(:,3),'r','LineWidth',1)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('psi-position [rad]','FontSize',15)
legend({'Psi-ref','psi'},'Location','southeast','FontSize',15)

hold on
subplot(3,1,2)
plot(t,X_ref,'--b','LineWidth',2)
hold on
plot(t,statesTotal(:,5),'r','LineWidth',1)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('X-position [m]','FontSize',15)
legend({'X-ref','X'},'Location','southeast','FontSize',15)

hold on
subplot(3,1,3)
plot(t,Y_ref,'--b','LineWidth',2)
hold on
plot(t,statesTotal(:,6),'r','LineWidth',1)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('Y-position [m]','FontSize',15)
legend({'Y-ref','Y'},'Location','southeast','FontSize',15)

% Plot the velocities
figure;
subplot(3,1,1)
plot(t,x_dot_ref,'--b','LineWidth',2)
hold on
plot(t,statesTotal(:,1),'r','LineWidth',1)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('dx [m/s]','FontSize',15)
legend({'dx-ref','dxt'},'Location','southeast','FontSize',15)

hold on
subplot(3,1,2)
plot(t,y_dot_ref,'--b','LineWidth',2)
hold on
plot(t,statesTotal(:,2),'r','LineWidth',1)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('dy [m/s]','FontSize',15)
legend({'dy-ref','dy'},'Location','southeast','FontSize',15)

hold on
subplot(3,1,3)
plot(t,statesTotal(:,4),'r','LineWidth',1)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('dpsi [rad/s]','FontSize',15)
legend({'dpsi'},'Location','southeast','FontSize',15)

% Plot the accelerations
figure;
subplot(3,1,1)
plot(t,accelerations_total(:,1),'b','LineWidth',2)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('ddx [m/s^2]','FontSize',15)
legend({'ddx'},'Location','southeast','FontSize',15)

hold on
subplot(3,1,2)
plot(t,accelerations_total(:,2),'b','LineWidth',2)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('ddy [m/s^2]','FontSize',15)
legend({'ddy'},'Location','southeast','FontSize',15)

hold on
subplot(3,1,3)
plot(t,accelerations_total(:,3),'b','LineWidth',2)
grid on
xlabel('t-time [s]','FontSize',15)
ylabel('ddpsi [rad/s^2]','FontSize',15)
legend({'ddpsi'},'Location','southeast','FontSize',15)

