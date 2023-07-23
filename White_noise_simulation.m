%% Simulation for White noise
%% Problem1.
dt = [0.01, 0.1];
mu = 0;
sigma = 1;

rng('shuffle');

data = struct;
n_sim = 1000;

for i = 1:length(dt)
    N = 10.0/dt(i); % time
    data(i).t = zeros(N+1,1);
    data(i).P_k = zeros(N+1,1); % initial condition : P(0)=0
    data(i).x_k = zeros(N+1,n_sim+1);
    data(i).mean = zeros(N+1,1); % initial condition : m(0)=0
    data(i).w = zeros(n_sim+1,1);
    
    for j=2:N+1
        data(i).t(j) = dt(i)*(j-1);
        for m = 1:n_sim
            w = sigma * randn(1) + mu;
            if j==2
                data(i).w(m) = w;
            end
            data(i).x_k(j, m) = data(i).x_k(j-1,m) + w*dt(i); % add the w value to previous time step value
            data(i).P_k(j) = data(i).P_k(j) + data(i).x_k(j,m)^2;
        end
        data(i).P_k(j) = 1/n_sim * data(i).P_k(j); % mean value of variance of simulation
        data(i).mean(j) = 1/n_sim * sum(data(i).x_k(j,:)); % mean value of x at t = j*dt
    end
end

% comparison between simulation and theory
P_theory = 0:dt(1):10;

% mean/variance plot
figure;
subplot(2,1,1);
plot(data(1).t, data(1).mean); grid; title('Pr1. mean of x_k after simulation = '+string(n_sim)+' dt='+string(dt(1)));
subplot(2,1,2);
plot(data(1).t, data(1).P_k); grid;
title('Pr1. variance of x $(P_k)$ after simulation = '+string(n_sim)+' dt='+string(dt(1)));
xlabel('t');
ylabel('P_{k}');


figure;
subplot(2,1,1);
plot(data(2).t, data(2).mean); grid; title('Pr1. mean of x_k after simulation = '+string(n_sim)+' dt='+string(dt(2)));
subplot(2,1,2);
plot(data(2).t, data(2).P_k); grid;
title('Pr1. variance of x (P_k) after simulation = '+string(n_sim)+' dt='+string(dt(2)));
xlabel('t');
ylabel('P_{k}');

% check white noise
figure;
h = histogram(data(1).w,'Normalization','probability');
title('Distribution of white noise term (simulation'+string(n_sim)+')');

%% Problem2.

for i = 3:length(dt)+2
    % initialization
    N = 10.0/dt(i-2); % time
    data(i).t = zeros(N+1,1);
    data(i).P_k = zeros(N+1,1); % initial condition : P(0)=0
    data(i).P_k(1) = 1;
    data(i).x_k = zeros(N+1,n_sim+1);
    for m = 1:n_sim
        data(i).x_k(1,:) = randn(1,n_sim+1); % mean 0, variance 1 at initial time
    end
    data(i).mean = zeros(N+1,1); % initial condition : m(0)=0
    
    for j=2:N+1
        data(i).t(j) = dt(i-2)*(j-1);
        for m = 1:n_sim
            w = sigma * randn(1) + mu;
            data(i).x_k(j, m) = data(i).x_k(j-1,m) + w*dt(i-2); % add the w value to previous time step value
            data(i).P_k(j) = data(i).P_k(j) + data(i).x_k(j,m)^2;
        end
        data(i).P_k(j) = 1/n_sim * data(i).P_k(j); % mean value of variance of simulation
        data(i).mean(j) = 1/n_sim * sum(data(i).x_k(j,:)); % mean value of x at t = j*dt
    end
end

% mean/variance plot
figure;
subplot(2,1,1);
plot(data(3).t, data(3).mean); grid; title('Pr2. mean of $x_k$ after simulation = '+string(n_sim)+' dt='+string(dt(1)));
subplot(2,1,2);
plot(data(3).t, data(3).P_k); grid;
title('Pr2. variance of x $(P_k)$ after simulation = '+string(n_sim)+' dt='+string(dt(1)));
xlabel('t');
ylabel('P_{k}');

figure;
subplot(2,1,1);
plot(data(4).t, data(4).mean); grid; title('Pr2. mean of x_k after simulation = '+string(n_sim)+' dt='+string(dt(2)));
subplot(2,1,2);
plot(data(4).t, data(4).P_k); grid;
title('Pr2. variance of x $(P_k)$ after simulation = '+string(n_sim)+' dt='+string(dt(2)));
xlabel('t');
ylabel('P_{k}');