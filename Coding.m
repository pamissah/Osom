%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coding Assignment (Part A) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Initial Variable and Parameter Setup and Preallocation

%%% Model Parameters
beta = 0.96;
gamma = 1.3;
r = 0.04;
sigma = 0.04;
rho = 0.9;             

%%% Income Grid Setup
Y_n = 5;                  
sd = Y_n/2-0.5;           
Y = linspace(-sd*sigma, sd*sigma, Y_n);      

%%% Asset Grid Setup
a_n = 1000;               
a_max = 4 * exp(max(Y));    
A = linspace(-exp(min(Y))/r, a_max, a_n)';

%% Calculations that do not need to be repeated
P = ones(Y_n) / Y_n; 
c_choice=zeros(a_n,Y_n,a_n);       
utility=c_choice;                  

for ap = 1:a_n
    c(:, :, ap)=(1 + r)*repmat(A, 1, Y_n)+exp(repmat(Y, a_n, 1))-repmat(A(ap), 1, Y_n);
end

c(c<0) = 0;
if gamma == 1
    utility = log(c);
else
    utility(c==0) = -inf;
    utility(c>0) = c(c>0).^(1-gamma)/(1-gamma);
end

%%% VFI Preallocations and Tolerances
tol = 10^(-9);                  
maxits = 10^4;            

V0 = zeros(a_n, Y_n);    
V1 = V0;                                      


%%  Main VFI Loop

count = 0;
dif = inf;
while dif>tol && count<maxits
    
    for y = 1:Y_n
        for ap = 1:a_n
            V_candidate = utility(:,y,ap)+beta*V0*P(y,:)';
            [V1(:,y), ~] = max(V_candidate, [], 2);
        end
    end
    dif = max(abs(V0(:)-V1(:)));
    V0 = V1;
    count = count + 1;
end

%%  Recovery of Consumption Policy Function

a_prime = ones(size(A));
c_policy=(1+r)*A +exp(Y(a_prime))-A(a_prime);

%%  Plots

figure(1)
plot(A, V1(:,1), A, V1(:,round(Y_n/2)), A, V1(:,Y_n))
xlabel('Assets'), ylabel('Value'), title('Value Function')
legend('Minimum Income', 'Medium Income', 'Maximum Income', 'Location', 'SouthOutside', 'Orientation', 'Horizontal')

figure(2)
plot(A, c_policy(:,1), A, c_policy(:,round(Y_n/2)), A, c_policy(:,Y_n))
xlabel('Assets')

figure(3)
plot(A, a_prime, A, a_prime, A, a_prime)
xlabel('Assets')
ylabel('Assets')
title('Optimal Savings')
legend('Minimum Income','Steady State Income','High Income','location','southoutside','orientation','horizontal')

%%  Simulations

sims = 1000;  
y_sim = zeros(sims, 1); 
c_sim = zeros(sims, 1);
a_sim = zeros(sims + 1, 1);

current_state = randi(Y_n, 1);

for t = 1:sims
    
    rnd = rand;
    next_state = find(cumsum(P(current_state, :)) >= rnd, 1, 'first');
    y_sim(t) = next_state;
    
    current_state = next_state;
end

figure(4)
subplot(3,1,1)
plot(sims/2+1:sims,exp(Y(y_sim(sims/2+1:sims))),sims/2+1:sims,ones(sims/2,1))
xlabel('Time')
ylabel('Income')

subplot(3,1,2)
plot(sims/2+1:sims,c_sim(sims/2+1:sims))
xlabel('Time')
ylabel('Consumption')

subplot(3,1,3)
plot(sims/2+1:sims,a_sim(sims/2+1:sims))
xlabel('Time')
ylabel('Assets')


%% Calculating a correlogram between the simulated income and consumption series to 4 lags with xcorr function .

%%% Transforming the simulated income and consumption into actual values
income = exp(Y(y_sim));
consumption = c_sim; 

%%% Calculating of the cross-correlation with 4 lags
lags = -4:4; 
[correlation, lag] = xcorr(income, consumption, 4); 

%%% Plotting the correlogram
figure (5);
stem(lag, correlation); 
xlabel('Lag'); ylabel('Cross-correlation'); title('Correlogram between Income and Consumption');










