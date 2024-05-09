%%%%%%%%%%%%%%%%%%%%%% REPLICATION ASSIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Defining the parameters
beta = 0.9;
delta = 0.1;
lambda = 0.75;
mu = 1;
rho = 1-delta;
F = 0.2;
N = 20;
R = 6;                % Initial Cutoff level for capital replacement
P=[0.9,0.1;0.1,0.9];  % Transition matrix 

A_L = 0.75;
A_H = 1.25;
A = [A_H,A_L];

%% Creating the Grid for idiosyncratic productivity (epsilon)

epsilon_min = 0.4;
epsilon_max = 1.6;
epsilon = transpose(linspace(epsilon_min, epsilon_max, N));

%% Defining the Capital Stock

ones_vector = ones(1, R - 1); 
multiplied_vector = rho * ones_vector; 
cumulative_product = cumprod(multiplied_vector);
k = [1, cumulative_product]; 


%%  Initialzing objects

%%%%%%%%%%%%%%%  Creating a 3D revenue function on n-gridpoints %%%%%%%%%%          

b = epsilon*k;
r = b*A(1);
r(:,:,2)= b*A(2);

%%  3D guess for the Value Functions %%%%%%%%%%%%%%%%%%%%%

V_R = ones(size(r));   % Guess for Replacement 
V_N = zeros(size(r));  % Guess for no Replacement 
V0 = max(V_R,V_N);     % V0 guess 
 
tol = 10^(-12);        % Tolerance level

%%  Value Function Iteration Loop %%%%%%%%%%%%%%%%%%%%


dif = 1;
count = 0;
maxits = 1e9;

tic


while dif > tol && count < maxits
   
    z0 = zeros(size(r));
    
    
    for i = 1:2
        V_R(:, :, i) = r(:, :, i) * lambda - F + beta * (P(i, 1) * repmat(V0(:, 1, 1), 1, R) + P(i, 2) * repmat(V0(:, 1, 2), 1, R));
        for j = 1:R-1
            V_N(:, j, i) = r(:, j, i) + beta * (P(i, 1) * V0(:, j+1, 1) + P(i, 2) * V0(:, j+1, 2));
        end
    end
    
    
    V1 = max(V_R, V_N);
    z0(V1 == V_R) = 1;
    
    
      
    dif = max(abs(V1 - V0), [], 'all');

    
    count = count + 1;
    
    
    V0 = V1;
end

toc



%%  New Cutoff level for Capital replacement is 7, %%%%%%%%%%%%%%

R_new = 7;

ones_vector = ones(1, R_new - 1); 
multiplied_vector = rho * ones_vector; 
cumulative_product = cumprod(multiplied_vector);
k_1 = [1, cumulative_product]; 


%%%%%%%%%%  Recreating a 3D revenue function on n-gridpoints  %%%%%%%%%%%%%          

b_1 = epsilon*k_1;
r_1 = b_1*A(1);
r_1(:,:,2)= b_1*A(2);

%%%%%%%%%%%%%%% 3D guess for the Value Functions %%%%%%%%%%%%%%%%%%%%%%%

V_R_1 = ones(size(r_1));        % Replacement V guess
V_N_1 = zeros(size(r_1));       % No Replacement V guess
V0_0 = max(V_R_1,V_N_1);        % V0 guess 
 
tol = 10^(-12);                % Tolerance level


dif=1;
count=0;
maxits=1e9;

tic


while dif > tol && count < maxits
    
    z01 = zeros(size(r_1));

    
    for i = 1:2
        V_R_1(:, :, i) = r_1(:, :, i) * lambda - F + beta * (P(i, 1) * repmat(V0_0(:, 1, 1), 1, R_new) + P(i, 2) * repmat(V0_0(:, 1, 2), 1, R_new));
        for j = 1:R_new - 1
            V_N_1(:, j, i) = r_1(:, j, i) + beta * (P(i, 1) * V0_0(:, j+1, 1) + P(i, 2) * V0_0(:, j+1, 2));
        end
    end

    
    V1_1 = max(V_R_1, V_N_1);
    z01(V1_1 == V_R_1) = 1;

   
    dif = max(abs(V1_1 - V0_0), [], 'all');

    
    count = count + 1;

    
    V0_0 = V1_1;
end
toc



%%%%%%%%%%%%%%%%%%%%%  Plotting the Spy Plot %%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,1,1)
spy(z0(:,:,1)')
xlabel('Firm Productivity')
ylabel('Time Since Last Replacement')
title('Replacement in Low Aggregate Productivity')

subplot(2,1,2)
spy(z0(:,:,2)')
xlabel('Firm Productivity')
ylabel('Time Since Last Replacement')
title('Replacement in High Aggregate Productvity')

%%%%%%%%% COMMENTS ON THE IMPORTANT FEATURES OF THE HAZARD FUNCTION %%%%%%%

%%% The simulation result for the policy function depicts that firms will 
%%% replace their old capital once the cutoff is reached. Again, in the
%%% event of  high idiosyncratic shocks, the firms have a higher
%%% probability of investing in a high aggregate productivity
%%% period than a low aggregate productivity period.



%%%%%%%%  Plotting Hazard Function for Capital Replacement %%%%%%%%%%%%%%%%

H = zeros(R,2);

for i=1:2
    for j=1:R
        H(j,i)=sum(z0(:,j,i))/N;
    end
end

time = 1:R ;

figure(2)
plot(time,H(:,1), 'b--',time,H(:,2),'m-')
title('THEORITICAL HAZARD FOR MACHINE REPLACEMENT')
xlabel('Time Since Last Replacement')
ylabel('Probability of Replacement')
legend('High State','Low State','Location','best')
xlim([1 R])
ylim([0 1.05])

%%%%%%%  Time Series for the evolution of one firm in the model %%%%%%%%%%%

ts=40;
T=160;
Time=(1:ts);
E=randi(N,1,ts);
esim=epsilon(E);

AT=zeros(1,T);
AT(1)=2;

for i=2:T
    if AT(i-1)==1
        AT(i)=randsample(1:2,1,true,P(1,:));
    else
        AT(i)=randsample(1:2,1,true,P(2,:));
    end
end

Asim=A(AT);

Y=zeros(1,ts);
K=zeros(1,ts+1);
K(1)=1;
sim_K=zeros(1,ts+1);
sim_K(1)=1;

for i=1:ts
    simoutput(i)=K(K(i))*Asim(i)*esim(i)*(1-z0(E(i),K(i),AT(i)))+z0(E(i),K(i),AT(i))*(lambda*K(K(i))*Asim(i)*esim(i)-F);
    if z0(E(i),K(i),AT(i))==1
        K(i+1)=1;
        sim_K(i+1)=K(K(i+1));
    else
        K(i+1)=K(i)+1;
        sim_K(i+1)=K(K(i+1));
    end
end

figure(3)
plot(Time,sim_K(1:ts),'b--',Time,simoutput,'m-')
title('Simulated Capital')
xlabel('Period')
ylabel('Evolution of the Firm')
ylim([0 max(simoutput)+0.1])
legend('Capital','Output','location','best')


%%%%%%  Replicating figure 3 with fixed Aggregate Productivity Shock %%%%%%


tt=50;
n=6;                     % Setting a specific number of firms
ir=zeros(tt,1);
w=zeros(N,tt+1);
w(:,1)=ones(N,1);        % 1 firm at each vintage of capital to start
A_index=1;               % Fixing Aggregate Productivity Shock

for i=1:tt
    for j=1:n
    ir(i)=ir(i)+H(j,A_index)*w(j,i);
    end
    w(1,i+1)=ir(i);
    for j=2:n
        w(j,i+1)=(1-H(j-1,A_index))*w(j-1,i);
    end
end

ir=ir/n;

figure(4)
plot(1:tt,ir,'magenta')
title('CONVERGENCE WITHOUT AGGREGATE SHOCKS - BASELINE PARAMETERS')
xlabel('Period')
ylabel('Investment Rate')

%%%%%%%%%  Replicating figure 4 when "A" follows a Markov process  %%%%%%%%

ir2=zeros(T-1,1);
w1=zeros(n,T);
w1(:,1)=ones(n,1); 

for i=1:T-1
    for j=1:n
    ir2(i)=ir2(i)+H(j, AT(i))*w1(j,i);
    end
    w1(1,i+1)=ir2(i);
    for j=2:n
        w1(j,i+1)=(1-H(j-1,AT(i)))*w1(j-1,i);
    end
end

ir2=ir2/n;

figure(5)
[y,line1,line2]=plotyy(1:T-1,ir2,1:T,Asim);
line1.Color = 'black'; 
line2.Color = 'magenta';
ylabel(y(1),'Investment Rate',Color="black")
ylabel(y(2),'Aggregate State',Color="black")
xlabel('period','FontWeight','bold')
line2.LineStyle=':';
line2.Marker='*';
title('AGGREGATE INVESTMENT FLUCTUATIONS - BASELINE SIMULATIONS')

%%  Implications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From figure 3 , we can see that over the 40 periods output is more
% volatile than capital. We can infer from figure 4, investement rate is
% volatile in the early periods of production but converges to a
% steady-state as time goes on when there no aggreagate shocks. This shows
% that firms engage in lumpy investment.








