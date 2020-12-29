%% Simulation Study: HATS VS SCM 

% 1. Set parameters to simulate the data, simulation and optimization settings
% 2. Plot true covariance, model function and example of realizations from this parametric setting
% 3. Inference: For different number of realizations (num_real_vec),
%    repeat num_sim simulations where the parameters and the covariance matrix are estimated with HATS and SCM
% 4-5-6. Compute distance from true covariance, confidence intervals for parameter estimates, and plot figures 

clear all
close all
addpath('functions')
    
%% 1. Set parameters to simulate the data
N = 256; 
T0 = 0; % initial time
Tf = 0.5; % final time
delta_t= abs(Tf-T0)/(N-1); % sampling interval
time_vec = T0 + [0:N-1]'* delta_t;

% Model parameters
a_q_true = 600;
b_q_true = 0.20;
c_q_true = 1000;
c_r_true = 70000;
noise_true = 100;

% Simulation settings
numReal_vec = [1 5 10 20 50]; % number of realizations considered
num_sim = 10; % number of simulations for every number of realizations in numReal_vec

% Note: The setting of the simulation study presented in the paper (which takes some time to run!) is the
% following: 
% numReal_vec = [1 5 10 20 30 40 50 60 70 80 90 100]; % number of realizations considered
% num_sim = 100; % number of simulations for every number of realizations in numReal_vec

%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for optimization

% Lower bounds:
noise_LB=0;
a_q_LB=1;
b_q_LB=T0;
c_q_LB=0;
c_r_LB=0;

% Upper bounds:
noise_UB=500;
a_q_UB=2000;
b_q_UB=Tf;
c_q_UB=10000;
c_r_UB=100000;

% Time vectors
t= time_vec;
s= t';
tau_Q = (t*ones(1,length(s))+ones(length(t),1)*s)./2;
tau_R = t*ones(1,length(s))-ones(length(t),1)*s;
tauQ = diag(tau_Q);
tauR = diag(flipud(tau_R));

%% 2. Example of realizations from this parametric setting
numReal_example = 3;
[X_example,C_true,R_true,Q_true] = lspdata_withNoise_simulation(numReal_example,a_q_true,b_q_true,c_q_true,c_r_true,noise_true,N,time_vec);

figure
contour(C_true)
view([90 90])
colorbar
title ('true covariance')

q_true =  noise_true+a_q_true*exp(-(c_q_true/2).*((tauQ-b_q_true).^2));
r_true = exp(-(c_r_true/8).*tauR.^2);

figure
subplot(221)
plot(tauQ,q_true)
axis 'tight'
xlabel('\eta (s)')
title('q')
subplot(222)
plot(tauR,r_true)
xlabel('\tau (s)')
title('r')
subplot(2,2,[3 4])
plot(time_vec,X_example)
xlabel('time (s)')
title('example realizations')
axis 'tight'
sgtitle('LSP model')


%% 3. Inference
% Allocate matrices to save distance from true covariance 
dist_HATS = zeros(length(numReal_vec),num_sim); 
dist_SCM = zeros(length(numReal_vec),num_sim);
dist_LS_SCM = zeros(length(numReal_vec),num_sim);

% and parameter estimates
par_est_HATS = zeros(5,length(numReal_vec),num_sim);
par_est_SCM  = zeros(5,length(numReal_vec),num_sim);

%%%%%%%%%%%%%%%%%%%%%%

for indReal = 1: length(numReal_vec)
    
    numReal = numReal_vec(indReal);
    
    for sim = 1: num_sim
        
        X = lspdata_withNoise_simulation(numReal,a_q_true,b_q_true,c_q_true,c_r_true,noise_true,N,time_vec);
        
        % Starting points for optimization
        
        a_q_0 = (a_q_UB+a_q_LB)/2; 
        b_q_0 = (T0+Tf)/2;
        % noise, c_q and c_r are uniform random in their range [LB, UB]
        noise_0 = noise_LB + (noise_UB-noise_LB)*rand; 
        c_q_0 = c_q_LB + (c_q_UB-c_q_LB)*rand;
        c_r_0 = c_r_LB + (c_r_UB-c_r_LB)*rand;
        
        all_AIP = X.^2; % all_AIP (:,i) contains the average istanteneous power of realization i
        % Average istantaneous power for the process 
        meanAIP = mean (all_AIP,2);
        
        % Estimate the lambda parameters by non-linear least square fitting of the function
        qfo = fitoptions('Method','NonlinearLeastSquares','Lower',[noise_LB a_q_LB b_q_LB c_q_LB],'Upper',[noise_UB a_q_UB b_q_UB c_q_UB],'StartPoint',[noise_0 a_q_0 b_q_0 c_q_0]);
        qft = fittype(' L + a* exp(-(c/2) * ((x-b).^2))','options',qfo);
        q_fit = fit(tauQ,meanAIP,qft); % the method is least square

        % Round off the estimates to the second decimal
        noise_HATS = round (q_fit.L * 100) / 100;
        a_q_HATS = round(q_fit.a * 100) / 100;
        b_q_HATS = round(q_fit.b * 100) / 100;
        c_q_HATS = round(q_fit.c * 100) / 100;
        
        % Generate the matrix Q 
        Q_est = noise_HATS + a_q_HATS * exp(-(c_q_HATS/2)*((tau_Q-b_q_HATS).^2));

        % Sample covariance estimate:

        if (numReal>1)
        C_SCM = cov(X');
        else
            C_SCM = X*X'; % for just one realization Matlab function cov would compute just the variance and not the cov. matrix 
        end
        
        % Estimate R using the fact that R is stationary
        R_HATS = zeros(N,N);
        R_SC = C_SCM./(Q_est);
        for k=1: N % we have only N diagonals to consider since the SCM is already symmetric,
            % and the index k represents how long is the diagonal
            pos = N-k; % position of the diagonal (i.e. shift from the main diagonal)
            mean_diag = mean(diag(R_SC,pos));    
            aux = diag(mean_diag*ones(1,k),pos); 
            R_HATS = R_HATS+aux;       
        clearvars pos mean_diag aux
        end
        % Make it symmetric
        aux = R_HATS';
        aux(logical(eye(size(aux)))) = 0;
        R_HATS = R_HATS+aux;

        r_antidiag_HATS = diag(flipud(R_HATS));
        R_fit_HATS = r_antidiag_HATS ;

        % Estimate the rho parameters by non-linear least square fitting of the function
        rfo = fitoptions('Method','NonlinearLeastSquares','Lower',[c_q_HATS],'Upper',[c_r_UB],'StartPoint',[c_r_0]);
        rft = fittype('exp(-(c/8).*x.^2)','options',rfo);
        r_HATS = fit(tauR,R_fit_HATS,rft);
        c_r_HATS = round(r_HATS.c * 100) / 100; % round off the estimates to the second decimal

        % Generate the matrix R and final covariance estimate
        final_R_HATS = exp(-(c_r_HATS/8).*(tau_R).^2); 
        C_tot_HATS = (final_R_HATS.*Q_est);
        
        % Inference using SCM
        
        fun_SCM = @(x) frobDist(x(1),x(2),x(3),x(4),x(5),time_vec,C_SCM)
        % x = [L, a_q, b_q, c_q, c_r]     
        A = [0,0,0,1,-1]; b = 0;  % costraint: c_q - c_r < 0, x(4)-x(5)<0
        Aeq = []; beq = [];
        theta0 = [noise_0,a_q_0,b_q_0,c_q_0,c_r_0];  % starting values for optimization
        lb = [noise_LB,a_q_LB,b_q_LB,c_q_LB,c_r_LB];
        ub = [noise_UB,a_q_UB,b_q_UB,c_q_UB,c_r_UB];
        nonlcon = [];
        temp_par_est_SCM = fmincon(fun_SCM,theta0,A,b,Aeq,beq,lb,ub,nonlcon);
        % round off the estimates to the second decimal:
        noise_SCM = round(temp_par_est_SCM(1)*100)/100;
        a_q_SCM = round(temp_par_est_SCM(2)*100)/100;
        b_q_SCM = round(temp_par_est_SCM(3)*100)/100;
        c_q_SCM = round(temp_par_est_SCM(4)*100)/100;
        c_r_SCM = round(temp_par_est_SCM(5)*100)/100;
        
        Q_LS_SCM = noise_SCM+a_q_SCM*exp(-(c_q_SCM/2)*((tau_Q-b_q_SCM).^2));
        R_LS_SCM = exp(-(c_r_SCM/8).*(tau_R).^2);  
        C_LS_SCM = (R_LS_SCM.*Q_LS_SCM); 
        
        % Frobenius distance between estimated covariance and true one
        S2_SCM = sqrt(sum(sum((C_SCM-C_true).^2)));
        S2_HATS = sqrt(sum(sum((C_tot_HATS-C_true).^2)));
        S2_LS_SCM = sqrt(sum(sum((C_LS_SCM-C_true).^2))); % in the paper, this is LS SCM
        
        dist_HATS(indReal,sim) = S2_HATS;
        dist_SCM(indReal,sim) = S2_SCM;
        dist_LS_SCM(indReal,sim)= S2_LS_SCM;
        
        par_est_HATS(:,indReal,sim) = [noise_HATS,a_q_HATS,b_q_HATS,c_q_HATS,c_r_HATS];
        par_est_SCM(:,indReal,sim) = [noise_SCM,a_q_SCM,b_q_SCM,c_q_SCM,c_r_SCM];
    
    end
    
end

%% 4. Distance from true covariance
% Normalize distance from true covariance, take mean and std 
FbN = sqrt(sum(diag(C_true'*C_true))); % Frobenius norm of the matrix C_true

dist_SCM_norm = dist_SCM./FbN;
dist_LS_SCM_norm = dist_LS_SCM./FbN;
dist_HATS_norm = dist_HATS./FbN;

norm_mean_dist_SCM = mean(dist_SCM_norm,2);
norm_mean_dist_LS_SCM = mean(dist_LS_SCM_norm,2);
norm_mean_dist_HATS = mean(dist_HATS_norm,2);

norm_std_dist_SCM = std(dist_SCM_norm,0,2);
norm_std_dist_LS_SCM = std(dist_LS_SCM_norm,0,2);
norm_std_dist_HATS = std(dist_HATS_norm,0,2);

SE_dist_SCM = norm_std_dist_SCM./sqrt(num_sim);
SE_dist_LS_SCM = norm_std_dist_LS_SCM./sqrt(num_sim);
SE_dist_HATS = norm_std_dist_HATS./sqrt(num_sim);

% error reduction 
m1_error = norm_mean_dist_LS_SCM;
m2_error = norm_mean_dist_HATS;
RED = ((m1_error-m2_error)./m1_error )*100.0 

%% 5. Parameter estimates

mean_par_SCM =  mean(par_est_SCM,3);
std_par_SCM = std(par_est_SCM,0,3);
SE_par_SCM = std_par_SCM./sqrt(num_sim);
c_par_SCM = std_par_SCM./mean_par_SCM;

LB_SCM = mean_par_SCM-2*SE_par_SCM;
UB_SCM =  mean_par_SCM+2*SE_par_SCM;

mean_par_HATS = mean(par_est_HATS,3); % average on the simulations
std_par_HATS= std(par_est_HATS,0,3);
SE_par_HATS = std_par_HATS./sqrt(num_sim);
c_par_HATS = std_par_HATS./mean_par_HATS;

LB_HATS = mean_par_HATS-2*SE_par_HATS;
UB_HATS =  mean_par_HATS+2*SE_par_HATS;

%% 6. Figures

figure

plot (numReal_vec,norm_mean_dist_SCM,'-o',numReal_vec,norm_mean_dist_LS_SCM,'-+',numReal_vec,norm_mean_dist_HATS,'->')
legend('SCM','LS SCM','HATS','Location','northeast');
grid;
ylabel('NMSE');
xlabel('realizations')

refline_L = ones(length(numReal_vec))*noise_true;
refline_a_q = ones(length(numReal_vec))*a_q_true;
refline_b_q = ones(length(numReal_vec))*b_q_true;
refline_c_q = ones(length(numReal_vec))*c_q_true;
refline_c_r = ones(length(numReal_vec))*c_r_true;
title ('Normalized Mean-Square Error of the covariance matrix estimates')


figure
subplot(5,2,1)
plot(numReal_vec,mean_par_HATS(1,:),'b-',numReal_vec,UB_HATS(1,:),'r:',numReal_vec,LB_HATS(1,:),'r:',numReal_vec,refline_L,'k--')
axis 'tight'
title('L')
subplot(5,2,3)
plot(numReal_vec,mean_par_HATS(2,:),'b-',numReal_vec,UB_HATS(2,:),'r:',numReal_vec,LB_HATS(2,:),'r:',numReal_vec,refline_a_q,'k--')
title('a_q')
axis 'tight'
subplot(5,2,5)
plot(numReal_vec,mean_par_HATS(3,:),'b-',numReal_vec,UB_HATS(3,:),'r:',numReal_vec,LB_HATS(3,:),'r:',numReal_vec,refline_b_q,'k--')
title('b_q')
axis 'tight'
subplot(5,2,7)
plot(numReal_vec,mean_par_HATS(4,:),'b-',numReal_vec,UB_HATS(4,:),'r:',numReal_vec,LB_HATS(4,:),'r:',numReal_vec,refline_c_q,'k--')
title('c_q')
axis 'tight'
subplot(5,2,9)
plot(numReal_vec,mean_par_HATS(5,:),'b-',numReal_vec,UB_HATS(5,:),'r:',numReal_vec,LB_HATS(5,:),'r:',numReal_vec,refline_c_r,'k--')
title('c_r')
xlabel('number of realizations')
axis 'tight'

subplot(5,2,2)
plot(numReal_vec,mean_par_SCM(1,:),'b-',numReal_vec, UB_SCM(1,:),'r:',numReal_vec,LB_SCM(1,:),'r:',numReal_vec,refline_L,'k--')
axis 'tight'
title('L')
subplot(5,2,4)
plot(numReal_vec,mean_par_SCM(2,:),'b-',numReal_vec, UB_SCM(2,:),'r:',numReal_vec,LB_SCM(2,:),'r:',numReal_vec,refline_a_q,'k--')
axis 'tight'
title('a_q')
subplot(5,2,6)
plot(numReal_vec,mean_par_SCM(3,:),'b-',numReal_vec, UB_SCM(3,:),'r:',numReal_vec,LB_SCM(3,:),'r:',numReal_vec,refline_b_q,'k--')
axis 'tight'
title('b_q')
subplot(5,2,8)
plot(numReal_vec,mean_par_SCM(4,:),'b-',numReal_vec, UB_SCM(4,:),'r:',numReal_vec,LB_SCM(4,:),'r:',numReal_vec,refline_c_q,'k--')
axis 'tight'
title('c_q')
subplot(5,2,10)
plot(numReal_vec,mean_par_SCM(5,:),'b-',numReal_vec, UB_SCM(5,:),'r:',numReal_vec,LB_SCM(5,:),'r:',numReal_vec,refline_c_r,'k--')
axis 'tight'
title('c_r')
xlabel('number of realizations')

sgtitle('HATS (left) and SCM (right) parameter estimates')


% take last estimates and plot covariances

noise_HATS = mean_par_HATS(1,end);
a_q_HATS = mean_par_HATS(2,end);
b_q_HATS = mean_par_HATS(3,end);
c_q_HATS = mean_par_HATS(4,end);
c_r_HATS = mean_par_HATS(5,end);

noise_SCM = mean_par_SCM(1,end);
a_q_SCM = mean_par_SCM(2,end);
b_q_SCM = mean_par_SCM(3,end);
c_q_SCM = mean_par_SCM(4,end);
c_r_SCM = mean_par_SCM(5,end);

Q_HATS = noise_HATS+a_q_HATS*exp(-(c_q_HATS/2)*((tau_Q-b_q_HATS).^2));
Q_SCM = noise_SCM+a_q_SCM*exp(-(c_q_SCM/2)*((tau_Q-b_q_SCM).^2));
R_HATS = exp(-(c_r_HATS/8).*(tau_R).^2);
R_SCM = exp(-(c_r_SCM/8).*(tau_R).^2);

C_HATS = (R_HATS.*Q_HATS);
C_LS_SCM = (R_SCM.*Q_SCM);

figure
subplot(221)
contour(C_true)
xlabel ('time (s)')
title('true')
view([90 90])
set(gca,'XTickLabel',0.1:0.1:0.5,'YTickLabel',0.1:0.1:0.5)

subplot(222)
contour(C_SCM)
xlabel ('time (s)')
title('SCM')
view([90 90])
set(gca,'XTickLabel',0.1:0.1:0.5,'YTickLabel',0.1:0.1:0.5)
 
subplot(223)
contour(C_LS_SCM)
view([90 90])
ylabel ('time (s)')
title('LS SCM')
set(gca,'XTickLabel',0.1:0.1:0.5,'YTickLabel',0.1:0.1:0.5)

subplot(224)
contour(C_HATS)
view([90 90])
title('HATS')
set(gca,'XTickLabel',0.1:0.1:0.5,'YTickLabel',0.1:0.1:0.5)
sgtitle('covariance estimates')


