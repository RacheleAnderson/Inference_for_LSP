function [X,C,R,Q]=lspdata_withNoise_simulation(num_sim, a_q, b_q, c_q, c_r,sigma_square,N,time);

% LSPDATA computes num_sim data-vector of a locally stationary process 
% with the model parameters specified, each trajectory of length N 
%   
% OUTPUT:
%   X: output data vector [N x num_sim]
%   X(:,i) represent the i-th realization
%   C: Covariance matrix
%   R: Toeplitz matrix (stationary part of covariance)
%   Q: Hankel matrix (time varying part of covariance)

% INPUT:
%   num_sim: number of simulated trajectories
%   a_q, b_q, c_q, c_r: Model parameters
%   sigma_square: constant noise level
%   N: one trajectory length
%   time: time vector of length N
%

t = time;
s=t';

tau_R = t*ones(1,length(s))-ones(length(t),1)*s;
tau_Q = (t*ones(1,length(s))+ones(length(t),1)*s)./2;

R = exp(-(c_r/8).*(tau_R).^2);
Q = sigma_square + a_q * exp(-(c_q/2).*((tau_Q-b_q).^2));

% Covariance matrix of the base-band lsp-process
C=(R.*Q);

% Realizations from the filtered white noise realization b
Noise= randn(N,num_sim); % Noise realization
c1= sqrtm(C);
X= real(c1)*Noise; %lsp-realization



    








    