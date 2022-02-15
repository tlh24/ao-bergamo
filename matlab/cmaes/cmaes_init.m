function [] = cmaes_init()

global cm;

cm.nVar=97;                % Number of Unknown (Decision) Variables

cm.VarSize=[1 cm.nVar];       % Decision Variables Matrix Size

cm.VarMin=-0.05;             % Lower Bound of Decision Variables
cm.VarMax= 0.05;             % Upper Bound of Decision Variables


% Maximum Number of Iterations
cm.MaxIt=2000;

% Population Size (and Number of Offsprings)
cm.lambda=(4+round(3*log(cm.nVar)))*10;

% Number of Parents
cm.mu=round(cm.lambda/2);

% Parent Weights
% default is log-distributed; we might want something else, given the noise
% in the evaluations.
cm.w = log(cm.mu+0.5)-log(1:cm.mu); 
cm.w = (cm.mu - (1:cm.mu)).^3.5; 
cm.w = cm.w/sum(cm.w);

% Number of Effective Solutions
cm.mu_eff = 1/sum(cm.w.^2);

% Step Size Control Parameters (c_sigma and d_sigma);
cm.sigma0 = 0.3*(cm.VarMax - cm.VarMin);
cm.cs = (cm.mu_eff + 2) / (cm.nVar + cm.mu_eff + 5);
cm.ds = 1 + cm.cs + 2 * max( sqrt((cm.mu_eff-1)/(cm.nVar+1))-1,0 );
cm.ENN = sqrt(cm.nVar) * (1-1/(4*cm.nVar)+1/(21 * cm.nVar^2)); % where does that 21 come from? 

% Covariance Update Parameters
cm.cc = (4 + cm.mu_eff / cm.nVar) / (4 + cm.nVar + 2 * cm.mu_eff / cm.nVar);
cm.c1 = 2 / ( (cm.nVar+1.3)^2 + cm.mu_eff);
alpha_mu = 2;
cm.cmu = min(1-cm.c1, alpha_mu * (cm.mu_eff - 2 + 1/cm.mu_eff) / ((cm.nVar + 2)^2 + alpha_mu * cm.mu_eff/2));
cm.hth = (1.4+2/(cm.nVar+1)) * cm.ENN;


cm.ps = cell(cm.MaxIt,1);
cm.pc = cell(cm.MaxIt,1);
cm.C = cell(cm.MaxIt,1);
cm.sigma = cell(cm.MaxIt,1);

cm.ps{1} = zeros(cm.VarSize);
cm.pc{1} = zeros(cm.VarSize);
cm.C{1} = eye(cm.nVar);
cm.sigma{1} = cm.sigma0;

cm.empty_individual.Position = [];
cm.empty_individual.Step = [];
cm.empty_individual.Cost = [];

cm.M = repmat(cm.empty_individual, cm.MaxIt,1);
cm.M(1).Position = unifrnd(cm.VarMin, cm.VarMax, cm.VarSize);
cm.M(1).Position = randn(cm.VarSize) * 0.01; % custom!! 
cm.M(1).Step = zeros(cm.VarSize);
cm.M(1).Cost = 1e6;

cm.BestSol = cm.M(1);

cm.BestCost = zeros(cm.MaxIt, 1);

cm.g = 0; % iteration number. 
