function [samples] = cmaes_run(Costs)
% this is pipelined -- 
% take cost, and apply it to the *last* sample generated. 
% then emit a new sample. 
% designed for interfacing with pipelined hardware, e.g the miroscope.

global cm; 

g = cm.g; 

if g > 0 && g <= cm.MaxIt
	for i = 1:cm.lambda
		cm.pop(i).Cost = Costs(i); 
		if Costs(i) < cm.BestSol.Cost
			cm.BestSol=cm.pop(i);
		end
	end
	
	% Sort Population
	[~, SortOrder]=sort(Costs);
	cm.pop = cm.pop(SortOrder);

	% Save Results
	cm.BestCost(g) = cm.BestSol.Cost;

	% Display Results
	disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(cm.BestCost(g))]);

	% Update Mean
	cm.M(g+1).Step = 0;
	for j=1:cm.mu
		cm.M(g+1).Step = cm.M(g+1).Step + cm.w(j) * cm.pop(j).Step;
	end
	cm.M(g+1).Position = cm.M(g).Position + cm.sigma{g} * cm.M(g+1).Step;

	% Update Step Size
	cm.ps{g+1} = (1-cm.cs) * cm.ps{g} + ...
		sqrt(cm.cs * (2-cm.cs) * cm.mu_eff) *...
		cm.M(g+1).Step / chol(cm.C{g})';
	cm.sigma{g+1} = cm.sigma{g} * ...
		exp(cm.cs/cm.ds * (norm(cm.ps{g+1})/cm.ENN-1))^0.3;

	% Update Covariance Matrix
	if norm(cm.ps{g+1}) / sqrt(1-(1-cm.cs)^(2*(g+1))) < cm.hth
	  hs = 1;
	else
	  hs = 0;
	end
	delta = (1-hs)*cm.cc*(2-cm.cc);
	cm.pc{g+1} = (1-cm.cc) * cm.pc{g} + ...
	  hs*sqrt(cm.cc * (2-cm.cc) * cm.mu_eff) * cm.M(g+1).Step;
	cm.C{g+1} = (1 - cm.c1 - cm.cmu) * cm.C{g} + ...
	  cm.c1 * (cm.pc{g+1}' * cm.pc{g+1} + delta * cm.C{g});
	for j=1:cm.mu
	  cm.C{g+1} = cm.C{g+1} + cm.cmu * cm.w(j) * ...
			cm.pop(j).Step' * cm.pop(j).Step;
	end

	% If Covariance Matrix is not Positive Definite or Near Singular
	[V, E] = eig(cm.C{g+1});
	if any(diag(E) < 0)
	  E = max(E,0);
	  cm.C{g+1} = V * E / V;
	end
end

% now, generate new samples. 
g = g+1; 
% clear the population. 
cm.pop=repmat(cm.empty_individual,cm.lambda,1);

% make new samples
% first one is just the mean
cm.pop(1).Step = zeros(cm.VarSize);
cm.pop(1).Position = cm.M(g).Position;
samples(1, :) = cm.pop(1).Position;
% now vary things
for i = 2:cm.lambda
	cm.pop(i).Step = mvnrnd(zeros(cm.VarSize), cm.C{g});
	cm.pop(i).Position = cm.M(g).Position + cm.sigma{g} * cm.pop(i).Step;
	samples(i, :) = cm.pop(i).Position; 
end

cm.g = g;


