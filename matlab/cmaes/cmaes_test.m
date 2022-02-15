function [] = cmaes_test()
global cm; 
cmaes_init; 
cost = zeros(180, 1); 
for u = 1:2000
	pop = cmaes_run(cost); 
	for j = 1:cm.lambda
		cost(j) = Ackley(pop(j,:)); 
	end
end

figure;
% plot(BestCost, 'LineWidth', 2);
semilogy(cm.BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;