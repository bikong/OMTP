       %%*************************************************************************%%
       %%        PARTICLE SWARM OPRIMIZATION FOR SUPPORT VECTOR MACHINE           %%
       %%*************************************************************************%%
                              %%  File_name:psoSVMcg.m  %%
                              %%  Author: Bikong        %%


%% type == 3 for regression
%% type == 1 for classfication
function [bestCVmse,bestc,bestg,pso_option] = psoSVMcg(train_label,train,type,pso_option)

if nargin == 2
    pso_option = struct('c1',1.5,'c2',1.7,'maxgen',200,'sizepop',50, ...
                    'k',0.6,'wV',1,'wP',1,'v',5, ...
                     'popcmax',10^4,'popcmin',10^(-1),'popgmax',10^3,'popgmin',10^(-2));
end
% c1:acceleration constants
% c2:acceleration constants
% maxgen:Evolutionary Generation, Initial 200;
% sizepop:Population Size, Initial 20;
% k: maximum/minimum of velocity equals k times maximum/minimum population
% wV:inertia weight 
% wP:inertia weight for velocity in population update
% v: SVM Cross Validation
% popcmax:SVM Panality C Maximum, Initial 1000;
% popcmin:SVM Panality C Minimum, Initial 0.1; 
% popgmax:SVM Kernel g Maximum, Initial 1000; 
% popgmin:SVM Kernel g Minimum, Initial 0.01; 

Vcmax = pso_option.k*pso_option.popcmax;
Vcmin = -Vcmax ;
Vgmax = pso_option.k*pso_option.popgmax;
Vgmin = -Vgmax ;

eps = 10^(-3);

%% Initiate particle's population and velocity 
for i=1:pso_option.sizepop
    
    % getting population and velocity randomly 
    pop(i,1) = (pso_option.popcmax-pso_option.popcmin)*rand+pso_option.popcmin;  
    pop(i,2) = (pso_option.popgmax-pso_option.popgmin)*rand+pso_option.popgmin;
    V(i,1)= Vcmax*rands(1,1);  
    V(i,2)= Vgmax*rands(1,1);
    
    % Fitness
    cmd = ['-v ',num2str(pso_option.v),' -c ',num2str( pop(i,1) ),' -g ',num2str( pop(i,2) ),' -s ',num2str(type), ' -p 0.1'];
    fitness(i) = svmtrain(train_label, train, cmd);
end


[global_fitness bestindex]=min(fitness); % initiate extremal global fitness 
local_fitness=fitness;   % initiate extremal individual fitness

global_x=pop(bestindex,:);   % global extreme value 
local_x=pop;    % individual extreme value 

% average fitness of each generation
avgfitness_gen = zeros(1,pso_option.maxgen); 


for i=1:pso_option.maxgen
    
    for j=1:pso_option.sizepop
        
        % velocity update
        V(j,:) = pso_option.wV*V(j,:) + pso_option.c1*rand*(local_x(j,:) - pop(j,:)) + pso_option.c2*rand*(global_x - pop(j,:));
        if V(j,1) > Vcmax
            V(j,1) = Vcmax;
        end
        if V(j,1) < Vcmin
            V(j,1) = Vcmin;
        end
        if V(j,2) > Vgmax
            V(j,2) = Vgmax;
        end
        if V(j,2) < Vgmin
            V(j,2) = Vgmin;
        end
        
        % population update
        pop(j,:)=pop(j,:) + pso_option.wP*V(j,:);
        if pop(j,1) > pso_option.popcmax
            pop(j,1) = pso_option.popcmax;
        end
        if pop(j,1) < pso_option.popcmin
            pop(j,1) = pso_option.popcmin;
        end
        if pop(j,2) > pso_option.popgmax
            pop(j,2) = pso_option.popgmax;
        end
        if pop(j,2) < pso_option.popgmin
            pop(j,2) = pso_option.popgmin;
        end
        
        % Mutation
        if rand>0.5
            k=ceil(2*rand);
            if k == 1
                pop(j,k) = (20-1)*rand+1;
            end
            if k == 2
                pop(j,k) = (pso_option.popgmax-pso_option.popgmin)*rand + pso_option.popgmin;
            end            
        end
        
        %Fitness Value 
        cmd = ['-v ',num2str(pso_option.v),' -c ',num2str( pop(j,1) ),' -g ',num2str( pop(j,2) ),' -s 3 -p 0.1'];
        fitness(j) = svmtrain(train_label, train, cmd);
        

        %%_________________________________________________________________________________________-
        %Individual Best Fitness
        if fitness(j) < local_fitness(j)
            local_x(j,:) = pop(j,:);
            local_fitness(j) = fitness(j);
        end

        if fitness(j) == local_fitness(j) && pop(j,1) < local_x(j,1)
            local_x(j,:) = pop(j,:);
            local_fitness(j) = fitness(j);
        end        
        
        %Global Best Fitness
        if fitness(j) < global_fitness
            global_x = pop(j,:);
            global_fitness = fitness(j);
        end

        if abs( fitness(j)-global_fitness )<=eps && pop(j,1) < global_x(1)
            global_x = pop(j,:);
            global_fitness = fitness(j);
        end
        %__________________________________________________________________________________________
    end
    
    fit_gen(i)=global_fitness;    
    avgfitness_gen(i) = sum(fitness)/pso_option.sizepop;
end




figure;
hold on;
plot(fit_gen,'r*-','LineWidth',1.2);
plot(avgfitness_gen,'o-','LineWidth',1.2);
legend('Best','Average');
xlabel('Generation','FontSize',10);
ylabel('Fitness','FontSize',10);
grid on;

bestc = global_x(1);
bestg = global_x(2);
bestCVmse = fit_gen(pso_option.maxgen);

line1 = 'Particle Swarm Optimization';
line2 = ['(c1=',num2str(pso_option.c1), ...
    ',c2=',num2str(pso_option.c2),',Max generation=', ...
    num2str(pso_option.maxgen),', Population=', ...
    num2str(pso_option.sizepop),')'];
line3 = ['Optimized c=',num2str(bestc),' Optimized g=',num2str(bestg), ...
    ' Minimum MSE=',num2str(bestCVmse),];
title({line1;line2;line3},'FontSize',12);


