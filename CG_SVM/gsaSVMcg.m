
       %%*************************************************************************%%
       %%     GRAVITATIONAL SEARCH ALGORITHM FOR SUPPORT VECTOR REGRESSION        %%
       %%*************************************************************************%%
                              %%  File_name:gsaSVMsg.m  %%
                              %%  Author: Bikong        %%


%% type == 3 for regression
%% type == 1 for classfication
function [bestMSE, bestc, bestg, gsa_option] = gsaSVMcg(train_label, train, type, gsa_option)


if nargin == 3
  gsa_option = struct('max_it',200, 'N', 50, ...
                      'dim',2, 'Rpower', 1,'Rnorm', 2, ...
                      'alfa',20,'GO',100, 'ElitistCheck', 1, ...
                      'v',5,'cmax',500,'cmin',0.1, ...
                      'gmax',10, 'gmin', 0.01);
end

% search/position boundary 
up = [gsa_option.cmax,gsa_option.gmax];
down = [gsa_option.cmin,gsa_option.gmin];

%  agent position random initialization 
for i = 1:gsa_option.dim
  high = up(i);
  low = down(i);
  X(:,i) = rand(gsa_option.N,1).*(high-low)+low;
end
	
% velocity initialization
V = zeros(gsa_option.N,gsa_option.dim);

%  iteration
for iteration=1:gsa_option.max_it


    %Checking allowable range. 
    for i=1:gsa_option.N 
    %%Agents that go out of the search space, are reinitialized randomly .
    	Tp=X(i,:)>up;Tm=X(i,:)<down;X(i,:)=(X(i,:).*(~(Tp+Tm)))+((rand(1,gsa_option.dim).*(up-low)+low).*(Tp+Tm));
    %%Agents that go out of the search space, are returned to the boundaries.
    %Tp=X(i,:)>up;Tm=X(i,:)<low;X(i,:)=(X(i,:).*(~(Tp+Tm)))+up.*Tp+low.*Tm;
    end

    %Evaluation of agents.
    for i = 1:gsa_option.N
 		%cmd = ['-v ',num2str(gsa_option.v),' -c ',num2str(X(i,1)),' -g ',num2str(X(i,2)),'-s 3 -p 0.1']; %damn
 		cmd = ['-v ',num2str(gsa_option.v),' -c ',num2str( X(i,1) ),' -g ',num2str( X(i,2) ),' -s ',num2str(type),' -p 0.01'];
      	fitness(i) = svmtrain(train_label, train, cmd);
    end

    [local_MSE,local_index] = min(fitness); %minimization.        
    
    % record bestMSE & bestindex every iteration 
    local_bestMSE(iteration) = local_MSE;
    mean_MSE(iteration) = mean(fitness);
    %local_bestindex(i,:) = X(local_index,:);
    
    if iteration == 1
      	bestMSE = local_MSE;
      	bestc = X(local_index,1);
      	bestg = X(local_index,2);
    end

    if local_MSE < bestMSE 
      	bestMSE = local_MSE;
      	bestc = X(local_index,1);
      	bestg = X(local_index,2);
    end

  
%BestChart=[BestChart bestMSE];
%MeanChart=[MeanChart mean(fitness)];


%This function calculates the mass of each agent. 

    Fmax=max(fitness); 
    Fmin=min(fitness); 
    Fmean=mean(fitness); 


    if Fmax == Fmin
       	M = ones(size(fitness,1),1);
    else
      	best=Fmin;worst=Fmax; 
      	M=(fitness-worst)./(best-worst); 
    end
    M=M./sum(M); 
    


%Calculation of Gravitational constant. 
%%%here, make your own function of 'G'
    G = gsa_option.GO*exp((-gsa_option.alfa)*(iteration/gsa_option.max_it)); 


%Calculation of accelaration in gravitational field. 
    final_per=2; %In the last iteration, only 2 percent of agents apply force to the others.

%%total force calculation
    if gsa_option.ElitistCheck==1 
      	kbest=final_per+(1-iteration/gsa_option.max_it)*(100-final_per); 
      	kbest=round(size(X,1)*kbest/100);
    else
      	kbest=gsa_option.N; 
    end
      	[Ms ds]=sort(M,'descend');

    for i=1:gsa_option.N
        E(i,:)=zeros(1,gsa_option.dim);
        for ii=1:kbest
            j=ds(ii);
            if j~=i
                R=norm(X(i,:)-X(j,:),gsa_option.Rnorm); %Euclidian distanse.
            for k=1:size(X,2) 
                E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^gsa_option.Rpower+eps));
                %note that Mp(i)/Mi(i)=1
            end
            end
        end
    end

%%acceleration
    a=E.*G; %note that Mp(i)/Mi(i)=1

%%Agent movement. 
    V=rand(size(X,1),size(X,2)).*V+a; 
    X=X+V; 

end 


%% plot
figure;
hold on;
plot(local_bestMSE,'r*-','LineWidth',1);
plot(mean_MSE,'o-','LineWidth',1);
legend('best','average');
xlabel('Iteration');
ylabel('Fitness');
grid on;
line1 = ' Gravitational Search Algorithm';
line2 = [' Optimized c = ',num2str(bestc),' Optimized g = ',num2str(bestg)....
		' Minimum MSE = ',num2str(bestMSE)];
title({line1;line2},'FontSize',11);


