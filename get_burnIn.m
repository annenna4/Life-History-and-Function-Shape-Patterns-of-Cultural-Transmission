function [pop] = get_burnIn(pDeath,nPop,pMut,copyAll,copyThresholdHigh,copyThresholdLow,b,localMode,binSize)

% BURN IN RULE: it is assumed that the stationary state has been reached when all initially present variant types have gone extinct, ie. when all variant types have undergone neutral
% dynamic. This is a relatively strict rule that potentially takes a large amount of time

% ageIni = ceil(rand(1,nPop)*50); % initial birth years
% t = max(ageIni);
% 
% % initialisation of population 
% numbTypeIni = 2; % number of types initially present
% value = numbTypeIni;
% pop = get_popIni_pub(nPop,numbTypeIni,ageIni);
% 
% while  min(pop(1,:))<numbTypeIni+1 % until all types have undergone neutral dynamics
%     t = t+1;
%     if localMode == 0
%         [pop,value] = get_dynamics_pub(t,pop,value,pDeath,nPop,pMut,b,copyAll,copyThresholdHigh,copyThresholdLow,0,[]);
%     elseif localMode == 1
%         [pop,value] = get_dynamics_local_pub(t,pop,value,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,0,[],binSize);
%     end    
% end

% ALTERNATIVE BURN IN RULE: here it is assumed that that stationary state has been reached if two the heterogeneity index of two populations, both initialised with opposing initial conditions
% (one with maximum heterogeneity, one with minimum heterogeneity), have crossed.

ageIni = ceil(rand(1,nPop)*50); % initial birth years
t = max(ageIni);

%initialisation of population 1
numbTypeIni = 2; % number of types initially present
value1 = numbTypeIni;
pop1 = get_popIni(nPop,numbTypeIni,ageIni);

%initialisation of population 2
numbTypeIni = nPop/10; % number of types initially present
value2 = numbTypeIni;
pop2 = get_popIni(nPop,numbTypeIni,ageIni);

diffPop1Pop2 = 5;
stepTime = 1000;

while diffPop1Pop2>0
    t = t+1;

    if localMode == 0
        [pop1,value1] = get_dynamics(t,pop1,value1,pDeath,nPop,pMut,b,copyAll,copyThresholdHigh,copyThresholdLow,0,[]);
        [pop2,value2] = get_dynamics(t,pop2,value2,pDeath,nPop,pMut,b,copyAll,copyThresholdHigh,copyThresholdLow,0,[]);
    elseif localMode == 1
        [pop1,value1] = get_dynamics_local(t,pop1,value1,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,0,[],binSize);
        [pop2,value2] = get_dynamics_local(t,pop2,value2,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,0,[],binSize);
    end

    if mod(t,stepTime) == 0
        type = unique(pop1(1,:));
        h = hist(pop1(1,:),type)./nPop;
        divPop1 = sum(h.^2);
        type = unique(pop2(1,:));
        h = hist(pop2(1,:),type)./nPop;
        divPop2 = sum(h.^2);
        diffPop1Pop2 = divPop1-divPop2;
        if diffPop1Pop2 > 0.1
            stepTime = 500;
        elseif diffPop1Pop2 <0.1 && diffPop1Pop2 > 0.0
            stepTime = 100;
        end
    end
end

pop = pop1;
