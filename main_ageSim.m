% Simulation of both versions of the age-structured neutral model described
% in the manuscript `Life History and Functionality in Cultural
% Transmission: Insight from Mismatches between Neutral Theory and Data' by
% Anne Kandler, Rafael D'Andrea, James O'Dwyer

clear all
close all

nPop = 10^5; % population size
pMut = 5*10^-3; % innovation rate (per transmission event)
pDeath = 0.03; % per capita death rate
b = 0; % strength of frequency-dependent transmssion, b<0 negative bias, b>0 positive bias

copyAll = 0; % if copyAll = 1 then copying happens from all age groups
copyThresholdHigh = 6; % upper bound of the age of the copying pool (ATTENTION: coincides with c_thresh+1 in the manuscript)
copyThresholdLow = 0; % lower bound of the age of the copying pool

tMax = 18; % time steps to be run after equilibrium has been reached
itMax = 10; % number of simulations

localMode = 1; % 0: age-structured neutral model as explained in section 2.1., 
               % 1: age-structured neutral model with local interactions as explained in section 2.2
binSize = 100; % size of local groups when localMode =1, MUST be divisor of nPop

PDmode = 0; % 0: do nothing, 1: calculate progeny distribution
freqMode = 0; % 0: do nothing, 1: records the frequencies of all cultural variant types in the interval [1,tMax] and their life times 
              % (ATTENTION: this option slows the simulation down greatly)

coeffPL = zeros(itMax,4); 
coeffPLpop = zeros(itMax,3);
mu_est = zeros(1,itMax);

for sim = 1:itMax

    fprintf('Simulation # %d\n',sim)

    % burn-in period for reaching stationarity
    fprintf('Burn-in period\n')
    [pop] = get_burnIn(pDeath,nPop,pMut,copyAll,copyThresholdHigh,copyThresholdLow,b,localMode,binSize);

    % re-name variant types for convenience
    names = unique(pop(1,:)); % unique names in pop
    for i = 1:numel(names)
        index = find(pop(1,:) == names(i)); % find indices of all instances of name names(i)
        pop(1,index) = ones(1,numel(index))*i; % re-name it with name i
    end
    value = numel(names);
    valueIni = value;
    tini = max(pop(2,:));

    % calculating the cultural composition of the population for tMax time
    % steps
    fprintf('Generating populations \n')

    namesFreq = zeros(1,value);
    freqTrait = zeros(tMax,value);
    ageVariant = [];

    for t = tini+1:tini+tMax

        if freqMode == 1 % recording the frequencies of all present variant types at t
            freqTrait = horzcat(freqTrait, zeros(tMax,value-size(freqTrait,2)));
            type = unique(pop(1,:)); % unique variant types in population
            h = hist(pop(1,:),type)./nPop; % and their frequencies
            freqTrait(t-tini,type) = h;
        end

        if localMode == 0 % age-structured neutral model as explained in section 2.1.
            [pop,value,namesFreq] = get_dynamics(t,pop,value,pDeath,nPop,pMut,b,copyAll,copyThresholdHigh,copyThresholdLow,PDmode,namesFreq);
        elseif localMode == 1 % age-structured neutral model with local interactions as explained in section 2.2
            [pop,value,namesFreq,addMut] = get_dynamics_local(t,pop,value,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,PDmode,namesFreq,binSize);
            effMut(t-tini) = addMut; % innovation rate in time step t
        end

    end

    % calculation of innovation rate of local model version
    if localMode == 1
        mu_est(sim) = mean(effMut);
    end

    % calculation life time of variants
    h = 1;
    if freqMode == 1
        for i = valueIni+1:size(freqTrait,2)
            idx = find(freqTrait(:,i));
            if idx(1)>1 && idx(end)<tMax % including only types where start and end date are known 
                ageVariant(h) = idx(end)-idx(1);
                h = h+1;
            end
        end
    end

    fprintf('Calculating VAD \n')
    types = unique(pop(1,:)); % variant types present in population
    freq = hist(pop(1,:),types); % frequencies of those types

    fprintf('Estimation of powerlaw - VAD \n')
    [alpha,xmin,L] = plfit(freq,'range',[1.01:0.01:3.01]); % using estimator developed by Clauset et al. (2009)
    coeffPLpop(sim,:) = [alpha,xmin,L];

    % plotting VAD
    type = unique(freq)';
    c = hist(freq,type)'./length(freq);
    c = [[type; type(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = [];
    maxFreqPop(sim) = c(end,1)./nPop; % frequency of most common variant type

    figure(1)
    loglog(c(:,1),c(:,2),'b--'); hold on;
    title('Variant abundance distribution')
    xlabel('Abundance x')
    ylabel('Probability P(x) of number of variants with abundance >= x')


    if PDmode == 1 % calculation of progeny distribution

        fprintf('Calculating PD\n')

        namesFreq = nonzeros(namesFreq);
        namesFreq = reshape(namesFreq,numel(namesFreq),1);
        sumProgeny = sum(namesFreq);

        fprintf('Estimation of powerlaw - PD \n')
        [alpha,xmin,L] = plfit(namesFreq,'range',[1.01:0.01:3.01]); % using estimator developed by Clauset et al. (2009)
        coeffPL(sim,1:3) = [alpha,xmin,L];
        coeffPL(sim,4) = sumProgeny;

        % plotting PD
        type = unique(namesFreq);
        c1 = hist(namesFreq,type)';

        c = c1./sum(c1);
        c = [[type; type(end)+1] 1-[0; cumsum(c)]];
        c(c(:,2)<10^-12,:) = [];

        figure(2)
        loglog(c(1:end,1),c(1:end,2),'b','LineWidth',1); hold on;
        title('Progeny distribution')
        xlabel('Abundance x')
        ylabel('Probability P(x) of number of variants with abundance >= x')

    end

    if freqMode == 1 

        fprintf('Plotting frequency time series and life time distribution\n')

        figure(3)
        for i = 1:size(freqTrait,2)
            plot([1:tMax],freqTrait(:,i)); hold on 
        end
        title('Frequency time series of cultural variants')
        xlabel('Time')
        ylabel('Frequencies of cultural variants')

        figure(4)
        [a,b] = ecdf(ageVariant);
        plot(b,a)
        title('Life time distribution of cultural variants')
        xlabel('Life time of cultural variant x')
        ylabel('Probability P(x) of number of variants with life time <= x')

    end
end

