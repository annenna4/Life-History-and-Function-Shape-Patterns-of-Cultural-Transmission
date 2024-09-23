function [pop,value,namesFreq] = get_dynamics(t,pop,value,pDeath,nPop,pMut,b,copyAll,copyThreshHigh,copyThreshLow,PDmode,namesFreq)

% DEATH
indexDeath = [];
nBirth = binornd(nPop,pDeath); % number of death = number of birth
indexDeath = randsample(nPop,nBirth); % generating indices of individuals to be removed

% REPRODUCTION
nMut = binornd(nBirth,pMut); % number of innovations
if copyAll == 0
    copyIndex = find(pop(2,:)>(t-copyThreshHigh) & pop(2,:)<(t-copyThreshLow) ); % defining copy pool
    i = 1;
    while isempty(copyIndex) % trouble shooting if there are no individuals of the considered age (happens, if at all, only at the beginning of the simulation)
        copyIndex = find(pop(2,:)>(t-(copyThreshHigh+i)) & pop(2,:)<(t-(copyThreshLow-i)));
        i = i+1;
    end
else % copy pool = population 
    copyIndex = [1:nPop];
end

if b==0 % unbiased transmission 
    h = numel(copyIndex);
    index = randsample(h,nBirth-nMut,'true');
    hAdd = pop(1,copyIndex(index));
else % frequency-dependent transmission of strength b
    % choosing role models from copy pool
    types = unique(pop(1,copyIndex)); % variant types present in copy pool
    if numel(types)>1
        h = hist(pop(1,copyIndex),types); % frequencies of all variant types currently present
        h = (h./numel(copyIndex)).^(1+b);
        hAdd = randsrc(nBirth-nMut,1,[types,h./sum(h)])
    else 
        hAdd = types*ones(1,nBirth-nMut);
    end
end

pop(1,indexDeath) = [hAdd value + [1:nMut]]; % adding copied types + innovations
pop(2,indexDeath) = t*ones(1,nBirth); % adding birth years

if PDmode == 1 % collecting progeny contained in hAdd
    names = unique(hAdd); 
    if numel(names)>1
        [progFreq] = hist(hAdd,names);
        namesFreq(names) = namesFreq(names) + progFreq; % updating progeny count
    else
        namesFreq(names) = namesFreq(names) + numel(hAdd); % updating progeny count
    end
    namesFreq = [namesFreq [ones(1,nMut)]]; % adding innovations to progeny count
end

value = value + nMut;

