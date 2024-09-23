function [pop,value,namesFreq,addMut] = get_dynamics_local(t,pop,value,pDeath,nPop,copyAll,copyThresholdHigh,copyThresholdLow,PDmode,namesFreq,binSize)

nBins = nPop/binSize; % number of bins

% DEATH
indexDeath = [];
nBirth = binornd(nPop,pDeath); % number of death = number of birth
indexDeath = randsample(nPop,nBirth); % generating indices of individuals to be removed

% REPRODUCTION
hAdd = [];
hAddPlus = 0;

h = floor((nBirth)/nBins); % ``equal" number of naive individuals per bin
nBirthBinV(1:nBins) = ones(1,nBins)*h;
rest = mod(nBirth,nBins); % number of remaining naive individuals
h = randperm(nBins); h = h(1:rest);
nBirthBinV(h) = nBirthBinV(h)+ones(1,rest); % distribution of the remaining individuals across bins

% dividing population into bins
v = randperm(nPop); % random permutation of population
binM = reshape(v,binSize,nBins)';

if copyAll == 0
    copyIndex = find(pop(2,:)>(t-copyThresholdHigh) & pop(2,:)<(t-copyThresholdLow) ); % defining copy pool
    i = 1;
    while isempty(copyIndex) % trouble shooting if there are no individuals of the considered age (happens, if at all, only at the beginning of the simulation)
        copyIndex = find(pop(2,:)>(t-(copyThreshHigh+i)) & pop(2,:)<(t-(copyThreshLow-i)));
        i = i+1;
    end
else % copy pool = population
    copyIndex = [1:nPop];
end
types = unique(pop(1,copyIndex)); % unique variant types in the copy pool
h = hist(pop(1,copyIndex),types); % and their frequencies in the copy pool

for i = 1:nBins
    mask = unique(pop(1,binM(i,:))); % variant types not to be copied
    if numel(types)>1
        h1 = randsrc(nBirthBinV(i),1,[types;h./sum(h)])'; % choose variants from copy pool
        [sharedvals,~] = ismember(h1,mask); % checking whether there are ``forbidden" types
        idx = find(sharedvals == 1); % find those
        h1(idx) = []; % delete variants present in the bin
        hAdd = [hAdd h1]; % collecting types of naive individuals not present in the bins
        hAddPlus = hAddPlus+numel(idx); % adding number of innvations needed to total count
    else
        hAddPlus = hAddPlus+nBirthBinV(i);
    end
end

pop(1,indexDeath) = [hAdd value + [1:hAddPlus]]; % adding copied + innovated types
pop(2,indexDeath) = t*ones(1,nBirth); % adding birth date

if PDmode == 1
    names = unique(hAdd); % copied types
    if numel(names)>1
        [progFreq] = hist(hAdd,names);
        namesFreq(names) = namesFreq(names) + progFreq; % updating progeny count
    elseif numel(names)==1
        namesFreq(names) = namesFreq(names) + numel(hAdd); % updating progeny count
    end
    namesFreq(value+1:value+hAddPlus) = ones(1,hAddPlus); % adding innovations to progeny count
end

value = value+hAddPlus;
addMut = (hAddPlus)/nBirth; % innovation rate in this time step


