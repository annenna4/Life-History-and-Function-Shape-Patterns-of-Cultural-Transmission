function [pop] = get_popIni(nPop,numbTypeIni,ageIni)

pop = zeros(2,nPop);
pop(2,:) = ageIni;

h = ones(1,nPop);
for i = 1:numbTypeIni
    h((i-1)*nPop/numbTypeIni+1:i*nPop/numbTypeIni) = i;
end
pop(1,:) = h; 

