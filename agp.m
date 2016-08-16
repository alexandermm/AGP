function agp
%A model of an adaptive growth process from the paper:
%Adaptive growth process: a model inspired by Pask's ear


%VARIABLES
w      = 50;   %Width of the system  (y direction)
l      = 500;  %Length of the system (x direction)
delta  = 0.1;  %Amount of pheromone used to reduce strength of weak trails
phe    = 1;    %Amount of pheromone left by ant after moving 
d      = 0.01; %Amount of pheromone subtracted from all cells
lambda = 2;    %Value used in leaky integrator
xMin   = 18;   
xMax   = 32;
%Pheromones
P  = zeros(w,l);
%Probabilites
probs = zeros(3);
%Reward value
R = 0;
%Array to store stats
sProbs   = 0;
stats    = zeros(151,1);
%Amount based on random prob to reach wanted cells: 12/50
stats(1) = 0.24; 


%MOVE ANTS
tic
for iter = 1:15000
    %Decay pheromones
    P = P.*(1 - 1/(495*R+5));
    %Subtract d from all cells
    P = (P - d).*((P - d) > 0);
    %Put ant in random initial position from 1 to w
    x = ceil(rand*w);
    
    for k = 1:l-1
        %Put some pheromone in 3 nearest position below ant
        P(mod(x-2,w)+1, k+1) = P(mod(x-2,w)+1, k+1) + delta;
        P(x           , k+1) = P(x           , k+1) + delta;
        P(mod(x  ,w)+1, k+1) = P(mod(x  ,w)+1, k+1) + delta;
        
        %Move ant based on prob. dist. given by above positions
        probs(1) = P(mod(x-2,w)+1, k+1);
        probs(2) = P(x           , k+1);
        probs(3) = P(mod(x  ,w)+1, k+1);  
        cumProb  = cumsum(probs/sum(probs));
        rNum     = rand;
        x = mod(x - (rNum < cumProb(1)) + (rNum > cumProb(2)) - 1, w) + 1;
        %Add pheromone to already travelled position  
        P(x,k) = P(x,k) + phe;
    end
    
    %Reward function 
    S = (x > xMin) && (x < xMax);
    sProbs = sProbs + S;
    R = R + (S-R)/lambda;
    
    %Store and compute statistics
    if mod(iter,100) == 0
        stats(iter/100+1,1) = sProbs/100;
        sProbs = 0;
    end
end
toc


%PLOT STATISTICS
plot(0:100:15000,stats)
axis([0 15000 0 1])
stats(end,1)


%FILTER AND PLOT PHEROMONES
maxP      = max(max(P))
d         = 1;
filteredP = P - (P - d).*(P > d);
figure
[dummy,h] = contourf(filteredP);
colormap jet
cbar_handle = colorbar('location','eastoutside');
set(get(cbar_handle,'ylabel'),'string','Pheromone level','FontSize',10)
daspect([2 1 1])
set(h,'EdgeColor','none') 
xlabel('Length','FontSize',12)
ylabel('Width' ,'FontSize',12)
title('Pheromone levels in system','FontSize',14)