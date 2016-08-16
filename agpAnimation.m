function agpAnimation
%A model of an adaptive growth process from the paper:
%Adaptive growth process: a model inspired by Pask's ear
%This implementation in MATLAB was done by Alex P. Martinez


%VARIABLES
w      = 50;   %Width of the system  (y direction)
l      = 400;  %Length of the system (x direction)
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
stats    = zeros(101,1);
%Amount based on random prob to reach wanted cells: 12/50
stats(1) = 0.24; 



%Define the area to be recorded
scrsz = get(0,'ScreenSize');
hf = figure;
set(hf,'Position',[scrsz(3)*0.1 scrsz(4)*0.2 scrsz(3)*0.8 scrsz(4)*0.5]);
frameCounter = 1;

%MOVE ANTS
tic
for iter = 1:10000
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
    
    %FILTER AND PLOT PHEROMONES
    if mod(iter,50) == 0
        filteredP = P - (P - 1).*(P > 1);
        [dummy,h] = contourf(filteredP);
        colormap jet; daspect([2 1 1]); set(h,'EdgeColor','none')
        
        cbar_handle = colorbar('location','eastoutside');
        set(get(cbar_handle,'ylabel'),'string','Pheromone level','FontSize',10)
        xlabel('Length','FontSize',12)
        ylabel('Width' ,'FontSize',12)
        title('Chemical levels in system','FontSize',14)
        
        F(:,frameCounter) = getframe(hf); 
        frameCounter = 1+frameCounter;
        pause(0.1)
    end
end
toc


%PLAY MOVIE
FPS = 10;
movie(hf,F,1,FPS);
movie2avi(F,'agp','fps',FPS,'compression','none')

%PLOT STATISTICS
figure
plot(0:100:10000,stats)
axis([0 10000 0 1])
xlabel('Time step','FontSize',12)
ylabel('Fraction of ants reaching target' ,'FontSize',12)
title('Ants reaching target','FontSize',14)
        
finalFraction = stats(end,1)
maxP = max(max(P))