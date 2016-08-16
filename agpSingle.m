function agpSingle
%Compute adaptive growth process


%VARIABLES
vars.w      = 40;    %Width of the system  (y direction)
vars.l      = 400;   %Length of the system (x direction)
vars.delta  = 0.1;   %Pheromone used to reduce strength of weak trails
vars.phe    = 1;     %Amount of pheromone left by ant after moving 
vars.d      = 0.01;  %Amount of pheromone subtracted from all cells
vars.lambda = 2;     %Value used in leaky integrator
vars.xMin   = 17;   
vars.xMax   = 24;
vars.maxIts = 10000;


tic
[P,Ss,Rs] = agpFunct(vars);
disp(sprintf('\n Computing time: %f',toc))


%PLOT STATISTICS
sProbs = 0;
sHist  = zeros(101,1);
sHist(1,1) = (vars.xMax - vars.xMin + 1)/vars.w;
%Compute stats
for iter = 1:vars.maxIts
    sProbs = sProbs + Ss(iter);
    if mod(iter,100) == 0
        sHist(iter/100+1,1) = sProbs/100;
        sProbs = 0;
    end
end
maxH = max(sHist,[],1);
disp(sprintf('\n Maximum prob of ants reaching target: %f\n',maxH))
%Plot of probability of reaching selected location 
figure
plot(0:100:vars.maxIts, sHist)
axis([0 vars.maxIts 0 1])
xlabel('Iteration','FontSize',12)
ylabel('Probability of ants reaching target in last 100 iterations',...
       'FontSize',12)
title('Ants reaching target over time','FontSize',14)
%plot reward value as a function of iteration 
figure
plot(Rs(1:1000,1))
xlabel('Iteration','FontSize',12)
ylabel('R value' ,'FontSize',12)
title('Rs for the first 1000 iterations','FontSize',14)


%PLOT PHEROMONES
%Plot pheromone levels as surface
figure
surf(P)
shading interp
xlabel('Length','FontSize',12)
ylabel('Width' ,'FontSize',12)
title('Pheromone levels of system in 3D','FontSize',14)
%Plot filtered pheromone levels in 2D
d         = 10;
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
title('Filtered pheromone levels of system in 2D','FontSize',14)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,Ss,Rs] = agpFunct(vars)


%GET VARIABLES AND DECLARE ARRAYS
w      = vars.w;
l      = vars.l;
delta  = vars.delta;
phe    = vars.phe;
d      = vars.d;
lambda = vars.lambda;
xMin   = vars.xMin;
xMax   = vars.xMax;
maxIts = vars.maxIts;

%Pheromones
P      = zeros(w,l);
%Probabilites
probs  = zeros(3);
%Reward value
R      = 0;
%Arrays to store stats
Ss     = zeros(maxIts,1);
Rs     = zeros(maxIts,1);


%MOVE ANTS
for iter = 1:maxIts
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
        x = mod(x - (rNum<cumProb(1)) + (rNum>cumProb(2)) - 1, w)+1;
        %Add pheromone to already travelled position  
        P(x,k) = P(x,k) + phe;
    end
    
    %Calculate reward function and store stats
    S        = (x > xMin) && (x < xMax);
    Ss(iter) = S;
    R        = R + (S-R)/lambda;
    Rs(iter) = R;
end