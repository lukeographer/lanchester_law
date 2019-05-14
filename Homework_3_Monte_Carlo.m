%Author: Luke Grantham
%REDID: 818681559
%Date: 2/10/2019
%For: Homework Assignment #2 CS558 Spring 2019.

%Setup
clear all;
close all;


dt = 0.01; %time delta
timeLim = 6; %time ending interval
numRuns = 10000;
simMat = zeros(timeLim * 100 + 2, 2, numRuns);
durationMat = zeros(1, numRuns);
xRemTroops = zeros(1, numRuns);
yRemTroops = zeros(1, numRuns);
alphaMin = 0.75;
alphaMax = 0.85;
betaMin = 0.85;
betaMax = 0.95;
x0 = 8000; %x starting troops
y0 = 7000; %y starting troops

xEventsSim = input("How many x reinforcement events? ");
txSim = x0 * input("Enter x % level which event occurs (between 0 and 1): ");

yEventsSim = input("How many y reinforcement events? ");
tySim = y0 * input("Enter y % level which event occurs (between 0 and 1): ");

%-----------------------------------------------------------------------

for i = 1:numRuns

    %Parameters
    time = 0; %time 
    x0 = 8000; %x starting troops
    y0 = 7000; %y starting troops
    alpha = (alphaMax-alphaMin).*rand(1,1) + alphaMin; %x lethality coefficient
    beta = (betaMax-betaMin).*rand(1,1) + betaMin; %y lethality coefficient
    x1 = 0.6; %x reinforcement percentage of original force
    y1 = 0.5; %y reinforcement percentage of original force
    alpha1 = (alphaMax-alphaMin).*rand(1,1) + alphaMin; %x reinforcement lethality coefficient
    beta1 = (betaMax-betaMin).*rand(1,1) + betaMin; %y reinforcement lethality coefficient
    
    %------------------------------------------------------------------------

    %Simulation variables
    t = time:dt:timeLim+dt; %simulation interval

    x = x0; %simulation variables
    y = y0; %simulation variables
    x1Troops = x1 * x; %number of x reinforcement troops 
    y1Troops = y1 * y; %number of y reinforcement troops

    xOccur = 0; %number of x reinforcements occured
    yOccur = 0; %number of y reinforcements occured
    
    xEvents = xEventsSim;
    tx = txSim;
    yEvents = yEventsSim;
    ty = tySim;

    mat = [0 0]; %matrix to hold simulation values
    numLoops = 1; %matrix indeces
    mat(numLoops, 1) = x0; %record x troops level
    mat(numLoops, 2) = y0; %record y troops level

    %Simulation 
    while time < timeLim %loop to simulate
      time = time + dt; %increment time
      numLoops = numLoops + 1; %index increment

      if (x <= tx && xEvents > 0) %evaluates when reinforcements are added
          xOccur = xOccur + 1; %number of x events incremented
          x = addTroops(x, x1Troops); %adds reinforcements
          alpha = updateLethality(alpha, alpha1, xOccur); %calculates lethality
          xEvents = xEvents - 1; %decremenets remaining events
      end
      if (y <= ty && yEvents > 0)%evaluates when reinforcements are added
          yOccur = yOccur + 1;%number of y events incremented
          y = addTroops(y, y1Troops); %adds reinforcements
          beta = updateLethality(beta, beta1, yOccur); %calculates lethality
          yEvents = yEvents - 1;%decremenets remaining events
      end

      if(x > 0 && y > 0) %model x troops
        dxdt = -beta * y; 
        dydt = -alpha * x;
        x = x + dxdt * dt; %calculate remaining X troops
        y = y + dydt * dt; %calculate remaining Y troops
      else
          if (x <= 0) 
              x = 0;
              dxdt = 0;
              break;
          end 
          if (y <= 0)
              y = 0;
              dydt = 0;
              break;
          end
      end
      


      mat(numLoops, 1) = x; %record x troops level
      mat(numLoops, 2) = y; %record y troops level
    end
    mat(numLoops:size(simMat), 1) = x;
    mat(numLoops:size(simMat), 2) = y;
    
    simMat(1:size(mat,1),:,i) = mat;
    durationMat(i) = time;
    xRemTroops(i) = x;
    yRemTroops(i) = y;
end

avg = mean(simMat,3);

figure, plot(t(1:uint16(max(durationMat)/dt)), avg(1:uint16(max(durationMat)/dt), :)); %plot results
legend('X Army','Y Army');
ylabel('Troop Count');
xlabel('Time');
dim = [.13 .7 .3 .3];
str = sprintf('x0 = %d, y0 = %d,',round(x0), round(y0));
annotation('textbox',dim,'String',str,'FitBoxToText','on');

figure,
subplot(2,2,1)
histogram(xRemTroops, 15, 'Normalization','probability');
title('X Remaining Troops PDF');
xlabel('Troops');
ylabel('Probability of Occurrence');

subplot(2,2,2)
cdfplot(xRemTroops);
title('X Remaining Troops CDF');
xlabel('Troops');
ylabel('Cummulative Probability');

subplot(2,2,3)
histogram(yRemTroops, 15, 'Normalization','probability');
title('Y Remaining Troops PDF');
xlabel('Troops');
ylabel('Probability of Occurrence');

subplot(2,2,4)
cdfplot(yRemTroops);
title('Y Remaining Troops CDF');
xlabel('Troops');
ylabel('Cummulative Probability');

figure, 
subplot(1,2,1)
histogram(durationMat, 15, 'Normalization','probability');
title('Duration PDF');
xlabel('Time');
ylabel('Probability of Occurrence');

subplot(1,2,2)
cdfplot(durationMat);
legend('Duration CDF','Location','best');
title('Duration CDF');



%------------------------------------------------------------------------

%Function Definitions
function[newTroops] = addTroops(existTroops, addedTroops)
    newTroops = existTroops + addedTroops;
end

function[newLethality] = updateLethality(existLethality, addedLethality, eventOccur)
    newLethality = (existLethality * eventOccur + addedLethality)/(2 + eventOccur - 1);
end