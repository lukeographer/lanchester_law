%Author: Luke Grantham
%REDID: 818681559
%Date: 4/13/2019
%For: Homework Assignment #2 MODIFIED FOR MONTE CARLO ANALYSIS CS558 Spring 2019.

clear all;
close all;

timeLim = 5; %time ending interval
numRuns = 10000;
simMat = zeros(timeLim * 100 + 2, 2, numRuns);
alphaMin = 0.75;
alphaMax = 0.85;
betaMin = 0.85;
betaMax = 0.95;
durationMat = zeros(1, numRuns);
xRemTroops = zeros(1, numRuns);
yRemTroops = zeros(1, numRuns);

for i = 1:numRuns
    time = 0; %time 
    dt = 0.01; %time delta

    x0 = 1000; %x starting troops
    y0 = 800; %y starting troops
    alpha = (alphaMax-alphaMin).*rand(1,1) + alphaMin; %x lethality coefficient
    beta = (betaMax-betaMin).*rand(1,1) + betaMin;%y lethality coefficient

    t = time:dt:timeLim+dt; %simulation interval

    x = x0; %simulation variables
    y = y0; %simulation variables
    mat = [0 0]; %matrix to hold simulation values
    numLoops = 1; %matrix indeces
    mat(numLoops, 1) = x0; %record x troops level
    mat(numLoops, 2) = y0; %record y troops level
    while time < timeLim %loop to simulate
      time = time + dt; %increment time
      numLoops = numLoops + 1; %index increment

      if(x > 0) %model x troops
        dxdt = -beta * y; 
        x = x + dxdt * dt; %calculate remaining X troops
      else
        dxdt = 0;
        x = 0;
        break;
      end

      if(y > 0) %model y troops
        dydt = -alpha * x;
        y = y + dydt * dt; %calculate remaining Y troops
      else
        dxdt = 0;
        y = 0;
        break;
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

