clc; clear; close all;

%data to post process
dataFileName = 'mcTestData';
filePath = ['..', filesep, 'lib', filesep];

load([filePath, dataFileName]);

numMC = size(data.X,3);

% Spaghetti Plot states
symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
title = ['Simulated System States | ', num2str(numMC), 'MC Runs'];
make_plots(data.time,data.X,symbols,title,false);

% Spaghetti Plot Measurements
symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
title = ['Simulated System Measurements | ', num2str(numMC), 'MC Runs'];
make_plots(data.time,data.Y,symbols,title,false);