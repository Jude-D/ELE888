%% runlab2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAB 2, Linear Discriminant Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ELE888-06
% Raymond Paras
% Jude D'Souza

% Attribute Information for IRIS data:
%    1. sepal length in cm
%    2. sepal width in cm
%    3. petal length in cm
%    4. petal width in cm

%    class label/numeric label: 
%       -- Iris Setosa / 1 
%       -- Iris Versicolour / 2
%       -- Iris Virginica / 3

%% Clean up
close all;
clear all;
clc;

%% Extract unique labels (class names)
load irisdata.mat
labels = unique(irisdata_labels);

%% Generate numeric labels
numericLabels = zeros(size(irisdata_features,1),1);
for i = 1:size(labels,1)
    numericLabels(find(strcmp(labels{i},irisdata_labels)),:)= i;
end

%% Build training data set for two class comparison
% merge feature samples with numeric labels for two class comparison (Iris
% Setosa vs. Iris Veriscolour vs. Iris Virginia
trainingSet = [irisdata_features(1:150,:) numericLabels(1:150,1) ];

%% Seperate into different datasets per class
datasetA = trainingSet(find(trainingSet(:,5)==1),2:3);
datasetB = trainingSet(find(trainingSet(:,5)==2),2:3);
datasetC = trainingSet(find(trainingSet(:,5)==3),2:3);

%% Question 1 and 2
disp('********** Question 1/2 **********');
% From datasets A and B
trainingSet30AB = [datasetA(1:15,:); datasetB(1:15,:)];
testingSet70AB = [datasetA(16:end,:); datasetB(16:end,:)];

% Initialization
eta = 0.01;
theta = 0;
initial_a = [0; 0; 1]; % Initial a(0)
max_iterations = 300;

disp('********** 30% Training Samples **********');
[a_training, Jp_a30] = lab2(trainingSet30AB,eta,theta,initial_a,max_iterations);

disp(' ');
disp('********** 70% Testing Samples **********');
[a_test, Jp_a70] = lab2(testingSet70AB,eta,theta,a_training,max_iterations);

% x2 vs x3
x2 = [datasetA(:,1); datasetB(:,1)];
x3 = (-1)*(a_test(2)*x2 + a_test(1))/a_test(3);

disp(' ');
disp(['H = ' num2str((-1)*a_test(2)/a_test(3)) '*x + ' num2str((-1)*a_test(1)/a_test(3))]);

figure;
plot(irisdata_features(find(numericLabels(:)==1),2),irisdata_features(find(numericLabels(:)==1),3),'rs');
hold on;
plot(irisdata_features(find(numericLabels(:)==2),2),irisdata_features(find(numericLabels(:)==2),3),'k.');
hold on;
plot(x2,x3,'-');
title('Sepal Width x_2 vs Petal Length x_3');
xlabel('Sepal Width (cm)');
ylabel('Petal Length (cm)');
legend('Iris Setosa', 'Iris Versicolour');

% Perception Criterion vs. Iterations (30% Training)
figure;
iterations = 1:length(Jp_a30);
plot(iterations, Jp_a30);
[y x] = min(Jp_a30);
hold on;
plot(x,y,'*','MarkerSize',15);
title('Question 2: Perception Criterion Function vs. Iterations');
legend('30% Training Samples (Dataset A&B)')
xlabel('Iterations');
ylabel('Perception Criterion');

%% Question 3
disp(' ');
disp('********** Question 3 **********');
% From datasets A and B
trainingSet70AB = [datasetA(1:35,:); datasetB(1:35,:)];
testingSet30AB = [datasetA(36:end,:); datasetB(36:end,:)];

% Initialization
eta = 0.01;
theta = 0;
initial_a = [0; 0; 1]; % Initial a(0)
max_iterations = 300;

disp('********** 70% Training Samples **********');
[a_training, Jp_a70] = lab2(trainingSet70AB,eta,theta,initial_a,max_iterations);

disp(' ');
disp('********** 30% Testing Samples **********');
[a_test, Jp_a30] = lab2(testingSet30AB,eta,theta,a_training,max_iterations);

% Perception Criterion vs. Iterations (70% Training)
figure;
iterations = 1:1:length(Jp_a70);
plot(iterations, Jp_a70);
[y x] = min(Jp_a70);
hold on;
plot(x,y,'*','MarkerSize',15);
title('Question 3: Perception Criterion Function vs. Iterations');
legend('70% Training Samples (Dataset A&B)')
xlabel('Iterations');
ylabel('Perception Criterion');

%% Question 4
% From datasets B and C
disp(' ');
disp('********** Question 4 **********');
trainingSet30BC = [datasetB(1:15,:); datasetC(1:15,:)];
testingSet70BC = [datasetB(16:end,:); datasetC(16:end,:)];

% Initialization
eta = 0.01;
theta = 0.01;
initial_a = [0; 0; 1]; % Initial a(0)
max_iterations = 50000;

disp('********** 30% Training Samples **********');
[a_training, Jp_a30] = lab2(trainingSet30BC,eta,theta,initial_a,max_iterations);

disp(' ');
disp('********** 70% Testing Samples **********');
[a_test,Jp_a70] = lab2(testingSet70BC,eta,theta,a_training,max_iterations);

% x2 vs x3
x2 = [datasetB(:,1); datasetC(:,1)];
x3 = (-1)*(a_test(2)*x2 + a_test(1))/a_test(3);

disp(' ');
disp(['H = ' num2str((-1)*a_test(2)/a_test(3)) '*x + ' num2str((-1)*a_test(1)/a_test(3))]);

figure;
plot(irisdata_features(find(numericLabels(:)==2),2),irisdata_features(find(numericLabels(:)==2),3),'rs');
hold on;
plot(irisdata_features(find(numericLabels(:)==3),2),irisdata_features(find(numericLabels(:)==3),3),'k.');
hold on;
plot(x2,x3,'-');
title('Sepal Width x_2 vs Petal Length x_3');
xlabel('Sepal Width (cm)');
ylabel('Petal Length (cm)');
legend('Iris Versicolour', 'Iris Virginica');

% Perception Criterion vs. Iterations (30% Training)
figure;
iterations = 1:1:length(Jp_a30);
plot(iterations, Jp_a30);
[y x] = min(Jp_a30);
hold on;
plot(x,y,'*','MarkerSize',15);
title('Question 4: Perception Criterion Function vs. Iterations');
legend('30% Training Samples (Dataset B&C)')
xlabel('Iterations');
ylabel('Perception Criterion');

% From datasets A and B
trainingSet70BC = [datasetB(1:35,:); datasetC(1:35,:)];
testingSet30BC = [datasetB(36:end,:); datasetC(36:end,:)];

disp(' ');
disp('********** 70% Training Samples **********');
[a_training, Jp_a70] = lab2(trainingSet70BC,eta,theta,initial_a,max_iterations);

disp(' ');
disp('********** 30% Testing Samples **********');
[a_test,Jp_a30] = lab2(testingSet30BC,eta,theta,a_training,max_iterations);

% Perception Criterion vs. Iterations (70% Training)
figure;
iterations = 1:1:length(Jp_a70);
plot(iterations, Jp_a70);
[y x] = min(Jp_a70);
hold on;
plot(x,y,'*','MarkerSize',15);
title('Question 4: Perception Criterion Function vs. Iterations');
legend('70% Training Samples (Dataset B&C)')
xlabel('Iterations');
ylabel('Perception Criterion');

%% Question 5
disp(' ');
disp('********** Question 5 **********');
% From datasets A and B
trainingSet30AB = [datasetA(1:15,:); datasetB(1:15,:)];
testingSet70AB = [datasetA(16:end,:); datasetB(16:end,:)];

% Initialization
theta = 0.01;
initial_a = [0; 0; 1]; % Initial a(0)
max_iterations = 300;

eta = 0.1;
disp(['Using eta = ' num2str(eta)]);
disp('********** 30% Training Samples **********');
[a_training, Jp_a30] = lab2(trainingSet30AB,eta,theta,initial_a,max_iterations);
disp('********** 70% Testing Samples **********');
[a_test, Jp_a70] = lab2(testingSet70AB,eta,theta,a_training,max_iterations);

% Eta = 0.1 (30% Training)
figure;
iterations = 1:1:length(Jp_a30);
plot(iterations, Jp_a30);
[y x] = min(Jp_a30);
hold on;
plot(x,y,'*','MarkerSize',15);
title('Question 5: Eta = 0.1');
legend('30% Training Samples (Dataset A&B)')
xlabel('Iterations');
ylabel('Perception Criterion');

disp(' ');
eta = 0.001;
disp(['Using eta = ' num2str(eta)]);
disp('********** 30% Training Samples **********');
[a_training, Jp_a30] = lab2(trainingSet30AB,eta,theta,initial_a,max_iterations);
disp('********** 70% Testing Samples **********');
[a_test, Jp_a70] = lab2(testingSet70AB,eta,theta,a_training,max_iterations);

% Eta = 0.001 (30% Training)
figure;
iterations = 1:1:length(Jp_a30);
plot(iterations, Jp_a30);
[y x] = min(Jp_a30);
hold on;
plot(x,y,'*','MarkerSize',15);
title('Question 5: Eta = 0.001');
legend('30% Training Samples (Dataset A&B)')
xlabel('Iterations');
ylabel('Perception Criterion');
