%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAB 1, Bayesian Decision Theory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Attribute Information for IRIS data:
%    1. sepal length in cm
%    2. sepal width in cm
%    3. petal length in cm
%    4. petal width in cm

%    class label/numeric label: 
%       -- Iris Setosa / 1 
%       -- Iris Versicolour / 2
%       -- Iris Virginica / 3
close all;
clear all;
clc;

%% this script will run lab1 experiments..
clear
load irisdata.mat

%% extract unique labels (class names)
labels = unique(irisdata_labels);

%% generate numeric labels
numericLabels = zeros(size(irisdata_features,1),1);
for i = 1:size(labels,1)
    numericLabels(find(strcmp(labels{i},irisdata_labels)),:)= i;
end

%% feature distribution of x1 for two classes
figure

subplot(1,2,1);hist(irisdata_features(find(numericLabels(:)==1),2),100), title('Iris Setosa, sepal width (cm)');
xlabel('Width (cm)');
ylabel('Count');
subplot(1,2,2), hist(irisdata_features(find(numericLabels(:)==2),2),100); title('Iris Veriscolour, sepal width (cm)');
xlabel('Width (cm)');
ylabel('Count');
figure

subplot(1,2,1), hist(irisdata_features(find(numericLabels(:)==1),1),100), title('Iris Setosa, sepal length (cm)');
xlabel('Length (cm)');
ylabel('Count');
subplot(1,2,2), hist(irisdata_features(find(numericLabels(:)==2),1),100); title('Iris Veriscolour, sepal length (cm)');
xlabel('Length (cm)');
ylabel('Count');   
figure

plot(irisdata_features(find(numericLabels(:)==1),1),irisdata_features(find(numericLabels(:)==1),2),'rs');
hold on;
plot(irisdata_features(find(numericLabels(:)==2),1),irisdata_features(find(numericLabels(:)==2),2),'k.');
axis([4 7 1 5]);
title('Sepal Length x_1 vs Sepal Width x_2');
xlabel('Sepal Length (cm)');
ylabel('Sepal Width (cm)');
legend('Iris Setosa', 'Iris Versicolour');
    

%% build training data set for two class comparison
% merge feature samples with numeric labels for two class comparison (Iris
% Setosa vs. Iris Veriscolour
trainingSet = [irisdata_features(1:100,:) numericLabels(1:100,1) ];


%% Lab1 experiments (include here)
disp('************EXPERIMENTS************');

%% Question 3
disp('************Feature x_2: Sepal Width************');
featureOfInterest = 2;  % sepal width
x1 = [3.3, 4.4, 5.0, 5.7, 6.3];
for k=1:numel(x1)
    x = x1(k);
    [posteriors_x,g_x] = lab1(x,trainingSet,featureOfInterest);
    disp(['x = ' num2str(x)]);
    disp(['Posterior Probability w1 = ' num2str(posteriors_x(1))]);
    disp(['Posterior Probability w2 = ' num2str(posteriors_x(2))]);
    disp(['Discriminant Function = ' num2str(g_x)]);
    
    if (g_x > 0)
        class_label = 1;
    else
        class_label = 2;
    end
    disp([labels(class_label)]);
end

%% Question 4
disp('************Threshold Per Feature************');
listOfFeatures = {'Sepal Length', 'Sepal Width', 'Petal Length', 'Petal Width'};
for k=1:4
    disp(['Threshold for ', listOfFeatures{k}]);
    th = threshold(trainingSet,k);
end

%% Question 6
disp('************Feature x_1: Sepal Length************');
featureOfInterest = 1;  % sepal length
x1 = [3.3, 4.4, 5.0, 5.7, 6.3];
for k=1:numel(x1)
    x = x1(k);
    [posteriors_x,g_x] = lab1(x,trainingSet,featureOfInterest);
    disp(['x = ' num2str(x)]);
    disp(['Posterior Probability w1 = ' num2str(posteriors_x(1))]);
    disp(['Posterior Probability w2 = ' num2str(posteriors_x(2))]);
    disp(['Discriminant Function = ' num2str(g_x)]);
    
    if (g_x > 0)
        class_label = 1;
    else
        class_label = 2;
    end
    disp([labels(class_label)]);
end