

function [x] = threshold(Training_Data, featureOfInterest)

% x = individual sample to be tested (to identify its probable class label)
% featureOfInterest = index of relevant feature (column) in Training_Data 
% Training_Data = Matrix containing the training samples and numeric class labels
syms x
D=Training_Data;

% D is MxN (M samples, N columns = N-1 features + 1 label)
[M,N]=size(D);    
 
f=D(:,featureOfInterest);  % feature samples
la=D(:,N); % class labels

%% %%%%Prior Probabilities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hint: use the commands "find" and "length"

%disp('Prior probabilities:');
p_w1 = length(find(la(:)==1))/M; % p(w1)
p_w2 = length(find(la(:)==2))/M; % p(w2)
%% %%%%%Class-conditional probabilities%%%%%%%%%%%%%%%%%%%%%%%

w1_indices = find(D(:,5)==1); % Indices of Iris Setosa
w1_set = D(w1_indices,featureOfInterest);   % Get specific column of values
mean_w1 = mean(w1_set); % mean of the class conditional density p(x|w1)
std_w1 = std(w1_set); % Standard deviation of the class conditional density p(x|w1)

w2_indices = find(D(:,5)==2);   % Indices of Iris Versicolour
w2_set = D(w2_indices,featureOfInterest);   % Get specific column of values
mean_w2 = mean(w2_set); % mean of the class conditional density p(x|w2)
std_w2 = std(w2_set); % Standard deviation of the class conditional density p(x|w2)

% use the above mean, std and the test feature to calculate p(x|w1)
cp_x_w1 = (1/(sqrt(2*pi)*std_w1))*exp(-0.5*(((x-mean_w1)/std_w1)^2));

% use the above mean, std and the test feature to calculate p(x|w2)
cp_x_w2 = (1/(sqrt(2*pi)*std_w2))*exp(-0.5*(((x-mean_w2)/std_w2)^2));

%% %%%%%%Compute the posterior probabilities%%%%%%%%%%%%%%%%%%%%

%disp('Posterior prob. for the test feature');
p_x = cp_x_w1*p_w1 + cp_x_w2*p_w2;

% p(w1/x) for the given test feature value
pos_w1_x = (cp_x_w1*p_w1)/p_x; % p(w1|x)

% p(w2/x) for the given test feature value
pos_w2_x = (cp_x_w2*p_w2)/p_x; % p(w2|x)

%% Compute the Threshold

equations = pos_w1_x == pos_w2_x
solx = solve(equations, x);
threshold = double(solx)

