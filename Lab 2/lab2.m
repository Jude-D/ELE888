%% lab2 function
% Solve for linear discriminant function

function [sol_vector, perception_criterion] = lab2(trainingSet,eta,theta,initial_a,max_iterations)
    a(:,1) = initial_a;

    % Augmentation
    y = [ones(1,length(trainingSet)); transpose(trainingSet)];

    % Normalize
    [~,col] = size(y);
    y_hat = [y(:,1:15) y(:,(col/2)+1:end)*(-1)];
    
    
    for k = 1:max_iterations+1
        general_ldf = transpose(a(:,k))*y_hat;
        misclassified_col = find(general_ldf <= 0);
        
        perception_criterion(k) = sum((-1)*transpose(a(:,k))*y_hat);
        gradient(:,k) = sum((-1)*y_hat(:,misclassified_col),2);
        
        termination_criterion = eta*gradient(:,k);
        if abs(termination_criterion) <= theta
            sol_vector = a(:,k);
            disp(['No significant change in weight vector at Iteration Step = ' num2str(k)]);
            disp(['w_0 = ' num2str(sol_vector(1))]);
            disp(['w_1 = ' num2str(sol_vector(2))]);
            disp(['w_2 = ' num2str(sol_vector(3))]);
            break
        end
        a(:,k+1) = a(:,k)-termination_criterion;
    end
    
    if abs(termination_criterion) > theta
        sol_vector = a(:,k); % Solution exceeds max iterations
        disp('No solution found');
    end
end