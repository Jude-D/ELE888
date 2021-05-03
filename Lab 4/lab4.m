
%% Set Up
%Following the steps from the lab manual
clc
clear all
close all
I= imread('house.tiff');
X = reshape(I, 256*256, 3);
X = double(X);
R = [0,0,1];G = [0,1,0];B = [1,0,0];
figure
plot3(X(:,1), X(:,2), X(:,3),'.','Color',[0.3, 0.7, 0])
title('All pixels in RGB');
xlabel('Red');
ylabel('Green');
zlabel('Blue');
figure
imshow(I)
title('Original image of house')
%% Part A when c = 2
D = [7.68 179.2 102.4; 128 204.8 51.2];
D_starting = [7.68 179.2 102.4; 128 204.8 51.2];
delta_D = zeros(size(D));
J = [];
D1 = [D_starting(1,:)];
D2 = [D_starting(2,:)];
iterations = 0;

while (and((delta_D ~= D),iterations < 100))
    delta_D = D;
    
    for q = 1:length(X)
       distanceToJ1(q,:) = X(q,:) - D(1,:);
       distanceToJ2(q,:) = X(q,:) - D(2,:);    
    end 
   
    cluster2 = sum(distanceToJ1.^2,2) > sum(distanceToJ2.^2,2);
    cluster1 = ~cluster2;

    D(1,:) = sum(X(cluster1, :)) / sum(cluster1);
    D(2,:) = sum(X(cluster2, :)) / sum(cluster2);
    
    J1 = sum((X - repmat(D(1,:), length(X), 1)).^2, 2);
    J2 = sum((X - repmat(D(2,:), length(X), 1)).^2,2);
    
    J = [J sum(min(J1, J2))];
    
    D1 = cat(1, D1, D(1,:));
    D2 = cat(1, D2, D(2,:));
    
    iterations = iterations+1;
    disp("iteration number: " + iterations)
    disp(D)
end

% i
close all
plot(J)
title('Error Criterion')
xlabel('Iteration')
ylabel('J')
grid;

% ii
figure
plot3(D1(:, 1), D1(:, 2), D1(:, 3), '-o')
hold on
plot3(D2(:, 1), D2(:, 2), D2(:, 3), '-o')
xlim([0 255]);
ylim([0 255]);
zlim([0 255]);
legend('D1','D2');
title('Cluster means');

grid;

% iii
figure
X1 = X(cluster1, :);
X2 = X(cluster2, :);
plot3(X1(:,1), X1(:,2), X1(:,3),'.','Color', D(1,:)/256)
hold on
plot3(X2(:,1), X2(:,2), X2(:,3),'.','Color', D(2,:)/256)
xlim([0 255]);
ylim([0 255]);
zlim([0 255]);
grid;
title('All pixels in RGB');
xlabel('Red');
ylabel('Green');
zlabel('Blue');


% iv
L = ones(length(X),3);

for i = 1 : length(X)
    if cluster2(i) == 1
       L(i,:) = D(2,:);
    else
       L(i,:) = D(1,:);
    end
end   

Ilabeled = reshape(L,256,256,3);
figure;
subplot(1,2,1);imshow(I)
title('original image')
subplot(1,2,2);imshow(uint8(Ilabeled))
title('labelled image')
sgtitle('Image in Labelled Form vs the Original Image');
%% Part B when c = 5, 1st case

c = 5;

recreated_image1=[];
Z = [90.9522485853271,55.1002680432840,71.7977762334955;119.999281006591,91.8594367884961,105.633672857999;166.502843339540,106.366664989180,96.4746615671078;163.744017526121,199.619691607685,220.143200202224;140.261463970379,154.625747650242,159.638279692395];

delta_Z = zeros(size(Z));
J = [];

while (and((delta_Z ~= Z),iterations < 100))
    delta_Z = Z;
    J = zeros(size(X,1), c);
        j1 = sum((X - repmat(Z(1,:), length(X), 1)).^2,2);
        j2 = sum((X - repmat(Z(2,:), length(X), 1)).^2,2);
        j3 = sum((X - repmat(Z(3,:), length(X), 1)).^2,2);
        j4 = sum((X - repmat(Z(4,:), length(X), 1)).^2,2);
        j5 = sum((X - repmat(Z(5,:), length(X), 1)).^2,2);
        J = [j1,j2,j3,j4,j5];
    [a, cluster] = min(J, [], 2);
    for i = [1:5]
        current_cluster = (cluster==i);
        Z(i, :) = sum(X(current_cluster, :)) / sum(current_cluster);
    end
    recreated_image1 = zeros(size(X));
    for i = [1:5]
        current_cluster = (cluster==i);
        recreated_image1 = (recreated_image1 + repmat(Z(i,:), length(X), 1) .* repmat(current_cluster, 1, width(X))/256);
    end
%         j1_cluster = (cluster == 1);
%         j2_cluster = (cluster == 2);
%         j3_cluster = (cluster == 3);
%         j4_cluster = (cluster == 4);
%         j5_cluster = (cluster == 5);
%         Q(1, :) = sum(X(j1_cluster, :)) / sum(j1_cluster);
%         Q(2, :) = sum(X(j2_cluster, :)) / sum(j2_cluster);
%         Q(3, :) = sum(X(j3_cluster, :)) / sum(j3_cluster);
%         Q(4, :) = sum(X(j4_cluster, :)) / sum(j4_cluster);
%         Q(5, :) = sum(X(j5_cluster, :)) / sum(j5_cluster);
%         Q = [j1_cluster, j2_cluster, j3_cluster, j4_cluster, j5_cluster]
    recreated_image1 = reshape(recreated_image1, 256, 256, 3);
    subplot(1,2,1);imshow(I)
    title('original image')
    subplot(1,2,2);imshow(recreated_image1)
    title('labelled image')
    sgtitle('Image in Labelled Form vs the Original Image');
end
disp('***Parameters for Part 1***')
disp(Z)
clusterZ = cluster;
%
figure
for i = [1:5]
    current_cluster = (cluster==i);
    Xi = X(current_cluster, :);
    plot3(Xi(:,1), Xi(:,2), Xi(:,3),'.','Color', Z(i,:)/256)
    hold all
    xlim([0 255]);
    ylim([0 255]);
    zlim([0 255]);
    grid;
    title('All pixels in RGB');
    xlabel('Red');
    ylabel('Green');
    zlabel('Blue');
end
%% Part B when c = 5, 2nd case
W = [173.087309789120,65.2262800888334,199.813031099228;74.0005303486661,57.3542478910000,172.885008831232;171.982890346039,170.965178115511,1.71912046553024;177.955967885245,216.164392070964,154.155644820940;17.4061487283227,88.1823772930668,99.0134257973721];
delta_W = zeros(size(W));
J = [];

while (and((delta_W ~= W),iterations < 100))
    delta_W = W;
    J = zeros(size(X,1), c);
        j1 = sum((X - repmat(W(1,:), length(X), 1)).^2,2);
        j2 = sum((X - repmat(W(2,:), length(X), 1)).^2,2);
        j3 = sum((X - repmat(W(3,:), length(X), 1)).^2,2);
        j4 = sum((X - repmat(W(4,:), length(X), 1)).^2,2);
        j5 = sum((X - repmat(W(5,:), length(X), 1)).^2,2);
        J = [j1,j2,j3,j4,j5];
    [a, cluster] = min(J, [], 2);
    for i = [1:5]
        current_cluster = (cluster==i);
        W(i, :) = sum(X(current_cluster, :)) / sum(current_cluster);
    end
    recreated_image2 = zeros(size(X));
    for i = [1:5]
        current_cluster = (cluster==i);
        recreated_image2 = (recreated_image2 + repmat(W(i,:), length(X), 1) .* repmat(current_cluster, 1, width(X))/256);
    end
%         j1_cluster = (cluster == 1);
%         j2_cluster = (cluster == 2);
%         j3_cluster = (cluster == 3);
%         j4_cluster = (cluster == 4);
%         j5_cluster = (cluster == 5);
%         Q(1, :) = sum(X(j1_cluster, :)) / sum(j1_cluster);
%         Q(2, :) = sum(X(j2_cluster, :)) / sum(j2_cluster);
%         Q(3, :) = sum(X(j3_cluster, :)) / sum(j3_cluster);
%         Q(4, :) = sum(X(j4_cluster, :)) / sum(j4_cluster);
%         Q(5, :) = sum(X(j5_cluster, :)) / sum(j5_cluster);
%         Q = [j1_cluster, j2_cluster, j3_cluster, j4_cluster, j5_cluster]
    recreated_image2 = reshape(recreated_image2, 256, 256, 3);
    subplot(1,2,1);imshow(I)
    title('original image')
    subplot(1,2,2);imshow(recreated_image2)
    title('labelled image')
    sgtitle('Image in Labelled Form vs the Original Image');
end
disp('***Parameters for Part 1***')
disp(Z)
clusterW = cluster;
%
figure
for i = [1:5]
    current_cluster = (cluster==i);
    Xi = X(current_cluster, :);
    plot3(Xi(:,1), Xi(:,2), Xi(:,3),'.','Color', Z(i,:)/256)
    hold all
    xlim([0 255]);
    ylim([0 255]);
    zlim([0 255]);
    grid;
    title('All pixels in RGB');
    xlabel('Red');
    ylabel('Green');
    zlabel('Blue');
end
%% Compare W and Z
    subplot(1,2,1);imshow(recreated_image1)
    title('Z image')
    subplot(1,2,2);imshow(recreated_image2)
    title('W image')
    sgtitle('Z image vs W image');

%% Part C

N = size(X,1);
XB1 = 0;
for i = [1:5]
    current_cluster = (clusterZ==i);
    Xk1 = X(current_cluster, :);
    mu_j = sort(sum((Z - repmat(Z(i,:), c, 1)).^2, 2).^.5);
    XB1 = XB1 + sum(sum((Xk1 - repmat(Z(i,:), length(Xk1), 1)).^2, 2).^.5) / mu_j(2);
end

XieBeni1 = XB1 / N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XB2 = 0;
for i = [1:5]
    current_cluster = (clusterW==i);
    Xk2 = X(current_cluster, :);
    mu_j = sort(sum((W - repmat(W(i,:), c, 1)).^2, 2).^.5);
    XB2 = XB2 + sum(sum((Xk2 - repmat(W(i,:), length(Xk2), 1)).^2, 2).^.5) / mu_j(2);
end

XieBeni2 = XB2 / N