
%background
% X- matrice that each colum represents input
% M - matrice including 4 vectors representing the original corners

% Neriya Mizrahi (315711697)
% Yair Lahad (205493018)

 
clear; close all; clc;

%% Get the data
load PCA_data;

% Training parameteres
eta         = 1e-3;	% learning rate
W = rand(2,3)*0.1;	% Weights_mat structure
n_epochs    = 100;	% number of training epochs


% setting variables
alfa = acosd(W(1,:)*W(2,:)'/norm(W(1,:))*norm(W(2,:))); %angle between vectors
yLen=2;
y=zeros(yLen,1); % output vector

%% Setting plot
figure(1)
scatter3(X(1,:),X(2,:),X(3,:),'b.') %plotting inputs
grid on
hold on
plot3([10*W(1,1),0],[10*W(1,2),0],[10*W(1,3),0],'r','LineWidth',3) %first eigenVector
plot3([10*W(2,1),0],[10*W(2,2),0],[10*W(2,3),0],'m','LineWidth',3) %second eigenVector
% Headlines and plot text description
title(['||W1|| = ' num2str(norm(W(1,:))) '       ||W2||= ' num2str(norm(W(2,:)))   '        \angle(W1,W2)='  num2str(alfa)])
text(11.5*W(1,1),11.5*W(1,2),11.5*W(1,3),'W1')
text(12*W(2,1),12*W(2,2),12*W(2,3),'W2')
xlabel('x1'), ylabel('x2'),zlabel('x3')
% fitting M borders
c=M(:,1); 
M(:,1)=M(:,2);
M(:,2)=c;
rec = fill3(M(1,:),M(2,:),M(3,:),'k');
alpha(rec,0.3) 
hold off
pause(0.3)

%% Sanger rule
for epoch= 1:n_epochs
    %looping input examples
    for i=1:size(X,2) 
        dw = zeros(3,2);
        sum = zeros(3,1);
        y=W*X(:,i); %updating output
        % looping output neurons
        for l=1:size(y,1)
              sum = sum + y(l,:)*W(l,:)';
              dw(:,l)=(eta/epoch)*y(l,:)*((X(:,i)-sum)); % updating learning rule
        end
        W = W + dw';     %updates weights
    end
    x_hat = y(l,:)*W(l,:)'; % reconstruction input
    alfa = acosd(W(1,:)*W(2,:)'/norm(W(1,:))*norm(W(2,:))); %angle between eigenVectors

%% Plotting
figure(1)
scatter3(X(1,:),X(2,:),X(3,:),'b.')
grid on
hold on
plot3([10*W(1,1),0],[10*W(1,2),0],[10*W(1,3),0],'r','LineWidth',3)  %first eigenVector
plot3([10*W(2,1),0],[10*W(2,2),0],[10*W(2,3),0],'m','LineWidth',3)  %second eigenVector
% Headlines and plot text description
xlabel('x1'),ylabel('x2'),zlabel('x3')
text(11*W(1,1),11*W(1,2),11*W(1,3),'W1')
text(12*W(2,1),12*W(2,2),12*W(2,3),'W2')
title(['||W1|| = ' num2str(norm(W(1,:))) '       ||W2||= ' num2str(norm(W(2,:)))   '         \angle(W1,W2)=' num2str(alfa)])
% fitting M borders
rec = fill3(M(1,:),M(2,:),M(3,:),'k');
alpha(rec,0.3) 
hold off
pause(0.01)
end

