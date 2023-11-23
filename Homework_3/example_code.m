


% ------------------- The instruction of pcode usage ----------------------
% This script implements a single step of forward and backward propagation 
% for a randomly created datapoint for the digit classificatioin problem 
% with the provided pcode files. Feel free to test your own functions with
% below code

%% Clear cache
clear;clc;close all

%% Create a random seed
rng(1)

%% Create a random input and output data
% Input X is a 784 by 1 vector
X = rand(784,1); 

% Output Y is a 10 by 1 one hot vector
Y = zeros(10,1);
Y(5) = 1; % This mean the fifth element is the category, i.e., digit 4

%% Specify the layer dimension
% The layer_dims variable indicate the structure of the neural network
% Input layer: 784 neurons
% Hidden layer 1: 64 neurons
% Hidden layer 2: 64 neurons
% Output layer: 10 neurons
layer_dims = [784, 64, 64, 10];

%% Initilize the neural network parameters
% The parameters variable stores all the weights and biases values as a cell
% structure. For example, parameters{1}.W is the weights from input layer
% to hidden layer 1, and parameters{3}.b is the biases from hidden layer 2 to
% output layer
parameters = initialize_parameters(layer_dims);

%% Implement forward propagation
% The forward_pass variable stores all the activations and output values as
% a cell structure. For example, forward_pass{2} is the values of hidden
% layer neurons, forward_pass{4} is the values of output layer
forward_pass = forward_propagation(X, parameters);

%% Compute the cost
% The cost variable store the cross entropy loss between the predicted 
% output value forward_pass{4} and actual output value Y
cost = compute_cost(forward_pass{4}, Y);

%% Implement backward propagation
% The gradients variable stores all the gradients values of all parameters
% as a cell structure. For example, gradients{1}.dW is the gradients value
% of weights from input layer to the hidden layer 1, and gradients{3}.db is 
% the gradients value of biases from hidden layer 2 to the output layer.
gradients = backward_propagation(X, Y, parameters, forward_pass);

%% Implement gradient descent
% The parameters variable stores all the updated weights and biases values 
% as a cell structure after the gradient descent with a learning rate of 0.1
lr = 0.1;
parameters = update_parameters(parameters, gradients, lr);

