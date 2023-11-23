function parameters = update_parameters(parameters, gradients, learning_rate)
% Updates the parameters of a feedforward neural network using gradient descent
% parameters: learned parameters, a struct containing W1, b1, W2, b2, etc.
% gradients: gradients of the cost with respect to each parameter, a struct containing dW1, db1, dW2, db2, etc.
% learning_rate: learning rate for gradient descent
% returns: updated parameters, a struct containing W1, b1, W2, b2, etc.

    L = length(parameters);
    for l = 1:L
        parameters{l}.W = parameters{l}.W - learning_rate * gradients{l}.dW;
        parameters{l}.b = parameters{l}.b - learning_rate * gradients{l}.db;
    end
end