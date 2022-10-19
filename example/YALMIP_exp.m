% Define A_x , b_x and symmetric variable E
A_x = [1 ,0; -1 ,0; 0 ,1; 0, -1; 1 ,1];
b_x = [1; 1; 2; 2; 2];
E = sdpvar(2, 2);

% Define objective
objective = -logdet(E);

% Define constraints
% The operator >= denotes matrix inequality
constraints = [E >= 0];
for i=1: size (A_x ,1)
    constraints =[ constraints , ...
        A_x(i ,:)*E*A_x(i ,:)' <= b_x (i) ^2];
end

% Solve the YALMIP Problem
optimize ( constraints , objective )

% Display the result
disp (inv ( value (E)))