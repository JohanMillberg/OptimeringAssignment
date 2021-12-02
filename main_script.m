% If we want to maximize the expected return, we choose w1 to be 1, since
% r1 is the greatest rate of return.

% Variance: 4.01

% max sum(wi*ri)
% sum(wi) = 1
% wi >= 0

CoV = 10^-2*[4.01 -1.19 0.6 0.74 -0.21;...
    0 1.12 0.21 -0.54 0.55;...
    0 0 3.31 0.77 0.29;...
    0 0 0 3.74 -1.04;...
    0 0 0 0 2.6];

lambda = 1;
prev_weights = [0 0 0 0 0];
Weights = [1/5 1/5 1/5 1/5 1/5]';

k = 0;
while abs((norm(Weights) - norm(prev_weights))/norm(Weights)) > 0.001
    k = k + 1;
    grad_f = zeros(5,1);
    for i = 1:5
        sum_ = 0;
        for j = 1:5
            if i == j
                sum_ = sum_ + 2*CoV(i,i)*Weights(i);
            else
                sum_ = sum_ + CoV(i,j)*Weights(j);
            end
        end
        grad_f(i) = sum_;
    end

    hess_f = zeros(5,5);
    for i = 1:5
        for j = 1:5
            if i == j
                hess_f(i,i) = 2*CoV(i,i);
            else
                hess_f(i,j) = CoV(i,j);
            end
        end
    end

    g = sum(Weights) - 1;
    grad_g = [1;1;1;1;1];
    hess_g = zeros(5,5);

    grad_L = grad_f - lambda.*grad_g;
    hess_L = hess_f - lambda.*hess_g;

    p_vec = [hess_L, -grad_g; -(grad_g'), 0] \ [-grad_L; g];
    prev_weights = Weights;
    Weights = Weights + p_vec(1:5,1);
    lambda = lambda + p_vec(6,1);
end
Variance = Weights'*CoV*Weights;
disp(Variance)

r = 1e-2.*[13.0, 4.4, 12.1, 7.1, 11.7];
H = zeros(5,5);

rho = 0.2;
Aeq = [r;ones(1,5)];
beq = [rho;1];
lb = zeros(5,1);
ub = ones(5,1);
f = zeros(5,1);

w = quadprog(2.*CoV,[],[],[],Aeq,beq,lb,ub);
%%% Short selling allowed

Aeq = [ones(1,5)];
beq = [1];

alphas = 0.05:0.05:1;
sigmas = zeros(1,length(alphas));
r_bars = zeros(1,length(alphas));
w_s = zeros(5,length(alphas));

for i = 1:length(alphas)
    H = alphas(i).*(2.*CoV);
    f = (alphas(i) - 1).*r;
    w = quadprog(H,f,[],[],Aeq,beq,[],[]);
    w_s(:,i) = w;
    sigmas(i) = w'*CoV*w;
    r_bars(i) = r*w;
end
figure(1)
hold on
plot(sigmas,r_bars)
xlabel('Variance')
ylabel('Estimated rate of return')
title('Short selling allowed')
hold off


%%% Short selling not allowed
lb = zeros(5,1);
ub = ones(5,1);

alphas = 0.05:0.05:1;
sigmas = zeros(1,length(alphas));
r_bars = zeros(1,length(alphas));
w_s = zeros(5,length(alphas));

for i = 1:length(alphas)
    H = alphas(i).*(2.*CoV);
    f = (alphas(i) - 1).*r;
    w = quadprog(H,f,[],[],Aeq,beq,lb,ub);
    w_s(:,i) = w;
    sigmas(i) = w'*CoV*w;
    r_bars(i) = r*w;
end
figure(2)
hold on
plot(sigmas,r_bars)
xlabel('Variance')
ylabel('Estimated rate of return')
title('Short selling not allowed')
hold off