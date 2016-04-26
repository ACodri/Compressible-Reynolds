function [f_prime] = diff_2o_generic(x,f)
%
% Compute the derivative of the vector x by applying a second order scheme,
% at the boundary the scheme is of the first order (forward and backward)

% it works also with non uniform grid

N = length(x);
siz = size(x);
f_prime = zeros(siz);

f_prime(1) = (f(2) - f(1)) / (x(2) - x(1));

for i = 2 : N-1
    f_prime(i) = ( f(i+1) - f(i-1) ) / ( x(i+1) - x(i-1) );
end

f_prime(N) = (f(N) - f(N-1)) / (x(N) - x(N-1));


% % upwind
% for i = 1 : N
%     if i == 1
%         f_prime(i) = (f(2) - f(1)) / (x(2) - x(1));
%     elseif i > 1 && i < N
%         f_prime(i) = ( f(i) - f(i-1) ) / ( x(i) - x(i-1) );
%     elseif i == N
%         f_prime(i) = (f(N) - f(N-1)) / (x(N) - x(N-1));
%     end
% end

