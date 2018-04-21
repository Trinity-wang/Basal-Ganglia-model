% =========================================================================
% Spike generator
% =========================================================================
function Spk = spkgenmat(t, N, r, alpha,n)
% N = 100;      % number of spkie trains
% r = 50        % firing rate (Hz)
% alpha = 0.5;  % correlation
dt=0.1;
% t = t./1000; % to use the correct units!
% dt = t(2)-t(1);
% T = zeros(1, N);
Spk = zeros(n,n);
% for i = 1:length(t)
    p=r*dt;
    U = rand;
    R = (U <= p); % the reference spike train
    
   for a=1:n
       for b=1:n
        if rand <= alpha
            T = R;
        else
            T = (rand <= p);
        end
      Spk(a,b) = sum(T);
       end
   end
end
