function [X_store, H_store] = MinimalRivalry_Simulator(p) 
    % numerical integration of minimal rivalry model using Euler-Maruyama
    
    % Christopher Whyte 02/10/25

    % initial conditions
    X = 0.1*rand(2,1); H = zeros(2,1); 
    
    % storage containers 
    X_store = zeros(2,p.sim_length);
    H_store = zeros(2,p.sim_length);

    for t = 1:p.sim_length
    
        DX = (-X + p.M.*max(p.W*X + p.input(:,t) - p.g.*H,0))/p.tau;
        DH = (-H + X)./p.tau_H;
        DW_H = p.sigma*sqrt(p.DT)*randn(2,1);

        X = X + p.DT*DX;
        H = H + p.DT*DH + DW_H;
        
        X_store(:,t) = X;
        H_store(:,t) = H;

    end 

end 

