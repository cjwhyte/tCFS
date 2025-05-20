function [X_store, H_store, input_store, percept] = tCFS_Simulator(p)

    % initial conditions
    X = zeros(2,1); H = zeros(2,1);
    L = .7; R = 1;

    % storage arrays
    X_store = zeros(2,p.sim_length);
    H_store = zeros(2,p.sim_length);
    percept = zeros(2,p.sim_length);
    input_store = zeros(2,p.sim_length);

    for t = 1:p.sim_length
        
        % stimulus 
        input = [L,R]';
    
        % numerical integration
        DW = p.sigma*sqrt(p.DT)*randn(2,1);
        DX = (-X + p.M.*max(p.W*X + input - p.g.*H,0))/p.tau;
        DH = (-H + X)./p.tau_H;
        X = X + p.DT*DX;
        H = H + p.DT*DH + DW;
    
        % dominance dependent target stimulus
        if X(1) > X(2) + p.percept_bound
           R = R + p.contrast_rate;
           percept(1,t) = 1;
        elseif X(2) > X(1) + p.percept_bound
           R = R - p.contrast_rate;
           percept(2,t) = 1;
        end 
        
        X_store(:,t) = X;
        H_store(:,t) = H;
        input_store(:,t) = input;
           
    end

end 