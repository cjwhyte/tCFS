%% Minimal rivalry model dominance duration calculator

% Christopher Whyte 02/10/25

function [m_durationL, m_durationR, durationR, durationL, percept] = MinimalRivalry_DominanceDurations(X)
    
    gauss_smooth = @(x,n) conv(x,gausswin(n)./sum(gausswin(n)),'same');

    error_bound = .1;
    rateL = gauss_smooth(X(1,:),250);
    rateR = gauss_smooth(X(2,:),250);

    perceiveL = rateL > (rateR + error_bound);
    perceiveR = rateR > (rateL + error_bound);

    vals = find(perceiveL>0);
    a = diff(vals);
    b = find([a inf]>1);
    durationL = diff([0 b]);

    vals = find(perceiveR>0);
    a = diff(vals);
    b = find([a inf]>1);
    durationR = diff([0 b]);

    m_durationL = mean(durationL(3:end-1));
    m_durationR = mean(durationR(3:end-1));

    percept = -1*perceiveL + perceiveR;
    
end 
