function [value,isterminal,direction] = poincare_neg_90L1(~,x)
    value = x(1)-1+0.012505;        % x = 0 crossing 
    isterminal = 1;
    direction = 1;       % detect crossing from positive to negative
end