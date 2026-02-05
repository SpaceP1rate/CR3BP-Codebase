function [value,isterminal,direction] = ode_event_xcrossneg(~,x)
    value = x(2);        % x = 0 crossing
    isterminal = 1;      % stop on first crossing
    direction = -1;       % detect crossing from positive to negative
end