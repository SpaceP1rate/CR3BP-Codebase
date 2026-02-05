function [value,isterminal,direction] = ode_event_xcross(~,x)
    value = x(2);        % y = 0, crosses x axis
    isterminal = 1;      % stop on first crossing
    direction = 0;       % any direction
end