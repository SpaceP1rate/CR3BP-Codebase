function [value,isterminal,direction] = ode_event_ycross(~,x)
    value = x(1);        % x = 0 crossing
    isterminal = 1;      % stop on first crossing
    direction = 0;       % any direction
end