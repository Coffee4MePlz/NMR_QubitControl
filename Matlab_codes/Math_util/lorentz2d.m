function y = lorentz2d(x, x0, y, y0, half_width_half_max)
    % A function to plot a Lorentzian (a.k.a. Cauchy) distribution given a
    % space vector 'x', a position and a half width at half maximum.
    % The distribution is then scaled to the specified height.

    if nargin <=2
        errordlg(["At least three input arguments are required.";...
            "These must be in the order x_space, position, half_width.";...
            "An extra height argument can be added to change scale."],...
            "Input Error")
    end

    %y = 1 ./ (half_width_half_max.*...
    %    (1+((x-position)./half_width_half_max).^2));
    
    y = 0.5 .* half_width_half_max .*...
        ( (x-x0).^2 + (y-y0).^2 + half_width_half_max.^2 ).^(-1.5);
    
%     if nargin == 3
         height = max(y);
%     end
    
    y = y.*(height/max(y))/(2*pi);
    
end