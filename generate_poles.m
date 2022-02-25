
function poles = generate_poles(number,option)
    
    if mod(number,2) ~= 0
        disp(['An odd number of poles was specified.  The number of ' ...
              'poles generated must be even.']);
        stop
    end
    
    if option == 1
        total = number/2;
        radius = 1;        
    elseif option == 2
        total = number;
        radius = 1;
    elseif option == 3
        total = number;
        radius = 0.7;
    elseif option == 4
        total = number;
        radius = 1;
    elseif option == 5
        total = number;
    end
    
    rs = zeros(1,total);
    thetas = zeros(1,total);
    
    if option < 4
        t0 = pi*(3-sqrt(5));
        
        for i=1:total
            rs(i) = radius*sqrt((i/(total+1)));
            thetas(i) = t0*i;
        end
        
        reals = zeros(1,total);
        imags = zeros(1,total);
        
        for i=1:total
            reals(i) = rs(i)*cos(thetas(i));
            imags(i) = rs(i)*sin(thetas(i));
        end
    elseif option == 4
        for k=1:(total/2)
           rk = sqrt(k/((total/2)+1));
           thetak = 2*sqrt(pi*k);
           reals(k) = rk*cos(thetak);
           imags(k) = rk*sin(thetak);
           %reals(2*k-1) = rk*cos(thetak);
           %imags(2*k-1) = rk*sin(thetak);
           %reals(2*k) = rk*cos(-thetak);
           %imags(2*k) = rk*sin(-thetak);
        end
    elseif option == 5
        reals = [];
        imags = [];
        for k=1:total
           rk = sqrt(k/(total+1));
           thetak = 2*sqrt(pi*k);
           if sin(thetak) > 0
               reals = [reals rk*cos(thetak)];
               imags = [imags rk*sin(thetak)];
               %reals = [reals;rk*cos(thetak);rk*cos(-thetak)];
               %imags = [imags;rk*sin(thetak);rk*sin(-thetak)];
           end
           k = k+1;
        end
    end
    
    if option == 1
        poles = reals + sqrt(-1)*imags;    
    elseif option == 2
        indeces = find(reals > 0);
        poles = reals(indeces) + sqrt(-1)*imags(indeces);
    elseif option == 3
        indeces = find(reals > 0);
        poles = reals(indeces) + sqrt(-1)*imags(indeces);
    elseif option == 4
        poles = reals + sqrt(-1)*imags;
    elseif option == 5
        poles = reals + sqrt(-1)*imags;
    end
    
    

    % for plotting:
    % scatter(real(poles),imag(poles));
    % hold on;
    % scatter(real(poles),-imag(poles));
    % hold off;
    % t = linspace(0,2*pi);
    % hold on;
    % plot(cos(t),sin(t),'k');
    % hold off;
    % stop

end