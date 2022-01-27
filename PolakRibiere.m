function [ ] = PolakRibiere ()
    f = @(x) (x(1)-2)^4 + (x(1) - 2*x(2))^2;
    grad = @(x) [4*(x(1)-2)^3 + 2*(x(1)-2*x(2)); 8*x(2) - 4*x(1)]; 
    
    
    x = [0;3];
    beta = 0; d=0;
    i = 0;
    vals = [];
    iterations = [];
    while (norm(grad(x)) > eps)
        g = grad(x);
        d = -g + beta*d;
        
        h = @(a) f(x + a*d);
        alpha = fminbnd(h, 0, 10); % 0<alpha<10
        vals = [vals f(x)];
        iterations = [iterations x];
        x = x + alpha*d;
   
        g_old = g;
        g = grad(x);
        beta = (g'*(g-g_old))/(g_old'*g_old); 
        
        %reset past gradient information
        if mod(i,10)==0
            beta=0;
        end
        i=i+1;
       
    end
      plot(1:length(vals),vals);
      title('Funktionsvärden plottat mot iterationssekvensen');
      figure;
      plot(iterations(1,:), iterations(2,:));
      title('Koordinaterna för iterationspunkterna');
      disp("STOPP");disp(x);
end