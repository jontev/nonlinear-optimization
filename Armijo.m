%Armijo-linjesökning som bygger på upprepade halveringar eller dubbleringar
% av ett initialt steg tills det aktuella steget är det maximala som
% uppfyller Armijo-kravet. (Givet att steget endast kan ökas med faktor 2

%Beräkna det initiala steget med en iteration av newtons metod 

%minimum antas i x*= [2,1] med f(x*) = 0'
%BL: 
%d_k = -grad(x_k)
%t_k från armijo
%x_k+1 = x_k + d_k*t_k
function [optimum] = Armijo( )
    beta=1/2; m=1; sigma=0.1;
    f = @(x) (x(1)-2)^4 + (x(1) - 2*x(2))^2;
    grad = @(x) [4*(x(1)-2)^3 + 2*(x(1)-2*x(2)); 8*x(2) - 4*x(1)];  
    H = @(x) [12*(x(1)-2)^2, -4; -4,8];
    
    armijocondition = @(x, d, m, s) f(x) - f(x + beta^m*s*d) >= -sigma*beta^m*s*grad(x)'*d;
     
    x = [0;3];
    t= 0;
    
    vals = [];
    iterations = [];
    while (norm(grad(x)) > (eps+0.001))
       d = -grad(x);
       %bestämma steglängd, initiala steget med en iteration av newtons
       s = t - (grad(x)'*d)/(d'*H(x)*d); 
       m=1;
       %om startsteget (s eller?) uppfyller armijovillkoret, loopa fram
       %steglängden genom att dubblera steglängden tills m+1 inte uppfyller
       %den
       
       if armijocondition(x,d,m,s)
           %startsteg uppfyllt, dubblera steglängden
           while armijocondition(x,d,m-1,s)
               m= m-1;
           end
       else
           %startsteg uppfyller ej armijovillkoret, 
           %halvera steglängden tills uppfyllt
           while ~armijocondition(x,d,m,s)
               m=m+1;
           end
       end
       t = beta^m*s;
       iterations = [iterations x];
       vals = [vals f(x)];
       %ny punkt
       x = x + t*d;
           
    end
    disp((vals));
    plot(1:length(vals), vals);
    title('Funktionsvärden plottat mot iterationsnummer');
    figure;
    plot(iterations(1,:), iterations(2,:));
    title('Koordinaterna för iterationspunkterna');
   disp(x);
    
end
