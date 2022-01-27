%Matlab-implementation av Broyden-DFP med Armijo-linjesökning
%
%f(x1,x2) = (x1-2)^4 + (x1-2x2)^2
%starta i x0 = [0,3]

function [optimum] = QuasiNewtonDFP( )
    %definitoner
    beta=1/2; m=1; sigma=0.1;
    f = @(x) (x(1)-2)^4 + (x(1) - 2*x(2))^2;
    grad = @(x) [4*(x(1)-2)^3 + 2*(x(1)-2*x(2)); 8*x(2) - 4*x(1)];  %kolumnvektor
   
    armijocondition = @(x, d, m, s) f(x) - f(x + beta^m*s*d) >= -sigma*beta^m*s*grad(x)'*d;
    
    
    x = [0;3];
    alpha= 1;
    D = eye(2);
    vals = [];
    iterations = [];
    for i=1:10^3 
       d = -D*grad(x);
       
       if armijocondition(x,d,m,1)
           %startsteg uppfyllt, dubblera steglängden
           while armijocondition(x,d,m-1,1)
               m= m-1;
           end
       else
           %startsteg uppfyller ej armijovillkoret, 
           %halvera steglängden tills uppfyllt
           while ~armijocondition(x,d,m,1)
               m=m+1;
           end
       end
       alpha = beta^m*1;
       
       vals = [vals f(x)];
       iterations = [iterations x];
       %ny punkt
       xold = x;
       x = x + alpha*d;
        
       %Uppdatera hessianinversen Dk+1
       p = x - xold;
       q = grad(x) - grad(xold);
     
       %DFP-uppdateringen
       D = D + (p*p')/(p'*q) - (D*(q*q')*D)/(q'*D*q);
        
       
    end
       plot(1:length(vals),vals);
      title('DFP Funktionsvärden plottat mot iterationssekvensen');
      figure;
      plot(iterations(1,:), iterations(2,:));
      title('DFP Koordinaterna för iterationspunkterna');
   optimum=x;
    
end
