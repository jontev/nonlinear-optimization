function [ ] = GaussNewton ()

   t = [1, 2, 4, 5, 8];
   y = [3, 4, 6, 11, 20];
   
   g = @(x) [ (x(1)*exp(x(2)*t(1))-y(1)), (x(1)*exp(x(2)*t(2))-y(2)), (x(1)*exp(x(2)*t(3))-y(3)), (x(1)*exp(x(2)*t(4))-y(4)), (x(1)*exp(x(2)*t(5))-y(5))]';

    function [grads] = gradg (x)
       grads = zeros(2,5); 
       for i=1:5
           grads(1,i) = exp(t(i)*x(2));
           grads(2,i) = t(i)*x(1)*exp(t(i)*x(1));
       end
    end
   x0 = [2.5, 0.25]';
    
   x = x0;
   for i=1:15
       disp(x);
      x = x-(gradg(x)*gradg(x)' ) \ (gradg(x)*g(x));
   end
   Res = zeros(1,5);
   for i=1:5
       Res(i) = x(1)*exp(x(2)*t(i));
   end
   scatter(t,y);
   hold on;
   
   plot(t, Res);
   title('Observationer och modell'); 
   disp("optimum"); disp(0.5*g(x)'*g(x));
   disp(x);
   f = @(x) 0.5*g(x)'*g(x);
   fminunc(f,x0)
end