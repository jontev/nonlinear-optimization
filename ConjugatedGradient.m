function [ ] = ConjugatedGradient ()
    
    Q = [1001, 1, 1, 1, 1; 1, 999, 1, 1, 1; 1,1,101, 1,1; 1,1,1,99,1;1,1,1,1,10];
    H = [1/1000, 0,0 ,0 ,0; 0, 1/1000, 0,0,0; 0,0,1/100,0,0;0,0,0,1/100,0;0,0,0,0,1/10];
    b = [1,1,1,1,1]';
    f = @(x)0.5*x'*Q*x - b'*x;
    grad = @(x) Q*x -b;
    
    
    
    %Lösning a), vanliga konjugerade gradientmetoden
    %ANalytisk linjesökning och plotta
    %1: f(xk) - f(x*) 2: norm(Qxk-b) logaritmisk skala
    
    %initialisering
    x = [0,0,0,0,0]';
    g = grad(x);
    B = 0;
    d = -g;
    alpha = -(g'*d)/(d'*Q*d);
    
    %Optimum måste vara känd på förhand...
    optimum = fminunc(f, x);
    %
    
    Results1 = zeros(1,6);
    Results2 = zeros(1,6);
    %iteration
    for i=1:6
       gold = grad(x);
       dold = d; 
       %fprintf('Just finished iteration #%d\n', i-1);disp(x);
       Results1(i) = f(x) - f(optimum);
       Results2(i) = norm(Q*x - b);
       x = x + alpha*d;
       %disp(x); disp(i);
       g = grad(x);
       B = (g'*g)/(gold'*gold);
       d = -g + B*dold;
       alpha = -(g'*d)/(d'*Q*d);
    end

    %Konvergerar på 5 iteration, dvs i x5
    plot(0:5, log(Results1)); %OPtimerar map de 3 största egenvärdena först
    title('Uppgift a) Logarimen av f(xk)-f(x*) för k=0,..,5');
    figure;
    plot(0:5, (Results2));
    title('Uppgift a)  ||Qxk-b|| för k=0,..,5');
    figure;
    disp("STOPP");disp(x);
    
    
    
    
    %Lösning b), konjugerade gradientmetoden med prekonditionering
    
     %initialisering
    x = [0,0,0,0,0]';
    g = grad(x);
    B = 0;
    d = -H*g;
    alpha = -(g'*d)/(d'*Q*d);
   
    Results1 = zeros(1,5);
    Results2 = zeros(1,5);
    %iteration
    for i=1:5
       gold = grad(x);
       dold = d;
       
       fprintf('Iteration #%d\n, close to optimal %d\n', i,norm(x-optimum ));disp(x);      
       Results1(i) = f(x) - f(optimum);
       Results2(i) = norm(Q*x - b);
       x = x + alpha*d;
      % disp(x); disp(i);
       g = grad(x);
       B = (g'*H*g)/(gold'*H*gold);
       d = -H*g + B*dold;
       alpha = -(g'*d)/(d'*Q*d);
    end
    plot(0:4, log(Results1)); %OPtimerar map de 3 största egenvärdena först
    title('Uppgift b) Logarimen av f(xk)-f(x*) för k=0,..,4');
    figure;
    plot(0:4, (Results2));
    title('Uppgift b) ||Qxk-b|| för k=0,..,4');
    figure;
    scatter(1:5,eig(Q));
    title('Egenvärdena för Q');
    figure;
    scatter(1:5,eig(H*Q));
    title('Egenvärdena för HQ');
    
    
    disp("STOPP b"); disp(x);
end