function [x] = SQP()
 %parametrar
 c = [0.063 5.04 0.035 10 3.36];
 dl = [99/100 99/100 9/10 99/100];
 du = [100/99 100/99 10/9 100/99];
 %max f(x) <=> min -f(x)
 f = @(x) -(c(1)*x(4)*x(7) - c(2)*x(1) - c(3)*x(2) - c(4)*x(3) - c(5)*x(5));
 %samma med bivillkor, h(x)>= 0 <=> -h(x) <= 0
 % h1(x)<= 0
 x = [1745, 12000, 110, 3048, 1974, 89.2, 92.8, 8, 3.6, 145]';
 x1=sym('x1'); x2=sym('x2'); x3=sym('x3'); x4=sym('x4'); x5=sym('x5');
 x6=sym('x6'); x7=sym('x7'); x8=sym('x8'); x9=sym('x9'); x10=sym('x10');
 %syms('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9' ,'x10');
 v = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10];
 
    function [constraint, grad, H] = symbolicManipulation(expression)
       symgrad = gradient(expression, v);
       temp = matlabFunction(symgrad, 'vars',v);
       grad = @(x) temp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10));
       temp = matlabFunction(expression,'vars', v);
       constraint = @(x) temp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10));
       hessian_symbolic = hessian(expression,v);
       temp = matlabFunction(hessian_symbolic, 'vars', v);
       H = @(x) temp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10));
       
    end
    all_constraints = [%Icke-linjära bivillkor <= 0
        -(x1*(1.12 +0.13167*x8 - 0.00667*x8^2) - dl(1)*x4);
        x1*(1.12 +0.13167*x8 - 0.00667*x8^2) - du(1)*x4;
    -(86.35 + 1.098*x8 - 0.038*x8^2 + 0.325*(x6-89) - dl(2)*x7);
    86.35 + 1.098*x8 - 0.038*x8^2 + 0.325*(x6-89) - du(2)*x7;
    -(35.82 - 0.222*x10 - dl(3)*x9);
    35.82 - 0.222*x10 - du(3)*x9;
    -(-133 + 3*x7 - dl(4)*x10);
    -133 + 3*x7 - du(4)*x10;
    1.22*x4 - x1 - x5; %likhetsvillkor index 9
    (98000*x3/(x4*x9+1000*x3)) - x6; %likhetsvillkor index 10
    ((x2+x5)/x1) - x8; %likhetsvillkor index 11
    %Bounds (<= 0) index 12 till 31
     x1-2000;
     -x1;
     x2-16000;
     -x2;
      x3-120;
      -x3;
      x4-5000;
      -x4;
      x5-2000;
      -x5;
      x6-93;
      85-x6;
      x7-95;
      90-x7;
      x8-12;
      3-x8;
      x9-4;
      1.2-x4;
      x10-162;
      145-x10];
 
     
     gradients = cell(1,length(all_constraints));
     constraints = cell(1,length(all_constraints)); 
     hessians = cell(1,length(all_constraints));
     for i=1:length(all_constraints)
        [con,grad, H] = symbolicManipulation(all_constraints(i));
        gradients{i} = grad;
        constraints{i} = con;
        hessians{i} = H;
     end
  

    function [c, ceq] = confun(x,d)
        %nonlinear inequality constraints
        %h(xk) + gradh(xk)'d <=0
        %index 1 till 8 och 12 till 31 
        
        %Funktioner.....
        olikhet1 = cellfun(@(f) f(x), constraints(1:8));
        olikhet2 = cellfun(@(f) f(x), constraints(12:end));
        
        likhet = cellfun(@(f) f(x), constraints(9:11));
        
        olikhet1grad = cellfun(@(f) f(x), gradients(1:8), 'UniformOutput', false);
        ineqgrads = cell2mat(olikhet1grad);
        olikhet2grad = cellfun(@(f) f(x), gradients(12:end), 'UniformOutput', false);
        ineqgrads2 = cell2mat(olikhet2grad);
        likhetgrad = cellfun(@(f) f(x), gradients(9:11), 'UniformOutput', false);
        eqgrad = cell2mat(likhetgrad);
        
        c = zeros(28,1);
        for i=1:28
           if i<= 8
               c(i) = olikhet1(i) + d'*ineqgrads(:,i);
           end
           if i>8
              c(i) = olikhet2(i-8) + d'*ineqgrads2(:,i-8);
           end
        end
        ceq = zeros(3,1);
        for i=9:11
    
           ceq(i-8) = likhet(i-8) + d'*eqgrad(:,i-8); 
        end
        
    end

gradf = @(x) -[-c(2),-c(3), -c(4), c(1)*x(7), -c(5), 0, c(1)*x(4), 0,0,0]';
hessianf = zeros(10);
hessianf(7,4) =c(1);
hessianf(4,7) = c(1);
Hf = @(x) hessianf;
%inequality constraints
    function [m] = Hh(x)
       m=1;
       for i=1:7
           temp = cell2mat(hessians(i));
           m = m*temp(x);
       end
       for i=12:31
           temp = cell2mat(hessians(i));
          m=m*temp(x);
       end
    end
%equality constraints
        function [val] = Hg(x)
            val=1;
            for i=1:31
                temp = cell2mat(hessians(i));
                val= val*temp(x);
            end
        end
    
    

B = 50*eye(10); %adderas till lagrange-hessianen 
lambda = -100*ones(10,1); my = 100*ones(10,1);
%grad2Lxx = @(x,lambda, my,B) Hf(x) + lambda'*Hh(x) + my'*Hg(x) + B; 
grad2Lxx = @(x, lambda, my, B) 50*eye(10);

x = [1745, 12000, 110, 3048, 1974, 89.2, 92.8, 8, 3.6, 145]';
d = ones(10,1);
q = @(d,x, lambda, my,B)  f(x) + gradf(x)'*d + 0.5*d'*grad2Lxx(x,lambda, my,B)*d;

   
while abs(f(x) - 1769) > 10
    A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
    d = fmincon(@(d) q(d,x, lambda, my,B),x,A,b,Aeq,beq,lb,ub,@(d) confun(x,d));
    xold = x;
    x = x+d;
  
    disp(f(x));
end
end

  