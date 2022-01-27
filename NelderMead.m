
%optimum x = [1,1] f = 0
%f2 = @(x) 10^5*(x(2) -x(1))^2 + (1- x(1))^2
%f3 = @(x) -4*exp(-((x(1)+2)^2+ (x(2) +1)^2)/10) + 4*exp(-((x(1)+2)^2+ (x(2) +1)^2)/100) + ((x(1)+2)^2 + (x(2)+1)^2 +x(1))/100
function [ Point ] = NelderMead( f )
    

    x = [3,8,-1; 4,0,-5];
    
    function [ x ] = ReplaceVertex (x, x_max, x_new)
       for i=1:length(x(1,:))
          if x(:,i)==x_max
              x(:,i)=x_new;
          end
       end
    end
    function [ c ] = ExemptMax (x, x_max)
       for i=1:length(x(1,:))
           if x(:,i)==x_max
               c = num2cell([x(:,1:i-1),x(:,i+1:end)],1);
           end
       end
    end
   
  
   
    centroids = [];
    values = [];
    
    for i=1:100
        c = num2cell(x,1);
        cellfun(f,c);
        [~, idx_min] = min(cellfun(f,c));
        x_min= cell2mat(c(idx_min));
        [~, idx_max] = max(cellfun(f,c));
        x_max = cell2mat(c(idx_max));
        x_hat = (-x_max + sum(x,2))/length(x(:,1));
        x_ref = 2*x_hat - x_max;
        
        c = ExemptMax(x,x_max);
        centroids = [centroids x_hat];
        values = [values f(x_hat)];
        
        
        %Fall (1), steg 2
        if f(x_min) > f(x_ref)
            x_exp = 2*x_ref - x_hat;
            x_new = (f(x_exp)<f(x_ref))*x_exp + (f(x_exp)>=f(x_ref))*x_ref;
            x = ReplaceVertex(x, x_max, x_new);
            continue;
        end
        
        %Fall (2), steg 3
        
        if max(cellfun(f,c))> f(x_ref) && f(x_ref) >= f(x_min)
            x_new = x_ref;
            x = ReplaceVertex(x, x_max, x_new);
        end
        
        %Fall (3), steg 4
        if f(x_ref)>= max(cellfun(f, c))
            x_new = (f(x_max)<=f(x_ref))*(x_max+x_hat)/2 + (f(x_max)>f(x_ref))*(x_ref + x_hat)/2;
            x= ReplaceVertex(x, x_max, x_new);
        end
        
        
    end
    plot(1:length(values), values);
    title('Funktionsvärde som funktion av iterationerna');
    figure;
    scatter(centroids(1,:), centroids(2,:)); 
    title('Centroidernas x1 resp x2 koordinater');
    disp(centroids);
    Point = x_min; 
end
