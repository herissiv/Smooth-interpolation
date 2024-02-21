R = @(x) 1/(1+25*x^2);
I = [-1,1];
ant_x = 9;
avstand = abs(I(length(I))-I(1)); 

x_er = [];
y_er = [];

for i = 1:(ant_x)
    x_er = [x_er,I(1)+ (avstand/(ant_x-1))*(i-1)]; 
    y_er = [y_er,R(x_er(i))];
end

hold on %Use this so that every "plot" gets plotted in the same plot
punkter(x_er,y_er)
[P8,EL] = Lagrange8(x_er,y_er,R)
[EL] = CuSplines(x_er,y_er,R)
legend('Punkter','Lagrange','Cubic splines')


function punkter(x_er,y_er)
    plot(x_er,y_er,'*')
end 


function [P8,EL] = Lagrange8(x_er,y_er,R)
    P8 = @(x) 0; %Set to zero, because I'm going to add to it 
    x_ledd = @(x) 1; %Set to one, because I'm going to multiply to it  
    for i = 1:length(x_er)
        if i ~= 1 %Divide it into 3 groups because the two next following 
            %"ranges" will fail at the endpoints. 
            for j = 1:(i-1) %Here i = 1 would have failed 
                x_ledd = @(x) x_ledd(x) * (x-x_er(j))/(x_er(i)-x_er(j));
                %Here we will get every x up to x_1, 
                %without x_i, from x_er in the numerator.
            end 
            for j = (i+1):length(x_er) %Here i = length(x_er) would have 
                %failed 
                x_ledd = @(x) x_ledd(x) * (x-x_er(j))/(x_er(i)-x_er(j));
                %Here we will get every x after x_i, 
                %without x_i, from x_er in the numerator. 
            end 
        elseif i == 1 %Special case
            for j = (i+1):length(x_er)
                x_ledd = @(x) x_ledd(x) * (x-x_er(j))/(x_er(i)-x_er(j));
            end
        elseif i == length(x_er) %Special case 
            for j = 1:(i-1)
                x_ledd = @(x) x_ledd(x) * (x-x_er(j))/(x_er(i)-x_er(j));
            end    
        end 
        P8 = @(x) P8(x) + y_er(i) * x_ledd(x); %Update P8 for every i
        x_ledd = @(x) 1; %Reset for every iteration 
    end 
    fplot(P8,[-1,1],'yellow')
    EL = []; %Errorlist  
    for i = 1:length(x_er)-1
        xe = (x_er(i+1)+x_er(i))/2; %Middle point between points 
        EL = [EL, abs(R(xe)-P8(xe))]; %List with errors for each "middlepoint" 
    end 
end 

function EL = CuSplines(x_er,y_er,R)
    %Information for I = [x_i,x_(i+1)]:
        %S_i(x_(i+1)) = S_(i+1)(x_(i+1))
        %S_i(x_i) = y_i
        %S_i'(x_(i+1)) = S_(i+1)'(x_(i+1))
        %S_i''(x_(i+1)) = S_(i+1)''(x_(i+1))
        %S_i(x) = y_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3
        %Number of S is equal to number of point minus 1.
    %Following function will take the following points to consideration 
    ant = length(x_er) - 1;
    dx = [];
    dy = [];
    n = length(x_er);
    A = zeros(n,n); %Take base in the equation A*c = r (3.24, page 170)
    r = zeros(n,1); % => c = A^(-1)*r given that A is invertible, as I assume 
    A(1,1) = 1;
    A(n,n) = 1;
    r(1) = 0; %Know that c_1 = 0
    r(n) = 0; %Know that c_n = 0
    
    for i = 1:ant %loop under use "i+1"
        dx = [dx, x_er(i+1) - x_er(i)];
        dy = [dy, y_er(i+1) - y_er(i)];
    end 
    
    for i = 2:ant %Already set the first and last row in A, to ant because dx use "i+1"
        A(i,i-1) = dx(i-1);
        A(i,i) = 2*dx(i-1)+2*dx(i);
        A(i,i+1) = dx(i);
        r(i) = 3*(dy(i)/dx(i) - dy(i-1)/dx(i-1));
    end 
    
    c = mtimes(inv(A),r); % c = A^(-1)*r
    d = []; %List with the d constants 
    b = []; %List with the b constants 
    
    for i = 1:ant %we use i+1 in the loop under 
        d = [d,(c(i+1)-c(i))/(3*dx(i))]; %Equation 3.22, page 170
        b = [b,dy(i)/dx(i) - c(i)*dx(i)-d(i)*(dx(i))^2]; %Equation 3.23, page 170
    end 
    EL = []; 
    for i = 1:ant
        S{i} = @(x) y_er(i) + b(i)*(x-x_er(i))+c(i)*(x-x_er(i))^2 + d(i)*(x-x_er(i))^3; %Equation 3.17, page 167
        fplot(S{i},[x_er(i),x_er(i+1)],'green') %Plot every function given their own interval [x_i,x_(i+1)]
        xe = (x_er(i+1)+x_er(i))/2; 
        EL = [EL, abs(R(xe) - S{i}(xe))];
    end  
end 


