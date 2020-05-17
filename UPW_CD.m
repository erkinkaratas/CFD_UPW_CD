
%November 11,2019 
%__________________________________________________________________________
clc;
clearvars;
close all;

%VARIABLES 
n=10;           % # of nodes
rho = 2;        % Density     
gamma = 0.5;    
L = 2;          %Length
u = 4;          %Velocity
delx = L/n;     %Distance between each node
phi_A = 0;      %Left boundary Condition
phi_B = 0;      %Right Boundary Condition

%Convective Terms
F_w = rho*u;    %West Face
F_e = F_w;      %East Face
F_A = rho * u;  %Left Boundary
F_B = rho * u;  %Right Boundary

%Diffusive Terms
D_e = gamma/delx;       %East Face
D_w = D_e;              %West Face
D_A = gamma/(delx/2);   %Left Boundary
D_B = D_A;              %Right Boundary
    

S = (rho * u / L)*delx;   % Source Term
B = S * ones(n,1);        % Right hand side of the matrix equation

phi_CD = zeros(n);          %Generating empty matrix for CD Scheme
phi_UP = zeros(n);          %Generating empty matrix for Upwind Scheme
T_analy_sol = zeros(1,n);   %Generating empty matrix for Analytical Solution
%__________________________________________________________________________


%******************** CENTRAL DIFFERENCING SCHEME *************************
            
for i = 2:length(phi_CD)-1 %Internal Nodes
            %Plugging the internal nodes into a matrix
            
    phi_CD(i,i-1:i+1) = [-D_e - F_e/2; D_A; -D_w + F_w/2];
end
phi_CD(1,1) = phi_CD(1,1) + F_e/2 + D_A + D_e;    %first row, first column
phi_CD(1,2) = phi_CD(1,2) + F_e/2 - D_e;          %first row, second column

phi_CD(n,n) = phi_CD(n,n) + D_B + D_w - F_w/2;     %last row, last column
phi_CD(n,n-1) = phi_CD(n,n-1) - F_w/2 - D_w;       %last row, (n-1)th column


Central_Diff_Result = phi_CD\B;      %Inverse multiplication
                                        %FINAL RESULT
%__________________________________________________________________________



%************************* UPWIND METHOD SCHEME ***************************

for i = 2:length(phi_UP)-1 %Internal Nodes
            %Plugging the internal nodes into a matrix

    phi_UP(i,i-1:i+1) = [-F_e-D_e; D_w + D_e + F_w; -D_e];
    
end

phi_UP(1,1) = phi_UP(1,1) + F_e + D_e + D_w; %first row, first column
phi_UP(1,2) = phi_UP(1,2) - D_A;             %first row, second column

phi_UP(n,n-1) = phi_UP(n,n-1) - (F_w + D_e); %last row, (n-1)th column
phi_UP(n,n) = phi_UP(n,n) + D_w + D_e + F_B; %last row, last column

Upwind_Result = phi_UP\B;           %Inverse multiplication
                                        %FINAL RESULT
%__________________________________________________________________________



%******************** ANALYTICAL SOLUTION *********************************
T_analy_sol(1) = ((delx/2)/L) - ((1 - exp(rho*u*delx/gamma)) / (1 - exp(rho*u*L/gamma)));

x = (3 * delx/2) : delx : (L - 3*delx/2);

for m = 2:n-1
    T_analy_sol(m) = (x(m-1)/L) - ((1 - exp(rho*u*x(m-1)/gamma)) / (1 - exp(rho*u*L/gamma)));
end
T_analy_sol(n) = ((L-(delx/2))/L) - ((1 - exp(rho*u*(L-delx/2)/gamma)) / (1 - exp(rho*u*L/gamma)));
%T_analy_sol = T_analy_sol*n;
T_analy_sol = T_analy_sol';     %FINAL RESULT
%__________________________________________________________________________


%******************** ERROR CALCULATION ***********************************
%Central Differencing Erreor %
error_Central_Dif = abs((Central_Diff_Result - T_analy_sol) / T_analy_sol) * 100;
error_Central_Dif = nonzeros(error_Central_Dif); %Removing zero terms from the matrix

%Upwind Scheme Error %
error_Upwind = abs((Upwind_Result - T_analy_sol) / T_analy_sol) * 100;
error_Upwind = nonzeros(error_Upwind);          %Removing zero terms from the matrix

%__________________________________________________________________________


%***********************Comparing the Results by Plotting******************
figure(1)
plot(T_analy_sol, 'g',"LineWidth", 3)
%delx/2:delx:L, 
hold on
grid on
plot(Upwind_Result, 'r x',"LineWidth", 2)
%delx/2:delx:L,
plot(Central_Diff_Result, 'b o',"LineWidth", 2)
xlabel('Node Numbers')
ylabel('Phi Values')
legend('Analytical Solution', 'Upwind Scheme', 'Central Differencing')
%__________________________________________________________________________


%************************* Error Plotting *********************************
figure(2)
plot(error_Central_Dif, 'c', "LineWidth", 3 )
hold on
grid on
plot(error_Upwind, 'm', "LineWidth", 3)
xlabel('Node Number')
ylabel('Percentage Error')
legend('Central Differencing Error', 'Upwind Scheme Error')
%__________________________________________________________________________


%************************** TABLE *****************************************
Node_Number = (1:n)';
Analytical_Results = T_analy_sol;
Error_CD = error_Central_Dif;
Error_UPW = error_Upwind;
Table1 = table(Node_Number,Analytical_Results,Central_Diff_Result,Error_CD,Upwind_Result,Error_UPW);
disp(Table1);
%__________________________________________________________________________
                   

