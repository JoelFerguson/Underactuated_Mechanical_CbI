%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Cart-pole CbI simulation example for "Total energy-shaping control for mechanical systems via Control-by-Interconnection"
% Authour: Joel Ferguson
% Date 22-Dec-2022
% Version: 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

% Useful function for indexing matrix elements of anonymous functions
mat_idx = @(A,i,j) A(i,j);

%% Define cart-pole system
% Define symbolic variables
syms q1 q2
syms p1 p2
q_sym = [q1; q2];
p_sym = [p1; p2];

% Input matrix
G = [1; 0];
Gp = [0 1];

% Mass matrix and components
mc=1; mp=1; lp=1; g = 9.8;
M = @(q) [mc+mp, mp*lp*cos(q(2));
        mp*lp*cos(q(2)), mp*lp^2];
Mi_int = matlabFunction(inv(M(q_sym)),'vars',q_sym);
Mi = @(q) Mi_int(q(1),q(2));
m11 = @(q) mat_idx(Mi(q),1,1);
m21 = @(q) mat_idx(Mi(q),2,1);
m22 = @(q) mat_idx(Mi(q),2,2);
% Gradient of mass matrix
T = @(q,p) 0.5*p.'*(M(q)\p);
dTdq_int = matlabFunction(simplify(jacobian(T(q_sym,p_sym),q_sym)).','vars',[q_sym; p_sym]);
dTdq = @(q,p) dTdq_int(q(1),q(2),p(1),p(2));
% Construction of 'E' matrix
E_int = matlabFunction(simplify(0.5*jacobian((M(q_sym)\p_sym),q_sym).'*M(q_sym)),'vars',[q_sym; p_sym]);
E = @(q,p) E_int(q(1),q(2),p(1),p(2));
% Check computations - expression elow should be equal to 0 if correct.
% dTdq(q_sym,p_sym) - E(q_sym,p_sym)*Mi(q_sym)*p_sym

% Potential energy
V = @(q) mp*g*lp*cos(q(2));
% Gradient of potential energy
dVdq_int = matlabFunction(jacobian(V(q_sym),q_sym).','vars',q_sym);
dVdq = @(q) dVdq_int(q(1),q(2));

% Hamiltonian
H = @(q,p) T(q,p) + V(q);

% Define open-loop dynamics --- x = [q; p]
dq = @(q,p) Mi(q)*p;
dp = @(q,p,u) - dTdq(q,p) - dVdq(q) + G*u;
dx = @(x,u) [dq(x(1:2),x(3:4)); dp(x(1:2),x(3:4),u)];

%% Define the added inverse mass at the origin
% Mai should be chosen to ensure that s1(), s3() > 0 at q2 = 0

% Choose the inverse added mass function ma11. This is a free function
ma11 = @(q) 0;
% Compute the gradient of ma11
dma11dq2_int = matlabFunction(diff(ma11(q_sym),q2),'vars',q_sym);
dma11dq2 = @(q) dma11dq2_int(q(1),q(2));

% Choose values for ma21, ma22 at the origin
ma21_0 = -2;
ma22_0 = 8;

% Define functions s1, s2, s3
s1_ma = @(q,ma21,ma22) (m22(q)+ma22) - ((m21(q)+ma21)/(m11(q)+ma11(q)))*(m21(q)+ma21);
s2_ma = @(q,ma21,ma22) ((m21(q)+ma21)/(m11(q)+ma11(q)))*m11(q) - m21(q);
s3_ma = @(q,ma21,ma22) ((m21(q)+ma21)/(m11(q)+ma11(q)))*m21(q) - m22(q);

% Verify that s1() > 0 at q = 0
q_0 = [0; 0];
s1_0 = s1_ma(q_0, ma21_0, ma22_0)

% Verify that s3() > 0 at q = 0
s3_0 = s3_ma(q_0, ma21_0, ma22_0)

% Verify that the minimum eigenvalue of Mi + Mai > 0. This is true if s1 >
% 0
Mdi_min_eig = min(eig(Mi(q_0) + [ma11(q_0) ma21_0; ma21_0 ma22_0]))


%% Solve kinetic energy PDE
% Construct symbolic expression for Mai and gradients
syms ma21 ma22
syms dma21dq2_sym dma22dq2_sym 
Mai_sym = @(q) [ma11(q) ma21; ma21 ma22];
dMaidq2_sym = @(q) [dma11dq2(q) dma21dq2_sym; dma21dq2_sym dma22dq2_sym];
dMaipdq_sym = @(q,p) [zeros(2,1) dMaidq2_sym(q)*p];

% Construct gradients of Mi
dMipdq_int = matlabFunction(jacobian(Mi(q_sym)*p_sym,q_sym),'vars',[q_sym; p_sym]);
dMipdq = @(q,p) dMipdq_int(q(1),q(2),p(1),p(2));

% Build Y matrix using symbolic definitions
Y_sym = @(q,p) 0.5*Mi(q)*dMaipdq_sym(q,p).' - 0.5*dMipdq(q,p)*Mai_sym(q);

% Define D matrix
D = @(q) [((m21(q)+ma21)/(m11(q)+ma11(q))) -1];

% Construct kinetic energy matching equation
constraint_full = D(q_sym)*(Y_sym(q_sym,p_sym) + Y_sym(q_sym,p_sym).')*D(q_sym).';
% Break matching equation into coefficients of p_i
constraint_component = jacobian(constraint_full,p_sym);
% Resolve matching equations into an ODE. Note that this could be solved
% implicitly rather then resolving explicitly
ODEs = solve([constraint_component(1); constraint_component(2)],[dma21dq2_sym; dma22dq2_sym]);
% Convert symbolic ODEs into functions
dma21dq2_int = matlabFunction(ODEs.dma21dq2_sym,'vars',[q2; ma21; ma22]);
dma22dq2_int = matlabFunction(ODEs.dma22dq2_sym,'vars',[q2; ma21; ma22]);
% Construct matching ODE with q2 as the independent variable and ma21, ma22
% as the states
dmadq2 = @(q2,ma) [dma21dq2_int(q2, ma(1), ma(2)); dma22dq2_int(q2, ma(1), ma(2))];

% Solve kinetic energy matching equations on the two half spaces [0,pi/2],
% [0,-pi/2]. In both cases, the initial conditions for ma21_0, ma22_0
% should be used.
options = odeset('RelTol',1e-6,'MaxStep',0.01);
ma0 = [ma21_0; ma22_0];
% Solve kinetic energy ODE on the half-space [0,pi/2]
q2_range = [0 pi/2];
[q2_num_1, ma_num_1] = ode23s(dmadq2, q2_range, ma0, options);
% Solve kinetic energy ODE on the half-space [0,-pi/2]
q2_range = [0 -pi/2];
[q2_num_2, ma_num_2] = ode23s(dmadq2, q2_range, ma0, options);

% Merge solutions from the two half-spaces. Note that the initial point is
% common between both datasets and should be removed from one.
q2_num_KE = [flip(q2_num_2); q2_num_1(2:end)];
ma21_num = [flip(ma_num_2(:,1)); ma_num_1(2:end,1)];
ma22_num = [flip(ma_num_2(:,2)); ma_num_1(2:end,2)];

% Generate functions for ma21, ma22 by interpolating numerical solutions
ma21 = @(q) interp1(q2_num_KE,ma21_num,q(2));
ma22 = @(q) interp1(q2_num_KE,ma22_num,q(2));

% Construct function for inverse added mass matrix
Mai = @(q) [ma11(q) ma21(q); ma21(q) ma22(q)];

% Create functions for gradients of Mai using the computed solutions
dMaidq2 = @(q) [dma11dq2(q) dma21dq2_int(q(2),ma21(q),ma22(q)); dma21dq2_int(q(2),ma21(q),ma22(q)) dma22dq2_int(q(2),ma21(q),ma22(q))];
dMaipdq = @(q,p) [zeros(2,1) dMaidq2(q)*p];

% Compute eigenvalues of Mi + Mai on domian to determine on which domain
% the closed-loop mass is positive
eig_vals = zeros(2,length(q2_num_KE));
ma11_val = zeros(1,length(q2_num_KE));
for i=1:length(q2_num_KE)
    ma11_val(i) = ma11([0;q2_num_KE(i)]);
    eig_vals(:,i) = eig(Mi([0;q2_num_KE(i)]) + [ma11([0;q2_num_KE(i)]), ma21([0;q2_num_KE(i)]); ma21([0;q2_num_KE(i)]), ma22([0;q2_num_KE(i)])]);
end

% Plot results
% Plot components of Mai
fig1 = figure(1);
subplot(1,2,1)
plot(q2_num_KE,ma11_val)
ylim([-2.5 8.5])
hold on
plot(q2_num_KE,ma21_num)
plot(q2_num_KE,ma22_num)
grid on
legend('m_{a11}','m_{a12}','m_{a22}')
xlabel('q_2')
% Plot closed-loop eigenvalues
subplot(1,2,2)
plot(q2_num_KE,min(eig_vals))
grid on
xlabel('q_2')
legend('\lambda_{min}[M^{-1}+M_a^{-1}]')
% Adjust line thicknesses
set(findall(fig1,'type','text'),'FontSize',11)
set(findall(fig1,'type','axes'),'FontSize',11)
set(findall(fig1,'type','line'),'linewidth',2.0)


%% Construct closed-loop potential energy function Vm
% Define functions s1, s2, s3 using solutions to kinetic energy PDE
s1 = @(q) s1_ma(q,ma21(q),ma22(q));
s2 = @(q) s2_ma(q,ma21(q),ma22(q));
s3 = @(q) s3_ma(q,ma21(q),ma22(q));

% Build potential energy ODE
dVmdq2 = @(q) -(s1(q)/s3(q))*Gp*dVdq(q);
dVmdt_wrap = @(q2,x) dVmdq2([0;q2]);
dVmdq = @(q) [0; dVmdq2(q)];

% Solve kinetic energy matching equations on the two half spaces [0,pi/2],
% [0,-pi/2]. In both cases, the initial conditions Vm(0) = 0 should be used
options = odeset('RelTol',1e-6,'MaxStep',0.01);
Vm0 = 0;
% Solve kinetic energy ODE on the half-space [0,pi/2]
q2_range = [0 pi/2];
[theta_Vm_1, Vm_num_1] = ode23s(dVmdt_wrap, q2_range, Vm0, options);
% Solve kinetic energy ODE on the half-space [0,-pi/2]
q2_range = [0 -pi/2];
[theta_Vm_2, Vm_num_2] = ode23s(dVmdt_wrap, q2_range, Vm0, options);

% Merge solutions from the two half-spaces. Note that the initial point is
% common between both datasets and should be removed from one.
q2_num_PE = [flip(theta_Vm_2); theta_Vm_1(2:end)];
Vm_num = [flip(Vm_num_2(:,1)); Vm_num_1(2:end)];

% Generate function for Vm by interpolating numerical solution
Vm = @(q) interp1(q2_num_PE,Vm_num,q(2));

% Plot results
fig2 = figure(2);
plot(q2_num_PE,Vm_num)
grid on
legend('V_m')
xlabel('q_2')
% Adjust line thicknesses
set(findall(fig2,'type','text'),'FontSize',11)
set(findall(fig2,'type','axes'),'FontSize',11)
set(findall(fig2,'type','line'),'linewidth',2.0)

%% Construct closed-loop potential energy function Vm
% Define beta functions for integrations
beta = @(q) G.'*(Mai(q) + Mi(q))*M(q);
beta1 = @(q) mat_idx(beta(q),1,1);
beta2 = @(q) mat_idx(beta(q),1,2);

% Define ODE to resolve for Gamma
dGammadq2 = @(q2) beta2([0;q2])/beta1([0;q2]);

% Solve Gamma ODE on the two half spaces [0,pi/2], [0,-pi/2]. In both cases, 
% the initial conditions Gamma(0) = 0 should be used
Gamma0 = 0;
options = odeset('RelTol',1e-6,'MaxStep',0.01);
% Solve kinetic energy ODE on the half-space [0,pi/2]
[gamma_num1,f_num1] = ode23s(@(q2,x) dGammadq2(q2), [0 q2_num_KE(end)], Gamma0, options);
% Solve kinetic energy ODE on the half-space [0,-pi/2]
[gamma_num2,f_num2] = ode23s(@(q2,x) dGammadq2(q2), [0 q2_num_KE(1)], Gamma0, options);

% Merge solutions from the two half-spaces. Note that the initial point is
% common between both datasets and should be removed from one.
q2_num_gamma = [flip(gamma_num2); gamma_num1(2:end)];
gamma_num = [flip(f_num2); f_num1(2:end)];

% Generate function for Gamma by interpolating numerical solution
Gamma = @(q) q(1) + interp1(q2_num_gamma,gamma_num,q(2));
dGammadq = @(q) [1, dGammadq2(q(2))];

% Define free potential function Vf in terms of argument Gamma
kappa = 5;
Vf = @(gamma) 0.5*kappa*gamma^2;
dVfdGamma = @(gamma) kappa*gamma;

% Define closed-loop potential energy and gradient
Vd = @(q) Vf(Gamma(q)) + Vm(q);
dVfdq = @(q) dGammadq(q).'*dVfdGamma(Gamma(q));
dVddq = @(q) dVmdq(q) + dVfdq(q);

% Plot results
% Define domain to evaluating Vd
s_num = [-0.5:0.01:0.5].';
% Evaluate Vd on specified domain
Vd_num = zeros(length(q2_num_KE), length(s_num));
for i=1:length(q2_num_KE)
    for j=1:length(s_num)
        Vd_num(i,j) = Vd([s_num(j) q2_num_KE(i)]);
    end
end
% Plot contour plot on log scale to show behaviour near zero
fig3 = figure(3);
contour(q2_num_KE, s_num, log(Vd_num.'), 'LineWidth',2)
xlabel('q_2')
ylabel('q_1')
legend('log(V_d)')
grid on
% Adjust line thicknesses
set(findall(fig3,'type','text'),'FontSize',11)
set(findall(fig3,'type','axes'),'FontSize',11)
set(findall(fig3,'type','line'),'linewidth',2.0)

%% Construct control law
% Define gradient of added energy function
dHadq = @(q,p) [0; 0.5*p.'*dMaidq2(q)*p] + dVddq(q);

% Defined closed-loop Hamiltonian function
Hd = @(q,p) 0.5*p.'*(Mi(q) + Mai(q))*p + Vd(q);

% Build Y matrix
Y = @(q,p) 0.5*Mi(q)*dMaipdq(q,p).' - 0.5*dMipdq(q,p)*Mai(q);

% Build J matrix. As J is 2x2, J11, J22 = 0.
J11 = 0;
J22 = 0;
Y11 = @(q,p) mat_idx(Y(q,p),1,1);
Y21 = @(q,p) mat_idx(Y(q,p),2,1);
J21 = @(q,p) ((ma21(q) + m21(q))/(ma11(q) + m11(q)))*(J11 - Y11(q,p)) + Y21(q,p);
J = @(q,p) [J11 -J21(q,p).'; J21(q,p) J22];

% Construct A, B, C
A = @(q,p) Mi(q)*p;
B = @(q,p) -E(q,p).'*Mai(q)*p + dHadq(q,p) - M(q)*J(q,p)*M(q)*A(q,p);
C = @(q,p) ((Mai(q) + Mi(q))\Mi(q))*B(q,p);

% Define damping injection term
v = @(q,p) -5*G.'*(Mai(q) + Mi(q))*p;

% Construct control signal
u = @(q,p) v(q,p) - G.'*(C(q,p) - dVdq(q));

%% Simulate system
% Construct ODE by substituting in the control signal
ode_wrap = @(t,x) dx(x,u(x(1:2),x(3:4)));

% Define simulation initial conditions
q0 = [0; 0.3];
p0 = [0; 0];
x0 = [q0; p0];

% Run simulation
t_sim = [0 5];
options = odeset('RelTol',1e-6);
[res.t,res.x] = ode23s(ode_wrap, t_sim, x0, options);

% Unpack results
res.q = res.x(:,1:2);
res.p = res.x(:,3:4);
% Compute closed-loop Hamiltonian
res.H = zeros(1,length(res.t));
res.Hd = zeros(1,length(res.t));
for i=1:length(res.t)
    res.H(i) = H(res.q(i,:).',res.p(i,:).');
    res.Hd(i) = Hd(res.q(i,:).',res.p(i,:).');
end

% Plot results
fig4 = figure(4);
% Plot state evolution
subplot(2,1,1)
plot(res.t,res.q)
hold on
plot(res.t,res.p)
legend('q_1','q_2','p_1','p_2')
grid on
xlabel('time (s)')
% Plot closed-loop Hamiltonian
subplot(2,1,2)
plot(res.t,res.Hd)
legend('H_d')
grid on
xlabel('time (s)')
% Adjust line thicknesses
set(findall(fig4,'type','text'),'FontSize',11)
set(findall(fig4,'type','axes'),'FontSize',11)
set(findall(fig4,'type','line'),'linewidth',2.0)
