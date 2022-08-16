%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Acrobot CbI simulation example for "Total energy-shaping control for mechanical systems via Control-by-Interconnection"
% Authour: Joel Ferguson
% Date 12-August-2022
% Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

% Useful function for indexing matrix elements of anonymous functions
mat_idx = @(A,i,j) A(i,j);

%% Define acrobot system
% Define symbolic variables
syms q2 q1
syms p2 p1
q_sym = [q1; q2];
p_sym = [p1; p2];

% Input matrix
G = [1; 0];
Gp = [0 1];

% Mass matrix and components
m1 = 1; Jl1 = 1.3333; l1 = 4; lc1 = 2;
m2 = 4; Jl2 = 0.3333; l2 = 1; lc2 = 0.5;
c1 = m2*lc2^2 + m1*l2^2 + Jl2;
c2 = m1*lc1^2 + Jl1;
c3 = m1*l2*lc1;
c4 = m2*lc2 + m1*l2;
c5 = m1*lc1;
g = 9.8;
M = @(q) [c2, c2+c3*cos(q(1));
        c2+c3*cos(q(1)), c1+c2+2*c3*cos(q(1))];
Mi_int = matlabFunction(inv(M(q_sym)),'vars',q_sym);
Mi = @(q) Mi_int(q(1),q(2));
m11 = @(q) mat_idx(Mi(q),1,1);
m12 = @(q) mat_idx(Mi(q),1,2);
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
V = @(q) g*(c4*cos(q(2)) + c5*cos(q(1)+q(2)));
% Gradient of potential energy
dVdq_int = matlabFunction(jacobian(V(q_sym),q_sym).','vars',q_sym);
dVdq = @(q) dVdq_int(q(1),q(2));

% Hamiltonian
H = @(q,p) T(q,p) + V(q);

% Define open-loop dynamics --- x = [q; p]
dq = @(q,p) Mi(q)*p;
dp = @(q,p,u) - dTdq(q,p) - dVdq(q) + G*u;
dx = @(x,u) [dq(x(1:2),x(3:4)); dp(x(1:2),x(3:4),u)];

%% Solve kinetic energy PDE
% Using the reported solution from Mahindrakar et al.
Mdi = [0.3385 -0.9997; -0.9997 5.9058];

% Define values for ma12, ma22 at the origin from reported solution
q0 = [0;0];
ma12_0 = Mdi(1,2) - mat_idx(Mi(q0),1,2);
ma22_0 = Mdi(2,2) - mat_idx(Mi(q0),2,2);

% Define expression for ma11 using reported solution
ma11 = @(q) Mdi(1,1) - mat_idx(Mi(q),1,1);
% Compute the gradient of ma11
dma11dq1_int = matlabFunction(diff(ma11(q_sym),q1),'vars',q_sym);
dma11dq1 = @(q) dma11dq1_int(q(1),q(2));

% Construct symbolic expression for Mai and gradients
syms ma12 ma22
syms dma12dq1_sym dma22dq1_sym 
Mai_sym = @(q) [ma11(q) ma12; ma12 ma22];
dMaidq1_sym = @(q) [dma11dq1(q) dma12dq1_sym; dma12dq1_sym dma22dq1_sym];
dMaipdq_sym = @(q,p) [dMaidq1_sym(q)*p zeros(2,1)];

% Construct gradients of Mi
dMipdq_int = matlabFunction(jacobian(Mi(q_sym)*p_sym,q_sym),'vars',[q_sym; p_sym]);
dMipdq = @(q,p) dMipdq_int(q(1),q(2),p(1),p(2));

% Build Y matrix using symbolic definitions
Y_sym = @(q,p) 0.5*Mi(q)*dMaipdq_sym(q,p).' - 0.5*dMipdq(q,p)*Mai_sym(q);

% Define D matrix
D = @(q) [((m12(q)+ma12)/(m11(q)+ma11(q))) -1];

% Construct kinetic energy matching equation
constraint_full = D(q_sym)*(Y_sym(q_sym,p_sym) + Y_sym(q_sym,p_sym).')*D(q_sym).';
% Break matching equation into coefficients of p_i
constraint_component = jacobian(constraint_full,p_sym);
% Resolve matching equations into an ODE. Note that this could be solved
% implicitly rather then resolving explicitly
ODEs = solve([constraint_component(1); constraint_component(2)],[dma12dq1_sym; dma22dq1_sym]);
% Convert symbolic ODEs into functions
dma12dq1_int = matlabFunction(ODEs.dma12dq1_sym,'vars',[q1; ma12; ma22]);
dma22dq1_int = matlabFunction(ODEs.dma22dq1_sym,'vars',[q1; ma12; ma22]);
% Construct matching ODE with q2 as the independent variable and ma12, ma22
% as the states
dmadq1 = @(q1,ma) [dma12dq1_int(q1, ma(1), ma(2)); dma22dq1_int(q1, ma(1), ma(2))];

% Solve kinetic energy matching equations on the two half spaces [0,pi],
% [0,-pi]. In both cases, the initial conditions for ma12_0, ma22_0
% should be used.
options = odeset('RelTol',1e-6,'MaxStep',0.01);
ma0 = [ma12_0; ma22_0];
% Solve kinetic energy ODE on the half-space [0,pi]
q1_range = [0 pi];
[q1_num_1, ma_num_1] = ode23s(dmadq1, q1_range, ma0, options);
% Solve kinetic energy ODE on the half-space [0,-pi]
q1_range = [0 -pi];
[q1_num_2, ma_num_2] = ode23s(dmadq1, q1_range, ma0, options);

% Merge solutions from the two half-spaces. Note that the initial point is
% common between both datasets and should be removed from one.
q1_num_KE = [flip(q1_num_2); q1_num_1(2:end)];
ma12_num = [flip(ma_num_2(:,1)); ma_num_1(2:end,1)];
ma22_num = [flip(ma_num_2(:,2)); ma_num_1(2:end,2)];

% Generate functinos for ma12, ma22 by interpolating numerical solutions
ma12 = @(q) interp1(q1_num_KE,ma12_num,q(1));
ma22 = @(q) interp1(q1_num_KE,ma22_num,q(1));

% Construct function for inverse added mass matrix
Mai = @(q) [ma11(q) ma12(q); ma12(q) ma22(q)];

% Create functions for gradients of Mai using the computed solutions
dMaidq1 = @(q) [dma11dq1(q) dma12dq1_int(q(1),ma12(q),ma22(q)); dma12dq1_int(q(1),ma12(q),ma22(q)) dma22dq1_int(q(1),ma12(q),ma22(q))];
dMaipdq = @(q,p) [dMaidq1(q)*p zeros(2,1)];

% Compute eigenvalues of Mi + Mai on domian to determine on which domain
% the closed-loop mass is positive
eig_vals = zeros(2,length(q1_num_KE));
ma11_val = zeros(1,length(q1_num_KE));
for i=1:length(q1_num_KE)
    ma11_val(i) = ma11([q1_num_KE(i);0]);
    eig_vals(:,i) = eig(Mi([q1_num_KE(i);0]) + [ma11([q1_num_KE(i);0]), ma12([q1_num_KE(i);0]); ma12([q1_num_KE(i);0]), ma22([q1_num_KE(i);0])]);
end

% Plot results
% Plot components of Mai
fig1 = figure(1);
subplot(1,2,1)
plot(q1_num_KE,ma11_val)
ylim([-2.5 8.5])
hold on
plot(q1_num_KE,ma12_num)
plot(q1_num_KE,ma22_num)
grid on
legend('m_{a11}','m_{a12}','m_{a22}')
xlabel('q_1')
% Plot closed-loop eigenvalues
subplot(1,2,2)
plot(q1_num_KE,min(eig_vals))
ylim([0 0.5])
grid on
xlabel('q_1')
legend('\lambda_{min}[M^{-1}+M_a^{-1}]')
% Adjust line thicknesses
set(findall(fig1,'type','text'),'FontSize',11)
set(findall(fig1,'type','axes'),'FontSize',11)
set(findall(fig1,'type','line'),'linewidth',2.0)

%% Construct closed-loop potential energy function Vm
% Define functions s1, s2, s3 using solutions to kinetic energy PDE
s1 = @(q) (m22(q)+ma22(q)) - ((m12(q)+ma12(q))/(m11(q)+ma11(q)))*(m12(q)+ma12(q));
s2 = @(q) ((m12(q)+ma12(q))/(m11(q)+ma11(q)))*m11(q) - m12(q);
s3 = @(q) ((m12(q)+ma12(q))/(m11(q)+ma11(q)))*m12(q) - m22(q);

% Define ODE describing the functinos f1, f2
dfdq1 = @(q1,f) s2([q1;0])\[g*c4*s1([q1;0]) + g*c5*cos(q1)*s1([q1;0]) + s3([q1;0])*f(2);
                g*c5*sin(q1)*s1([q1;0]) - s3([q1;0])*f(1)];

% Define value of functinos f1, f2 at the origin
f1_0 = 0;
f2_0 = -50;

% Solve for f1, f2 on the two half spaces [0,pi],
% [0,-pi]. In both cases, the initial conditions for f1_0, f2_0
% should be used.
f_0 = [f1_0; f2_0];
options = odeset('RelTol',1e-6,'MaxStep',0.005);
% Solve ODE on the half-space [0,pi]
q1_range = [0 pi];
[q1_fNum_1, f_num_1] = ode23s(dfdq1, q1_range, f_0, options);
% Solve ODE in second half-space  [0,-pi]
q1_range = [0 -pi];
[q1_fNum_2, f_num_2] = ode23s(dfdq1, q1_range, f_0, options);

% Merge solutions from the two half-spaces. Note that the initial point is
% common between both datasets and should be removed from one.
q1_num_f = [flip(q1_fNum_2); q1_fNum_1(2:end)];
f1_num = [flip(f_num_2(:,1)); f_num_1(2:end,1)];
f2_num = [flip(f_num_2(:,2)); f_num_1(2:end,2)];

% Generate functions for f1, f2 and gradients by interpolating numerical solutions
f1 = @(q) interp1(q1_num_f,f1_num,q(1));
f2 = @(q) interp1(q1_num_f,f2_num,q(1));
df1dq1 = @(q) [1 0]*dfdq1(q(1),[f1(q); f2(q)]);
df2dq1 = @(q) [0 1]*dfdq1(q(1),[f1(q); f2(q)]);

% Construct functinos for Vm and gradient
Vm = @(q) f1(q)*sin(q(2)) + f2(q)*cos(q(2));
dVmdq = @(q) [df1dq1(q)*sin(q(2)) + df2dq1(q)*cos(q(2));
                f1(q)*cos(q(2)) - f2(q)*sin(q(2))];

% Plot results
% Define domain for plotting
q1_num = [-pi:0.2:pi];
q2_num = [-pi:0.2:pi];
% Evaluate Vm on domain
Vm_num_2d = zeros(length(q1_num),length(q2_num));
for i=1:length(q1_num)
    for j=1:length(q2_num)
        Vm_num_2d(i,j) = Vm([q1_num(i); q2_num(j)]);
        qi = [q1_num(i); q2_num(j)];
    end
end
fig2 = figure(2);
surf(q1_num,q2_num,Vm_num_2d.')
grid on
xlabel('q_1')
ylabel('q_2')
zlabel('V_m')
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
dGammadq1 = @(q1) beta1([q1;0])/beta2([q1;0]);

% Solve Gamma ODE on the two half spaces [0,pi], [0,-pi]. In both cases, 
% the initial conditions Gamma(0) = 0 should be used
Gamma0 = 0;
options = odeset('RelTol',1e-6,'MaxStep',0.01);
% Solve kinetic energy ODE on the half-space [0,pi/2]
[gamma_num1,f_num1] = ode23s(@(q1,x) dGammadq1(q1), [0 q1_num_KE(end)], Gamma0, options);
% Solve kinetic energy ODE on the half-space [0,-pi/2]
[gamma_num2,f_num2] = ode23s(@(q1,x) dGammadq1(q1), [0 q1_num_KE(1)], Gamma0, options);

% Merge solutions from the two half-spaces. Note that the initial point is
% common between both datasets and should be removed from one.
q1_num_gamma = [flip(gamma_num2); gamma_num1(2:end)];
gamma_num = [flip(f_num2); f_num1(2:end)];

% Generate function for Gamma by interpolating numerical solution
Gamma = @(q) q(2) + interp1(q1_num_gamma,gamma_num,q(1));
dGammadq = @(q) [dGammadq1(q(2)), 1];

% Define free potential function Vf in terms of argument Gamma
kappa = 250;
Vf = @(gamma) 0.5*kappa*gamma^2;
dVfdGamma = @(gamma) kappa*gamma;

% Define closed-loop potential energy and gradient
Vd = @(q) Vf(Gamma(q)) + Vm(q);
dVfdq = @(q) dGammadq(q).'*dVfdGamma(Gamma(q));
dVddq = @(q) dVmdq(q) + dVfdq(q);

% Plot results
% Define domain for plotting
q1_num = [-pi:0.05:pi];
q2_num = [-pi:0.05:pi];
% Evaluate Vm on domain
Vd_num = zeros(length(q1_num),length(q2_num));
for i=1:length(q1_num)
    for j=1:length(q2_num)
        Vd_num(i,j) = Vd([q1_num(i) q2_num(j)]);
    end
end
% Plot contour plot on log scale to show behaviour near zero. Note that the
% minimum has been shifted for ease of viewing
fig3 = figure(3);
% surf(q1_num, q2_num, Vd_num.', 'LineWidth',2)
contour(q1_num, q2_num, log(Vd_num - min(min(Vd_num)) + 0.5).', 'LineWidth',2)
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
dHadq = @(q,p) [0.5*p.'*dMaidq1(q)*p; 0] + dVddq(q);

% Defined closed-loop Hamiltonian function
Hd = @(q,p) 0.5*p.'*(Mi(q) + Mai(q))*p + Vd(q);

% Build Y matrix
Y = @(q,p) 0.5*Mi(q)*dMaipdq(q,p).' - 0.5*dMipdq(q,p)*Mai(q);

% Build J matrix. As J is 2x2, J11, J22 = 0.
J11 = 0;
J22 = 0;
Y11 = @(q,p) mat_idx(Y(q,p),1,1);
Y12 = @(q,p) mat_idx(Y(q,p),1,2);
J12 = @(q,p) ((ma12(q) + m12(q))/(ma11(q) + m11(q)))*(J11 - Y11(q,p)) + Y12(q,p);
J = @(q,p) [J11 -J12(q,p).'; J12(q,p) J22];

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
q0 = [0; 0.5];
p0 = [0; 0];
x0 = [q0; p0];

% Run simulation
t_sim = [0 20];
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
% Plot configuration evolution
subplot(3,1,1)
plot(res.t,res.q)
legend('q_1','q_2')
grid on
xlabel('time (s)')
% Plot momentum evolution
subplot(3,1,2)
plot(res.t,res.p)
legend('p_1','p_2')
grid on
xlabel('time (s)')
% Plot closed-loop Hamiltonian
subplot(3,1,3)
plot(res.t,res.Hd - min(res.Hd))
legend('H_d')
grid on
xlabel('time (s)')
% Adjust line thicknesses
set(findall(fig4,'type','text'),'FontSize',11)
set(findall(fig4,'type','axes'),'FontSize',11)
set(findall(fig4,'type','line'),'linewidth',2.0)

