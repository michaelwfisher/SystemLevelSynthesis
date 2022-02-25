
% This is a numerical example motivated by power converter control
% design, and solved using system level synthesis with the finite
% impulse response (FIR), which is referred to here as deadbeat
% control (DBC), and by the simple pole approximation (SPA).

%% Setup

% simple pole approximation (SPA) approximate number of controller poles
% number of poles specified 
Ps = [6 14];

% deadbeat control (DBC) lengths of finite impule response (FIR)
Ts = [31 50 100];
% note: DBC is solved using golden section search for each T in Ts

% sampling time
h = 0.001;

% length of impulse and step responses
Tlength = 10000;

% state space model
A =  [0.988 0 0 0 0;
      0 0 0 0 0;
      1 0 0 0 0;
      0 0 0 0.995 0;
      0 0 0 0 0.9];

B =  [0 0;
      0 0;
      0 0;
      0.005 0;
      0 0.1];

Bh = [-0.0001 0;
       0 1;
       0.0066 0;
       0 0;
       0 0];

C =  [0 0.829 -0.428 1.02 0;
      0 0.428 0.829 0 -1.02];

n = size(A,2);
m = size(B,2);

% weight matrices in objective
Q = eye(n);
QQ = eye(2);
R = 0.01*eye(m);

% impulse response of desired transfer function
tauDes = 1;
tauqDes = 1;

dDes = 0.053;
dqDes = 100;

Ades = [(1-(h/tauDes)) (1-(h/tauqDes))]';

z = tf('z',h);
Tdes_tf = [-((h*dDes)/tauDes)/(z-(1-(h/tauDes))) 0;
           0 -((h*dqDes)/tauqDes)/(z-(1-(h/tauqDes)))];

T = Tlength;
impDesP = zeros(1,T);
for index = 1:T
    impDesP(index) = (-(h*dDes)/tauDes)*((1-(h/tauDes))^(index-1));
end
impDesQ = zeros(1,T);
for index = 1:T
    impDesQ(index) = (-(h*dqDes)/tauqDes)*((1-(h/tauqDes))^(index-1));
end
gDesImp = zeros(2,2*T);
for index=1:T
    gDesImp(1,2*(index-1)+1) = impDesP(index);
    gDesImp(2,2*(index-1)+2) = impDesQ(index);
end
Tdes_imp = gDesImp;

%% Finding True Optimal Solution

Aaug = [A zeros(5,2);
        zeros(2,5) diag([1-h/tauDes 1-h/tauqDes])];

Baug = [Bh B;
        diag([-h*dDes/tauDes -h*dqDes/tauqDes]) zeros(2,2)];

Caug = [QQ*([C zeros(2,2)] - [zeros(2,5) eye(2)]);
        R*zeros(2,7);
        eye(5) zeros(5,2)];

Daug = [zeros(2,2) zeros(2,2);
        zeros(2,2) R*eye(2);
        zeros(5,2) zeros(5,2)];

Plant = ss(Aaug,Baug,Caug,Daug,h);
[KK,TT,gg] = h2syn(Plant,5,2);

TTdes = ss(diag([1-h/tauDes 1-h/tauqDes]),...
           diag([-h*dDes/tauDes -h*dqDes/tauqDes]),...
           [QQ*eye(2);
            zeros(2,2)],...
           zeros(4,2),...
           h);

Topt = TT+TTdes;

%% Deadbeat Control Method (i.e. Finite Impulse Response)


% initialize structs for storing output of DBC design
Phix_temps = cell(1,length(Ts));
Phiu_temps = cell(1,length(Ts));    
    
Phixs = cell(1,length(Ts));
Phius = cell(1,length(Ts));

sls_inverses = cell(1,length(Ts));
sls_phixs = cell(1,length(Ts));
sls_phius = cell(1,length(Ts));

sls_objs = cell(1,length(Ts));
sls_gammas = cell(1,length(Ts));
sls_iterations = cell(1,length(Ts));

sls_costs = zeros(1,length(Ts));

    
for ii=1:length(Ts)
    
T = Ts(ii);    

% initialize variables
Phi_x = sdpvar(n,n*T,'full');
Phi_u = sdpvar(m,n*T,'full');

% parameters for setup
converged = 0;
lower = 0;
upper = 1;
convTol = 0.01;
iteration = 0;

% golden ratio
phig = (1+sqrt(5))/2;

midlow = upper - (upper-lower)/phig;
midup = lower + (upper-lower)/phig;

% run golden section search
while ~converged
    
for jj=1:2
    
if jj==1
    gamma = midlow;
elseif jj==2
    gamma = midup;
end

disp('Making objective');
tic;

Mobj = [QQ*(C*Phi_x*kron(eye(T),Bh) - Tdes_imp(:,1:m*T));
        R*Phi_u*kron(eye(T),Bh)];
Objective = norm(Mobj,'fro');    

disp('Objective made');
toc;

disp('Making constraints');
tic;

Constraints = [Phi_x(:,1:n) == eye(n)];

for i=1:(T-1)
    Constraints = [Constraints, ...
        Phi_x(:,(i*n+1):((i+1)*n)) - A*Phi_x(:,((i-1)*n+1):(i*n)) ...
        - B*Phi_u(:,((i-1)*n+1):(i*n)) == 0];
end

% slack variable constraint for DBC:
V = -A*Phi_x(:,((T-1)*n+1):(T*n)) -B*Phi_u(:,((T-1)*n+1):(T*n));
Constraints = [Constraints,norm(V,2) <= gamma];


disp('Constraints made');
toc;

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

% Solve the problem
sol = optimize(Constraints,Objective,options);

obj = (1/(1-gamma))*value(Objective);

if sol.problem == 0
    if jj==1
        lowerSol = obj;
        lowerSolved = 1;
    elseif jj==2
        upperSol = obj;
        upperSolved = 1;
    end
else
    if jj==1
        lowerSol = 1e10;
        lowerSolved = 0;
        obj = lowerSol;
    elseif jj==2
        upperSol = 1e10;
        upperSolved = 0;
        obj = upperSol;
    end
end

end

iteration
lower
upper
lowerSol
upperSol


if lowerSol < upperSol
    upper = midup;
else
    lower = midlow;
end

midlow = upper - (upper-lower)/phig;
midup = lower + (upper-lower)/phig;

lower
upper
midlow
midup

if ~upperSolved && ~lowerSolved && (upper-lower < convTol)
    converged = 1;
elseif upperSolved && lowerSolved && abs(upperSol-lowerSol)/...
        ((upperSol+lowerSol)/2) < convTol
    converged = 1;
else
    iteration = iteration+1;
end

end

sls_objs{ii} = obj;
sls_gammas{ii} = gamma;
sls_iterations{ii} = iteration;

x = value(Phi_x);
u = value(Phi_u);

sls_phixs{ii} = x;
sls_phius{ii} = u;

z = tf('z',h);

disp('making transfer functions')
tic

Phix = x(:,1:n)/z;
for i=2:T
   Phix = Phix + x(:,((i-1)*n+1):(i*n))/(z^i);
end

Phiu = u(:,1:n)/z;
for i=2:T
   Phiu = Phiu + u(:,((i-1)*n+1):(i*n))/(z^i);
end

toc

Phix_temps{ii} = Phix;
Phiu_temps{ii} = Phiu;

end


%% Simple Pole Approximation Method

T = Tlength;

Phixs_pfd = cell(1,length(Ps));
Phius_pfd = cell(1,length(Ps));

% add poles from plant (eig(A)) and desired transfer function (Ades)
evals = unique(eig(A));
evals = [evals;Ades];

for ii=1:length(Ps)
    
poles = Ps(ii)-length(evals);

% select poles from the Archimedes spiral
pu = generate_poles(poles,5);
px = [evals' pu];

% define several constants
Pu = length(pu);
Pe = length(evals);
Pc = 2*Pu;
Px = Pe + Pc;
Pcu = Pe + Pc;

% initialize variables
Phi_x = sdpvar(n,n*Px,'full');
Phi_u = sdpvar(m,n*Pcu,'full');

% construct matrices for impulse response in objective
N = zeros(Pc,T);
Mx = zeros(Px,T);
Mu = zeros(Pcu,T);

for i=1:Pu
    for j=1:T
        N(2*(i-1)+1,j) = 2*real(pu(i)^(j-1));
        N(2*(i-1)+2,j) = -2*imag(pu(i)^(j-1));
    end
end

for i=1:T
    for j=1:Pe
        Mx(j,i) = evals(j)^(i-1);
    end
end
Mx((Pe+1):end,:) = N;

Mu = Mx;

Mx = kron(Mx,Bh);
Mu = kron(Mu,Bh);

coeffs_step = zeros(Px,T);
for i=1:Pe
    for j=1:T
        coeffs_step(i,j) = (1-evals(i)^(j-1))/(1-evals(i));
    end
end

for i=1:Pu
    for j=1:T
        coeffs_step(Pe+2*(i-1)+1,j) = 2*real((1-pu(i)^(j-1))/(1-pu(i)));
        coeffs_step(Pe+2*i,j) = -2*imag((1-pu(i)^(j-1))/(1-pu(i)));
    end
end

coeffs_step_u = coeffs_step;

coeffs_step = kron(coeffs_step,Bh);
coeffs_step_u = kron(coeffs_step_u,Bh);

disp('Making objective');
tic;

Mobj = [QQ*(C*Phi_x*Mx - Tdes_imp);
        R*Phi_u*Mu];
Objective = norm(Mobj,'fro');

disp('Objective made');
toc;


disp('Making constraints');
tic;

temp = 2*ones(2*Pu,1);
temp(2:2:end) = 0;
temp = [ones(Pe,1);temp];
temp = kron(temp,eye(n));

pp = zeros(Px);
for j=1:Pe
    pp(j,j) = evals(j);
end

for i=1:Pu
    pp(Pe+2*(i-1)+1,Pe+2*(i-1)+1) = real(pu(i));
    pp(Pe+2*(i-1)+1,Pe+2*(i-1)+2) = imag(pu(i));
    pp(Pe+2*(i-1)+2,Pe+2*(i-1)+1) = -imag(pu(i));
    pp(Pe+2*(i-1)+2,Pe+2*(i-1)+2) = real(pu(i));
end

pp = kron(pp,eye(n));

Constraints = [Phi_x*pp - A*Phi_x - B*Phi_u == 0, ...
               Phi_x*temp == eye(n)];

disp('Constraints made');
toc

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

% Solve the problem
sol = optimize(Constraints,Objective,options);

obj = value(Objective);

x = value(Phi_x);
u = value(Phi_u);

z = tf('z',h);

Phix = x(:,1:n)/(z-evals(1));
for j=2:Pe
    Phix = Phix+ x(:,((j-1)*n+1):(j*n))/(z-evals(j));
end
for i=1:Pu
    Gx = x(:,((Pe+2*(i-1))*n+1):((Pe+2*(i-1)+1)*n));
    Gy = x(:,((Pe+2*(i-1)+1)*n+1):((Pe+2*i)*n));
    G = Gx + sqrt(-1)*Gy;

    Phix = Phix ...
           + (2*real(G)*z - 2*real(G*conj(pu(i))))/...
           (z^2 - 2*real(pu(i))*z + pu(i)*conj(pu(i)));
    
end

Phiu = u(:,1:n)/(z-evals(1));
for j=2:Pe
    Phiu = Phiu + u(:,((j-1)*n+1):(j*n))/(z-evals(j));
end
add = Pe;

for i=1:Pu
    Hx = u(:,((2*(i-1))*n+1+add):((2*(i-1)+1)*n+add));
    Hy = u(:,((2*(i-1)+1)*n+1+add):((2*i)*n+add));
    H = Hx + sqrt(-1)*Hy;
    
    Phiu = Phiu ...
           + (2*real(H)*z - 2*real(H*conj(pu(i))))/...
           (z^2 - 2*real(pu(i))*z + pu(i)*conj(pu(i)));
end

pfd_cost = norm([QQ*(C*Phix*Bh-Tdes_tf);
                 R*Phiu*Bh],2);

Phixs_pfd{ii} = Phix;
Phius_pfd{ii} = Phiu;

end


%% Plotting Results

% figure 1 shows the impulse responses of the control designs
figure(1)
subplot(2,2,1);
[x,t] = impulse(Tdes_tf(1,1),h*T);
plot(t,x,'k--');
[x,t] = impulse(C(1,:)*Phix_temps{1}*Bh(:,1),h*T);
plot(t,x);
hold on;
[x,t] = impulse(C(1,:)*Phix_temps{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(Topt(1,1),h*T);
plot(t,x);
[x,t] = impulse(Tdes_tf(1,1),h*T);
plot(t,x);
[x,t] = impulse(Topt(1,1),h*T);
x1 = -0.1120;
x2 = 4.2715;
y1 = -0.0585;
y2 = 0.0384;
plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k');
[x,t] = impulse(Tdes_tf(1,1),h*T);
plot(t,x,'k--');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution',...
       'Location','southeast');
axis([-0.2 10 -2.9 0.3]);
xlabel('Time [s]');
ylabel('T_{yw}^{(1,1)}');
fixPlot(1);
subplot(2,2,2);
[x,t] = impulse(Tdes_tf(1,2),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(1,:)*Phix_temps{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phix_temps{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(Topt(1,2),h*T);
plot(t,x);
[x,t] = impulse(Tdes_tf(1,2),h*T);
plot(t,x);
[x,t] = impulse(Topt(1,2),h*T);
plot(t,x);
x1 = -0.1;
x2 = 4;
y1 = -2;
y2 = 10;
plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k');
[x,t] = impulse(Tdes_tf(1,2),h*T);
plot(t,x,'k--');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -30 900]);
xlabel('Time [s]');
ylabel('T_{yw}^{(1,2)}');
fixPlot(1);
subplot(2,2,3);
[x,t] = impulse(Tdes_tf(2,1),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(2,:)*Phix_temps{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phix_temps{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(Topt(2,1),h*T);
plot(t,x);
[x,t] = impulse(Tdes_tf(2,1),h*T);
plot(t,x);
[x,t] = impulse(Topt(2,1),h*T);
plot(t,x);
x1 = -0.2584;
x2 = 4.6950;
y1 = -0.0667;
y2 = 0.0151;
plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k');
[x,t] = impulse(Tdes_tf(2,1),h*T);
plot(t,x,'k--');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -0.2 6]);
xlabel('Time [s]');
ylabel('T_{yw}^{(2,1)}');
fixPlot(1);
subplot(2,2,4);
[x,t] = impulse(Tdes_tf(2,2),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(2,:)*Phix_temps{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phix_temps{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(Topt(2,2),h*T);
plot(t,x);
[x,t] = impulse(Tdes_tf(2,2),h*T);
plot(t,x);
[x,t] = impulse(Topt(2,2),h*T);
plot(t,x);
x1 = -0.1327;
x2 = 4.6648;
y1 = -115.8164;
y2 = 44.2302;
plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'k');
[x,t] = impulse(Tdes_tf(2,2),h*T);
plot(t,x,'k--');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -125  440]);
xlabel('Time [s]');
ylabel('T_{yw}^{(2,2)}');
fixPlot(1);
set(gcf,'Position',[200 100 1174 874]);

% figure 2 shows zoomed-in views of the impulse responses of the control designs
figure(2)
subplot(2,2,1);
[x,t] = impulse(Tdes_tf(1,1),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(1,:)*Phix_temps{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phix_temps{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(Topt(1,1),h*T);
plot(t,x,'b:');
[x,t] = impulse(Tdes_tf(1,1),h*T);
plot(t,x,'k--');
[x,t] = impulse(Topt(1,1),h*T);
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution',...
       'Location','southeast');
axis([-0.1120    4.2715   -0.0585    0.0384]);
xlabel('Time [s]');
ylabel('T_{yw}^{(1,1)}');
fixPlot(2);
subplot(2,2,2);
[x,t] = impulse(Tdes_tf(1,2),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(1,:)*Phix_temps{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phix_temps{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(1,:)*Phixs_pfd{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(Topt(1,2),h*T);
plot(t,x,'b:');
[x,t] = impulse(Tdes_tf(1,2),h*T);
plot(t,x,'k--');
[x,t] = impulse(Topt(1,2),h*T);    
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.1 4 -2 10]);
xlabel('Time [s]');
ylabel('T_{yw}^{(1,2)}');
fixPlot(2);
subplot(2,2,3);
[x,t] = impulse(Tdes_tf(2,1),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(2,:)*Phix_temps{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phix_temps{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = impulse(Topt(2,1),h*T);
plot(t,x,'b:');
[x,t] = impulse(Tdes_tf(2,1),h*T);
plot(t,x,'k--');
[x,t] = impulse(Topt(2,1),h*T);    
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution',...
       'Location','southeast');
axis([-0.2584    4.6950   -0.0667    0.0151]);
xlabel('Time [s]');
ylabel('T_{yw}^{(2,1)}');
fixPlot(2);
subplot(2,2,4);
[x,t] = impulse(Tdes_tf(2,2),h*T);
plot(t,x,'k--');
hold on;
[x,t] = impulse(C(2,:)*Phix_temps{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phix_temps{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(C(2,:)*Phixs_pfd{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = impulse(Topt(2,2),h*T);
plot(t,x,'b:');
[x,t] = impulse(Tdes_tf(2,2),h*T);
plot(t,x,'k--');
[x,t] = impulse(Topt(2,2),h*T);    
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution',...
       'Location','southeast');
axis([-0.1327    4.6648 -115.8164   44.2302]);
xlabel('Time [s]');
ylabel('T_{yw}^{(2,2)}');
fixPlot(2);
set(gcf,'Position',[203 314 1280 957]);

% figure 3 shows the step responses of the control designs
figure(3)
subplot(2,2,1);
[x,t] = step(Tdes_tf(1,1),h*T);
plot(t,x,'k--');
hold on;
[x,t] = step(C(1,:)*Phix_temps{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(C(1,:)*Phix_temps{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(C(1,:)*Phixs_pfd{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(C(1,:)*Phixs_pfd{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(Topt(1,1),h*T);
plot(t,x,'b:');
[x,t] = step(Tdes_tf(1,1),h*T);
plot(t,x,'k--');
[x,t] = step(Topt(1,1),h*T);
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -0.06 0.01]);
xlabel('Time [s]');
ylabel('T_{yw}^{(1,1)}');
fixPlot(3);
subplot(2,2,2);
[x,t] = step(Tdes_tf(1,2),h*T);
plot(t,x,'k--');
hold on;
[x,t] = step(C(1,:)*Phix_temps{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(C(1,:)*Phix_temps{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(C(1,:)*Phixs_pfd{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(C(1,:)*Phixs_pfd{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(Topt(1,2),h*T);
plot(t,x,'b:');
[x,t] = step(Tdes_tf(1,2),h*T);
plot(t,x,'k--');
[x,t] = step(Topt(1,2),h*T);    
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -5 50]);
xlabel('Time [s]');
ylabel('T_{yw}^{(1,2)}');
fixPlot(3);
subplot(2,2,3);
[x,t] = step(Tdes_tf(2,1),h*T);
plot(t,x,'k--');
hold on;
[x,t] = step(C(2,:)*Phix_temps{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(C(2,:)*Phix_temps{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(C(2,:)*Phixs_pfd{1}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(C(2,:)*Phixs_pfd{2}*Bh(:,1),h*T);
plot(t,x);
[x,t] = step(Topt(2,1),h*T);
plot(t,x,'b:');
[x,t] = step(Tdes_tf(2,1),h*T);
plot(t,x,'k--');
[x,t] = step(Topt(2,1),h*T);    
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -0.01 0.06]);
xlabel('Time [s]');
ylabel('T_{yw}^{(2,1)}');
fixPlot(3);
subplot(2,2,4);
[x,t] = step(Tdes_tf(2,2),h*T);
plot(t,x,'k--');
hold on;
[x,t] = step(C(2,:)*Phix_temps{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(C(2,:)*Phix_temps{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(C(2,:)*Phixs_pfd{1}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(C(2,:)*Phixs_pfd{2}*Bh(:,2),h*T);
plot(t,x);
[x,t] = step(Topt(2,2),h*T);
plot(t,x,'b:');
[x,t] = step(Tdes_tf(2,2),h*T);
plot(t,x,'k--');
[x,t] = step(Topt(2,2),h*T);
plot(t,x,'b:');
hold off;
legend('T_{desired}','SLS - 31 poles',...
       'SLS - 300 poles',...
       'SPA - 7 poles','SPA - 15 poles','Optimal Solution');
axis([-0.2 10 -105 10]);
xlabel('Time [s]');
ylabel('T_{yw}^{(2,2)}');
fixPlot(3);
set(gcf,'Position',[203 314 1280 957]);
