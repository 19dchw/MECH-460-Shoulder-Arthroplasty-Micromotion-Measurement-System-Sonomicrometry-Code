close all;
clear;
format long;

% Define static sensor positions (B and C) in millimetres, scaled to micrometers
scale_factor = 1e3; % Scale to micrometers

% Radius of baseplate

Rg = 22.5*scale_factor;

N = 100;

errs = zeros([N,1]);

Ris = zeros([N,1]);

Ros = zeros([N,1]);

Alphas = zeros([N,1]);

Betas = zeros([N,1]);

hs = zeros([N,1]);

data_basic = zeros([N,14]);
data_const1 = zeros([N,14]);
data_const2 = zeros([N,14]);
data_const3 = zeros([N,14]);

max_perc = zeros(N,1);

for i = 1:N

    % Inner radius of receiver disk
    
    Ri = 32*scale_factor;

    %Ri = Rg + scale_factor*i*10/N;

    Ris(i) = Ri;
    
    % Outer radius of receiver disk

    Ro = Ri + 6*scale_factor;
    
    %Ro = 38*scale_factor;

    %Ro = Ri + 6*scale_factor*i*3/N;

    Ros(i) = Ro;
    
    % Angular offset between inner receivers (deg)
    
    Alpha = atand(0.25);

    %Alpha = Alpha*i/N;

    Alphas(i) = Alpha;

    % Angular offset between outer receivers (deg)

    Beta = atand(4/38);

    %Beta = Beta*i/N;

    Betas(i) = Beta;

    % z-distance from top reference plane to baseplate

    h = 20*scale_factor;

    %h = 8*scale_factor + 32*scale_factor*i/N;

    hs(i) = h;

    % A-receivers
    
    A1 = [Ro;-Ro*tand(Beta);0];
    A2 = [Ri;-Ri*tand(Alpha);0];
    A3 = [Ri;Ri*tand(Alpha);0];
    A4 = [Ro;Ro*tand(Beta);0];
    
    % B-receivers
    
    B1 = [Ro*tand(Beta);Ro;0];
    B2 = [Ri*tand(Alpha);Ri;0];
    B3 = [-Ri*tand(Alpha);Ri;0];
    B4 = [-Ro*tand(Beta);Ro;0];
    
    % C-receivers
    
    C1 = [-Ro;Ro*tand(Beta);0];
    C2 = [-Ri;Ri*tand(Alpha);0];
    C3 = [-Ri;-Ri*tand(Alpha);0];
    C4 = [-Ro;-Ro*tand(Beta);0];
    
    % Initial baseplate
    
    A0_initial = [Rg;0;h];
    B0_initial = [0;Rg;h];
    C0_initial = [-Rg;0;h];
    
    centre_initial = [0;0;h];
    
    % Vectors from centre to sensors
    
    vcAi = A0_initial - centre_initial;
    vcBi = B0_initial - centre_initial;
    vcCi = C0_initial - centre_initial;

    % Set translation and rotation
    
    max_trans = 300;
    
    max_rot = deg2rad(2);
    
    %trans = (rand([3,1])*2 - 1)*max_trans;
    
    %angs = (rand([1,3])*2 - 1)*max_rot;

    trans = [50;70;180];

    angs = deg2rad([1,0.7,0.5]);
    
    rotmat = eul2rotm(angs,'XYZ');
    
    % Apply translation and rotation
    
    centre = trans + centre_initial;
    
    vcA = rotmat*vcAi;
    vcB = rotmat*vcBi;
    vcC = rotmat*vcCi;
    
    A0 = centre + vcA;
    B0 = centre + vcB;
    C0 = centre + vcC;
    
    % OR SET POSITIONS MANUALLY
    
    % A0 = [];
    % B0 = [];
    % C0 = [];

    % Get measured distances with noise
    
    %err = 3 + 0.25*i;

    err = 7;
    
    [rA1,rA2,rA3,rA4] = distances(A0,A1,A2,A3,A4,err);
    
    [rB1,rB2,rB3,rB4] = distances(B0,B1,B2,B3,B4,err);
    
    [rC1,rC2,rC3,rC4] = distances(C0,C1,C2,C3,C4,err);

    rA1 = 12*round(rA1/12);
    rA2 = 12*round(rA2/12);
    rA3 = 12*round(rA3/12);
    rA4 = 12*round(rA4/12);

    rB1 = 12*round(rB1/12);
    rB2 = 12*round(rB2/12);
    rB3 = 12*round(rB3/12);
    rB4 = 12*round(rB4/12);

    rC1 = 12*round(rC1/12);
    rC2 = 12*round(rC2/12);
    rC3 = 12*round(rC3/12);
    rC4 = 12*round(rC4/12);


    % OR SET DISTANCES MANUALLY
    
    % rA1 =  + (rand*2-1)*e;
    % rA2 =  + (rand*2-1)*e;
    % rA3 =  + (rand*2-1)*e;
    % rA4 =  + (rand*2-1)*e;
    % 
    % rB1 =  + (rand*2-1)*e;
    % rB2 =  + (rand*2-1)*e;
    % rB3 =  + (rand*2-1)*e;
    % rB4 =  + (rand*2-1)*e;
    % 
    % rC1 =  + (rand*2-1)*e;
    % rC2 =  + (rand*2-1)*e;
    % rC3 =  + (rand*2-1)*e;
    % rC4 =  + (rand*2-1)*e;
    
    % Find the points on the baseplate using each of the four possible
    % combinations of sensors for each point, then average them
    
    A0_sol1 = trilat(A1,rA1,A2,rA2,A3,rA3,A0_initial);
    A0_sol2 = trilat(A4,rA4,A2,rA2,A3,rA3,A0_initial);
    A0_sol3 = trilat(A1,rA1,A4,rA4,A3,rA3,A0_initial);
    A0_sol4 = trilat(A1,rA1,A2,rA2,A4,rA4,A0_initial);
    
    A0_sol_avg = 0.25*(A0_sol1 + A0_sol2 + A0_sol3 + A0_sol4);
    
    B0_sol1 = trilat(B1,rB1,B2,rB2,B3,rB3,B0_initial);
    B0_sol2 = trilat(B4,rB4,B2,rB2,B3,rB3,B0_initial);
    B0_sol3 = trilat(B1,rB1,B4,rB4,B3,rB3,B0_initial);
    B0_sol4 = trilat(B1,rB1,B2,rB2,B4,rB4,B0_initial);
    
    B0_sol_avg = 0.25*(B0_sol1 + B0_sol2 + B0_sol3 + B0_sol4);
    
    C0_sol1 = trilat(C1,rC1,C2,rC2,C3,rC3,C0_initial);
    C0_sol2 = trilat(C4,rC4,C2,rC2,C3,rC3,C0_initial);
    C0_sol3 = trilat(C1,rC1,C4,rC4,C3,rC3,C0_initial);
    C0_sol4 = trilat(C1,rC1,C2,rC2,C4,rC4,C0_initial);
    
    C0_sol_avg = 0.25*(C0_sol1 + C0_sol2 + C0_sol3 + C0_sol4);
    
    centre_sol = 0.5*(A0_sol_avg+C0_sol_avg);
    
    % Calculate translation
    
    trans_sol = centre_sol - centre_initial;
    
    % Get vectors from centre to sensors. C isn't necessary to determine
    % orientation
    
    vcA_sol = A0_sol_avg - centre_sol;
    vcB_sol = B0_sol_avg - centre_sol;
    
    % Convert to unit vectors in new orientation, and cross product them to get
    % the z unit vector
    
    x_sol_norm = vcA_sol/norm(vcA_sol);
    y_sol_norm = vcB_sol/norm(vcB_sol);
    z_sol_norm = cross(x_sol_norm,y_sol_norm);
    
    % Make them the columns of a matrix: this is now a rotation matrix
    
    rotmat_sol = [x_sol_norm,y_sol_norm,z_sol_norm];
    
    % Get the angles
    
    angs_sol = rad2deg(rotm2eul(rotmat_sol,"XYZ"));
    
    max_error = max([norm(A0_sol_avg-A0),norm(B0_sol_avg-B0),norm(C0_sol_avg-C0)]);
    
    fprintf('Actual translation: [%.6f, %.6f, %.6f] µm\n',trans);
    fprintf('Actual rotation: [%.6f, %.6f, %.6f] deg\n',rad2deg(angs));
    fprintf('Solved translation: [%.6f, %.6f, %.6f] µm\n',trans_sol);
    fprintf('Solved rotation: [%.6f, %.6f, %.6f] deg\n',angs_sol);
    fprintf('Maximum sensor position error: %.6f µm', max_error);
    
    data_basic(i,1) = max_error;
    
    % NOW TRYING TO CONSTRAIN TO BASEPLATE SIZE
    
    low = [A0_sol_avg,B0_sol_avg,C0_sol_avg] - [30,30,30];
    high = [A0_sol_avg,B0_sol_avg,C0_sol_avg] + [30,30,30];
    
    rAB = norm(A0_initial-B0_initial);
    rAC = norm(A0_initial-C0_initial);
    rBC = norm(B0_initial-C0_initial);
    
    [trans_new, fval] = fmincon(@(x)baseplate_dims(x,rAB,rAC,rBC),[A0_sol_avg,B0_sol_avg,C0_sol_avg],[],[],[],[],low,high);
    
    A0_sol_new = trans_new(:,1);
    B0_sol_new = trans_new(:,2);
    C0_sol_new = trans_new(:,3);
    
    centre_sol_new = 0.5*(A0_sol_new+C0_sol_new);
    
    % Calculate translation
    
    trans_sol_new = centre_sol_new - centre_initial;
    
    % Get vectors from centre to sensors. C isn't necessary to determine
    % orientation
    
    vcA_sol_new = A0_sol_new - centre_sol_new;
    vcB_sol_new = B0_sol_new - centre_sol_new;
    
    % Convert to unit vectors in new orientation, and cross product them to get
    % the z unit vector
    
    x_sol_new_norm = vcA_sol_new/norm(vcA_sol_new);
    y_sol_new_norm = vcB_sol_new/norm(vcB_sol_new);
    z_sol_new_norm = cross(x_sol_new_norm,y_sol_new_norm);
    
    % Make them the columns of a matrix: this is now a rotation matrix
    
    rotmat_sol_new = [x_sol_new_norm,y_sol_new_norm,z_sol_new_norm];
    
    % Get the angles
    
    angs_sol_new = rad2deg(rotm2eul(rotmat_sol_new,"XYZ"));
    
    max_error_new = max([norm(A0_sol_new-A0),norm(B0_sol_new-B0),norm(C0_sol_new-C0)]);
    
    %max_error_new = max([norm(A0_sol1-A0),norm(B0_sol1-B0),norm(C0_sol1-C0)]);
    
    fprintf('Solved translation PARTIAL CONSTRAINTS: [%.6f, %.6f, %.6f] µm\n',trans_sol_new);
    fprintf('Solved rotation PARTIAL CONSTRAINTS: [%.6f, %.6f, %.6f] deg\n',angs_sol_new);
    fprintf('Maximum sensor position error PARTIAL CONSTRAINTS: %.6f µm', max_error_new);
    
    data_const1(i,1) = max_error_new;
    data_const1(i,8) = fval;
    
    % MultiStart options and setup
    options = optimoptions('lsqnonlin', 'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);
    initial_guess = [A0_sol_avg, B0_sol_avg, C0_sol_avg]; % Start from average of points
    problem = createOptimProblem('lsqnonlin', ...
        'objective', @(vars) all_distances(vars, [A1,A2,A3,A4], [B1,B2,B3,B4], [C1,C2,C3,C4], [rA1,rA2,rA3,rA4], [rB1,rB2,rB3,rB4], [rC1,rC2,rC3,rC4], rAB, rAC, rBC), ...
        'x0', initial_guess, ...
        'lb', initial_guess - 60, 'ub', initial_guess + 60, ...
        'options', options);
    
    % Run MultiStart with 50 random starting points within a 300 µm radius
    ms = MultiStart('UseParallel', true);
    num_start_points = 300;
    [best_solution, best_fval, exitflag, output, solutions] = run(ms, problem, num_start_points);
    
    % Identify the best solution (smallest residual)
    [~, idx] = min([solutions.Fval]); % Find index of the smallest residual solution
    most_probable_solution = solutions(idx).X; % Best solution coordinates
    
    A0_lsq = most_probable_solution(:,1);
    B0_lsq = most_probable_solution(:,2);
    C0_lsq = most_probable_solution(:,3);
    
    centre_lsq = 0.5*(A0_lsq+C0_lsq);
    
    % Calculate translation
    
    trans_lsq = centre_lsq - centre_initial;
    
    % Get vectors from centre to sensors. C isn't necessary to determine
    % orientation
    
    vcA_lsq = A0_lsq - centre_lsq;
    vcB_lsq = B0_lsq - centre_lsq;
    
    % Convert to unit vectors in new orientation, and cross product them to get
    % the z unit vector
    
    x_lsq_norm = vcA_lsq/norm(vcA_lsq);
    y_lsq_norm = vcB_lsq/norm(vcB_lsq);
    z_lsq_norm = cross(x_lsq_norm,y_lsq_norm);

    % x_init_norm = vcAi/norm(vcAi);
    % y_init_norm = vcBi/norm(vcBi);
    % z_init_norm = cross(x_init_norm,y_init_norm);
    
    % Make them the columns of a matrix: this is now a rotation matrix
    
    rotmat_lsq = [x_lsq_norm,y_lsq_norm,z_lsq_norm];

    %rotmat_lsq_test = [x_lsq_norm,y_lsq_norm,z_lsq_norm]/[x_init_norm,y_init_norm,z_init_norm];
    
    % Get the angles
    
    angs_lsq = rad2deg(rotm2eul(rotmat_lsq,"XYZ"));
    
    max_error_lsq = max([norm(A0_lsq-A0),norm(B0_lsq-B0),norm(C0_lsq-C0)]);

    max_perc(i) = max([abs(1-norm(A0_lsq-A0_initial)./norm(A0-A0_initial)),abs(1 - norm(B0_lsq-B0_initial)./norm(B0-B0_initial)),abs(1 - norm(C0_lsq-C0_initial)./norm(C0-C0_initial))]);

    trans_error_lsq = centre_lsq - centre;

    %rot_error_lsq = angs_lsq - rad2deg(angs);

    rot_error_lsq = rad2deg(rotm2eul(rotmat_lsq / rotmat,'XYZ'));
    
    disp(best_fval)
    
    fprintf('Solved translation BETTER CONSTRAINTS: [%.6f, %.6f, %.6f] µm\n',trans_lsq);
    fprintf('Solved rotation BETTER CONSTRAINTS: [%.6f, %.6f, %.6f] deg\n',angs_lsq);
    fprintf('Maximum sensor position error BETTER CONSTRAINTS: %.6f µm', max_error_lsq);
    
    data_const2(i,:) = [max_error_lsq, trans_error_lsq(1), trans_error_lsq(2), trans_error_lsq(3), rot_error_lsq(1), rot_error_lsq(2), rot_error_lsq(3), best_fval, trans_lsq(1), trans_lsq(2), trans_lsq(3), angs_lsq(1), angs_lsq(2), angs_lsq(3)];  

    % MultiStart options and setup
    options = optimoptions('lsqnonlin', 'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);
    initial_guess = [A0_sol_new, B0_sol_new, C0_sol_new]; % Start from average of points
    problem = createOptimProblem('lsqnonlin', ...
        'objective', @(vars) all_distances(vars, [A1,A2,A3,A4], [B1,B2,B3,B4], [C1,C2,C3,C4], [rA1,rA2,rA3,rA4], [rB1,rB2,rB3,rB4], [rC1,rC2,rC3,rC4], rAB, rAC, rBC), ...
        'x0', initial_guess, ...
        'lb', initial_guess - 60, 'ub', initial_guess + 60, ...
        'options', options);
    
    % Run MultiStart with 50 random starting points within a 300 µm radius
    ms = MultiStart('UseParallel', true);
    num_start_points = 300;
    [best_solution, best_fval, exitflag, output, solutions] = run(ms, problem, num_start_points);
    
    % Identify the best solution (smallest residual)
    [~, idx] = min([solutions.Fval]); % Find index of the smallest residual solution
    most_probable_solution = solutions(idx).X; % Best solution coordinates
    
    A0_lsq = most_probable_solution(:,1);
    B0_lsq = most_probable_solution(:,2);
    C0_lsq = most_probable_solution(:,3);
    
    centre_lsq = 0.5*(A0_lsq+C0_lsq);
    
    % Calculate translation
    
    trans_lsq = centre_lsq - centre_initial;
    
    % Get vectors from centre to sensors. C isn't necessary to determine
    % orientation
    
    vcA_lsq = A0_lsq - centre_lsq;
    vcB_lsq = B0_lsq - centre_lsq;
    
    % Convert to unit vectors in new orientation, and cross product them to get
    % the z unit vector
    
    x_lsq_norm = vcA_lsq/norm(vcA_lsq);
    y_lsq_norm = vcB_lsq/norm(vcB_lsq);
    z_lsq_norm = cross(x_lsq_norm,y_lsq_norm);
    
    % Make them the columns of a matrix: this is now a rotation matrix
    
    rotmat_lsq = [x_lsq_norm,y_lsq_norm,z_lsq_norm];
    
    % Get the angles
    
    angs_lsq = rad2deg(rotm2eul(rotmat_lsq,"XYZ"));
    
    max_error_lsq = max([norm(A0_lsq-A0),norm(B0_lsq-B0),norm(C0_lsq-C0)]);
    
    disp(best_fval)
    
    fprintf('Solved translation EVEN BETTER CONSTRAINTS?: [%.6f, %.6f, %.6f] µm\n',trans_lsq);
    fprintf('Solved rotation EVEN BETTER CONSTRAINTS?: [%.6f, %.6f, %.6f] deg\n',angs_lsq);
    fprintf('Maximum sensor position error EVEN BETTER CONSTRAINTS?: %.6f µm', max_error_lsq);
       
    data_const3(i,1) = max_error_lsq;
    data_const3(i,8) = best_fval;

    errs(i) = err;

    Ris(i) = Ri;

    Ros(i) = Ro;

    Alphas(i) = Alpha;

    Betas(i) = Beta;

end

% Get distances from receivers to transmitter and add error
function [r1,r2,r3,r4] = distances(p0,p1,p2,p3,p4,e)
    r1 = norm(p1-p0) + (rand*2-1)*e;
    r2 = norm(p2-p0) + (rand*2-1)*e;
    r3 = norm(p3-p0) + (rand*2-1)*e;
    r4 = norm(p4-p0) + (rand*2-1)*e;
end

% Trilateration function

function F = trilat(pos1,r1,pos2,r2,pos3,r3,pos0)
    c = optimvar('c',3);
    
    eq1 = (c(1)-pos1(1))^2 + (c(2)-pos1(2))^2 + (c(3)-pos1(3))^2 == r1^2;
    eq2 = (c(1)-pos2(1))^2 + (c(2)-pos2(2))^2 + (c(3)-pos2(3))^2 == r2^2;
    eq3 = (c(1)-pos3(1))^2 + (c(2)-pos3(2))^2 + (c(3)-pos3(3))^2 == r3^2;
    
    prob = eqnproblem;
    prob.Equations.eq1 = eq1;
    prob.Equations.eq2 = eq2;
    prob.Equations.eq3 = eq3;
    
    c0.c = pos0';
    [sol,~,~] = solve(prob,c0);
    
    F(1,:) = sol.c(1);
    F(2,:) = sol.c(2);
    F(3,:) = sol.c(3);
end

function F = baseplate_dims(vars, a0b0, a0c0, b0c0)
    a0 = vars(:,1);
    b0 = vars(:,2);
    c0 = vars(:,3);

    % Set up minimization of deviation from baseplate dimensions

    F = sum([(norm(a0 - b0) - a0b0)^2,(norm(a0 - c0) - a0c0)^2,(norm(b0 - c0) - b0c0)^2]);
end

function G = all_distances(vars,a,b,c,ra,rb,rc,a0b0,a0c0,b0c0)

    a0 = vars(:,1);
    b0 = vars(:,2);
    c0 = vars(:,3);

    a1 = a(:,1);
    a2 = a(:,2);
    a3 = a(:,3);
    a4 = a(:,4);

    b1 = b(:,1);
    b2 = b(:,2);
    b3 = b(:,3);
    b4 = b(:,4);

    c1 = c(:,1);
    c2 = c(:,2);
    c3 = c(:,3);
    c4 = c(:,4);

    a1a0 = ra(1);
    a2a0 = ra(2);
    a3a0 = ra(3);
    a4a0 = ra(4);

    b1b0 = rb(1);
    b2b0 = rb(2);
    b3b0 = rb(3);
    b4b0 = rb(4);

    c1c0 = rc(1);
    c2c0 = rc(2);
    c3c0 = rc(3);
    c4c0 = rc(4);

    G = [
        norm(a0 - b0) - a0b0;
        norm(a0 - c0) - a0c0;
        norm(b0 - c0) - b0c0;

        norm(a0 - a1) - a1a0;
        norm(a0 - a2) - a2a0;
        norm(a0 - a3) - a3a0;
        norm(a0 - a4) - a4a0;

        norm(b0 - b1) - b1b0;
        norm(b0 - b2) - b2b0;
        norm(b0 - b3) - b3b0;
        norm(b0 - b4) - b4b0;

        norm(c0 - c1) - c1c0;
        norm(c0 - c2) - c2c0;
        norm(c0 - c3) - c3c0;
        norm(c0 - c4) - c4c0;
    ];

end

%%
figure
plot(data_basic(:,1)*50)
hold on
plot(data_basic(:,2))

figure
plot(data_const1(:,1)*50)
hold on
plot(data_const1(:,2))

figure
plot(data_const2(:,1)*50)
hold on
plot(data_const2(:,2))

figure
plot(data_const3(:,1)*50)
hold on
plot(data_const3(:,2))

legend()

%%

% plot([data_basic(:,1) data_const1(:,1)])
% hold on

data_basicS = data_basic;
data_const1S = data_const1;
data_const2S = data_const2;
data_const3S = data_const3;

for i = 1:N-1
    for j = 1:N-1
        if data_basicS(j,1) > data_basicS(j+1,1)

            temp = data_basicS(j);
            data_basicS(j) = data_basicS(j+1);
            data_basicS(j+1) = temp;

            temp = data_const1S(j);
            data_const1S(j) = data_const1S(j+1);
            data_const1S(j+1) = temp;

            temp = data_const2S(j);
            data_const2S(j) = data_const2S(j+1);
            data_const2S(j+1) = temp;

            temp = data_const3S(j);
            data_const3S(j) = data_const3S(j+1);
            data_const3S(j+1) = temp;
        end
    end
end

figure
plot(data_basicS(:,1),'.','MarkerSize',10)
hold on
%plot(data_const1S(:,1),'LineWidth',2)
plot(data_const2S(:,1),'.','MarkerSize',10)
%plot(data_const3S(:,1),'LineWidth',2)

title('Comparison of error. Maximum measurement error = ±10µm')
xlabel('Sorted Trial Index')
ylabel('Maximum absolute sensor position error (µm)')

%legend('Averaged','Baseplate Constraint','Least Squares, Averaged Start','Least Squares, Baseplate Start','Location','northwest')
legend('Averaged','Least Squares','Location','northwest')

x = ["Least Squares, Baseplate Start" "Least Squares, Averaged Start" "Baseplate Constraint" "Averaged" ];
y = [mean(data_const3(:,1)) mean(data_const2(:,1)) mean(data_const1(:,1)) mean(data_basic(:,1))];

figure
b = barh(x,y);

b.FaceColor = 'flat';

b.CData(1,:) = [0.4940 0.1840 0.5560];
b.CData(2,:) = [0.9290 0.6940 0.1250];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0 0.4470 0.7410];

title('Comparison of error. Maximum measurement error = ±10µm')
xlabel('Average absolute sensor position error (µm)')
ylabel('Step/Method of Calculation')

set(gcf,'position',[400,200,750,450])

%%

figure

plot(errs,data_const2(:,1),'r.','MarkerSize',15)

hold on

plot(errs,vecnorm(data_const2(:,2:4),2,2),'b.','MarkerSize',15)

%lsline

%plot(tline)

xlabel('Measurement error (µm)')
ylabel('Error of calculated positions (µm)')

legend('Maximum absolute sensor position error', 'Absolute centre position error','Location','northwest')

tline1 = polyfit(errs,abs(data_const2(:,1)),2);
tline2 = polyfit(errs,vecnorm(data_const2(:,2:4),2,2),2);

trendvals1 = polyval(tline1,errs);
trendvals2 = polyval(tline2,errs);

plot(errs,trendvals1,'r-','LineWidth',1)
plot(errs,trendvals2,'b-','LineWidth',1)

figure

plot(errs,abs(data_const2(:,5)),'color','r','LineStyle','none','Marker','.','MarkerSize',15)
hold on
plot(errs,abs(data_const2(:,6)),'color','b','LineStyle','none','Marker','.','MarkerSize',15)
plot(errs,abs(data_const2(:,7)),'Color',[0.4660 0.6740 0.1880],'LineStyle','none','Marker','.','MarkerSize',15)

legend('x-rotation','y-rotation','z-rotation','Location','northwest')

%lsline

tlinex = polyfit(errs,abs(data_const2(:,5)),2);
tliney = polyfit(errs,abs(data_const2(:,6)),2);
tlinez = polyfit(errs,abs(data_const2(:,7)),2);

trendvalsx = polyval(tlinex,errs);
trendvalsy = polyval(tliney,errs);
trendvalsz = polyval(tlinez,errs);

plot(errs,trendvalsx,'r-','LineWidth',1)
plot(errs,trendvalsy,'b-','LineWidth',1)
plot(errs,trendvalsz,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1)


xlabel('Measurement error (µm)')
ylabel('Rotational error (absolute value) (deg)')

%%

figure

ax = axes;

xlabel('Degree of Freedom')

yyaxis left

boxchart([data_const2(:,9:11),NaN(N,3)])

hold on

plot([trans,NaN(3,1)],'k_','MarkerSize',50,'LineWidth',1)

%ylim([-20 210])

%set(gca,'XTickLabel',{'x-translation','y-translation','z-translation'})

ylabel('Translation (µm)')

yyaxis right

boxchart([NaN(N,3),data_const2(:,12:14)])

plot([NaN(1,3),rad2deg(angs)],'k_','MarkerSize',50,'LineWidth',1)

%ylim([0.3 1.2])

set(gca,'XTickLabel',{'x-translation','y-translation','z-translation','x-rotation','y-rotation','z-rotation'})

ylabel('Rotation (deg)')

ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%%

figure

plot(Ris*0.001,data_const2(:,1),'r.')

hold on

plot(Ris*0.001,vecnorm(data_const2(:,2:4),2,2),'b.')

xlabel('Inner receiver distance (mm)')
ylabel('Error of calculated positions (µm)')

legend('Maximum absolute sensor position error', 'Absolute centre position error','Location','northwest')

figure

plot(Ris,data_const2(:,5),'color','r','LineStyle','none','Marker','.')
hold on
plot(Ris,data_const2(:,6),'color','b','LineStyle','none','Marker','.')
plot(Ris,data_const2(:,7),'Color',[0.4660 0.6740 0.1880],'LineStyle','none','Marker','.')

xlabel('Inner receiver distance (µm)')
ylabel('Rotational error (deg)')

legend('x-rotation','y-rotation','z-rotation','Location','northwest')

%%

figure

plot(Ros*0.001,data_const2(:,1),'r.')

hold on

plot(Ros*0.001,vecnorm(data_const2(:,2:4),2,2),'b.')

xlabel('Outer receiver distance (mm)')
ylabel('Error of calculated positions (µm)')

legend('Maximum absolute sensor position error', 'Absolute centre position error','Location','northwest')

figure

plot(Ros,data_const2(:,5),'color','r','LineStyle','none','Marker','.')
hold on
plot(Ros,data_const2(:,6),'color','b','LineStyle','none','Marker','.')
plot(Ros,data_const2(:,7),'Color',[0.4660 0.6740 0.1880],'LineStyle','none','Marker','.')

xlabel('Outer receiver distance (µm)')
ylabel('Rotational error (deg)')

legend('x-rotation','y-rotation','z-rotation','Location','northwest')

%%

figure

plot(Alphas*2,data_const2(:,1),'r.')

hold on

plot(Alphas*2,vecnorm(data_const2(:,2:4),2,2),'b.')

xlabel('Inner receiver angle (deg)')
ylabel('Error of calculated positions (µm)')

legend('Maximum absolute sensor position error', 'Absolute centre position error','Location','northwest')

figure

plot(Alphas*2,data_const2(:,5),'color','r','LineStyle','none','Marker','.')
hold on
plot(Alphas*2,data_const2(:,6),'color','b','LineStyle','none','Marker','.')
plot(Alphas*2,data_const2(:,7),'Color',[0.4660 0.6740 0.1880],'LineStyle','none','Marker','.')

xlabel('Inner receiver angle (deg)')
ylabel('Rotational error (deg)')

legend('x-rotation','y-rotation','z-rotation','Location','northwest')

%%

figure

plot(Betas*2,data_const2(:,1),'r.')

hold on

plot(Betas*2,vecnorm(data_const2(:,2:4),2,2),'b.')

xlabel('Outer receiver angle (deg)')
ylabel('Error of calculated positions (µm)')

legend('Maximum absolute sensor position error', 'Absolute centre position error','Location','northwest')

figure

plot(Betas*2,data_const2(:,5),'color','r','LineStyle','none','Marker','.')
hold on
plot(Betas*2,data_const2(:,6),'color','b','LineStyle','none','Marker','.')
plot(Betas*2,data_const2(:,7),'Color',[0.4660 0.6740 0.1880],'LineStyle','none','Marker','.')

xlabel('Outer receiver angle (deg)')
ylabel('Rotational error (deg)')

legend('x-rotation','y-rotation','z-rotation','Location','northwest')

%%

figure

plot(hs*0.001,data_const2(:,1),'r.')

hold on

plot(hs*0.001,vecnorm(data_const2(:,2:4),2,2),'b.')

xlabel('Vertical distance between transmitter and receiver planes (mm)')
ylabel('Error of calculated positions (µm)')

legend('Maximum absolute sensor position error', 'Absolute centre position error','Location','northwest')

figure

plot(hs,data_const2(:,5),'color','r','LineStyle','none','Marker','.')
hold on
plot(hs,data_const2(:,6),'color','b','LineStyle','none','Marker','.')
plot(hs,data_const2(:,7),'Color',[0.4660 0.6740 0.1880],'LineStyle','none','Marker','.')

xlabel('Vertical distance between transmitter and receiver planes (µm)')
ylabel('Rotational error (deg)')

legend('x-rotation','y-rotation','z-rotation','Location','northwest')

%%

disp(mean(max_perc))