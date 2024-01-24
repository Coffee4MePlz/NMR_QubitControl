function res = PlotBloch(rho_in, numfig, npts, cor_ini, animate)
if ~iscell(rho_in)
    disp("ERROR: rho needs to be a cell")
    return
end

% Definicoes uteis
id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
% Rotacoes
Rotx = @(thet1) expm(-1i*thet1*sx/2.0);
Roty = @(thet2) expm(-1i*thet2*sy/2.0);
Rotz = @(thet3) expm(-1i*thet3*sz/2.0);


figure(numfig)
% Plota a esfera
h = plot_bloch_sphere();

% get rho size for looping
[p1,q1] = size(rho_in);

for n_state = 1:p1
    if isempty(cor_ini)
        cor_ini = rand(1,3); % cor inicial RGB
    end
    for n_state_evol = 1:q1
        if n_state_evol==1
            is_ini=1;
        elseif n_state_evol ==q1
            is_ini=2;
        else 
            is_ini=3;
        end
        % Plota vetor inicial
        rho = rho_in{n_state,n_state_evol};
        plot_bloch_vector(rho,cor_ini,is_ini);
        if n_state_evol < q1
            link_bloch_vectors(rho_in{n_state,n_state_evol}, rho_in{n_state,n_state_evol+1},cor_ini, npts)
        end
    end
end

hold on

if animate==true
    % Identifica o estado inicial
    b_v = [trace(sx*rho) trace(sy*rho) trace(sz*rho)];
    text(b_v(1),b_v(2),b_v(3)+sign(b_v(3))*0.1,...
      '$\rho(0)$','Interpreter','latex','color',cor_ini,'FontSize',13)
    
    % Hamiltoniano de 2 qubtis em termos de matrizes de Pauli
    %H_imput = [1 (0.2-0.2i);(0.2+0.2i) 0];
    %r_H = Hamiltonian_2_Pauli(H_imput)
    
    % Evolucao
    r_H = 0.5*[0,0,1,0]; % definine H em termos das matrizes de Pauli | hbar = 1
    H_s = (r_H(4)*id + r_H(1)*sx + r_H(2)*sy + r_H(3)*sz);
    t_final  = 5*pi; % Instante de tempo final
    
    % Operador evolucao
    U_t = @(t) expm(-1i*t*H_s); % hbar = 1
    
    % Plota evolucao
    color_trace = [1,0.6,0.6];
    npontos = 40;
    plot_volution_bloch_vector(rho,color_trace,npontos,t_final,H_s,2)
    
    % Estado Final
    %ket_final = U_t(t_final)*ket;
    rho = U_t(t_final)*rho*U_t(t_final)';
    cor_fin = [0,0,1]; % cor inicial RGB
    % Plota vetor final
    %rho = kron(ket_final,ket_final');
    plot_bloch_vector(rho,cor_fin,3);
    % Identifica estado final
    b_v = [trace(sx*rho) trace(sy*rho) trace(sz*rho)];
    text(b_v(1),b_v(2),b_v(3)+sign(b_v(3))*0.1,...
      '$\rho(t)$','Interpreter','latex','color',cor_fin,'FontSize',13);


    % --------------------------
    % rho_in, t_final, H_s
    
    nfig = numfig+10;
    npontos = 400;
    flag = 1; % flag = 1 -> unitary evolution | flag = 2 -> dephasing 
    t_final = 5*2*pi;
    
    Bloch_animation(nfig,rho,H_s,npontos,t_final,flag)
    % --------------------------
end
end

function Bloch_animation(nfig,rho_in_mat,H_s,npontos,t_final,flag)
% flag = 1 -> unitary evolution | flag = 2 -> dephasing 
id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];

% Setting up the Plot
figure(nfig)
hold on;

% Creating Data to Animate
% Time array
aat = linspace(0, t_final, npontos);
%rho_in_mat = rho_in{1};
% Evolved Bloch vectoor coordinates
for j = 1:length(aat)
  if flag == 2
    rho_t = expm(-1i*aat(j)*H_s)*rho_in_mat*expm(1i*aat(j)*H_s);
    p = 1-exp(-aat(j)*4/t_final);
    rho_t = dephasing_ch(rho_t, p);
  else
    rho_t = expm(-1i*aat(j)*H_s)*rho_in_mat*expm(1i*aat(j)*H_s);
  end
  ax(j) = real(trace(rho_t*sx));
  ay(j) = real(trace(rho_t*sy));
  az(j) = real(trace(rho_t*sz));

  Tr_rho_2(j) = trace(rho_t^2); 

end


% Plot Bloch sphere
    title(sprintf('Trajectory\nTime: %0.2f msec  Purity = %0.2f ', aat(1), Tr_rho_2(1)), 'Interpreter', 'Latex');
h1 = plot_bloch_sphere();

% Create file name variable
filename = 'Block_animation.gif';

% Plotting with no color to set axis limits
plot3(ax,ay,az,'Color','none');

%view(-37.5,30);  % Setting viewing angle

% Plotting the first iteration
p = plot3(ax(1),ay(1),az(1),'-','Color','b'); % line
m = scatter3(ax(1),ay(1),az(1),'o','filled','r'); % points
v = plot3([0,ax(1)],[0,ay(1)],[0,az(1)],'Color',[1,0.1,0.1],'LineWidth',1.5); % Block vector
scatter3(ax(1),ay(1),az(1),'o','filled','Color',[0.6,0.6,0.6]); % mark the starting point
% Iterating through the length of the time array
for k = 1:length(aat)
    % Updating the line
    p.XData = ax(1:k);
    p.YData = ay(1:k);
    p.ZData = az(1:k);
    % Updating the point
    m.XData = ax(k); 
    m.YData = ay(k);
    m.ZData = az(k);
    % Updating the Block vector
    v.XData = [0,ax(k)]; 
    v.YData = [0,ay(k)];
    v.ZData = [0,az(k)];
    % Updating the title
    title(sprintf('Trajectory\nTime: %0.2f msec  Purity = %0.2f ', aat(k), Tr_rho_2(k)), 'Interpreter', 'Latex');
    % Delay
    pause(0.01)
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end
end

end

% --------------------------

function h = plot_bloch_sphere()
% PLOT_BLOCH_SPHERE  Bloch sphere plot.
%  h = plot_bloch_sphere()
%  Plots a Bloch sphere, a geometrical representation of the state space of a single qubit.
%  Pure states are on the surface of the sphere, nonpure states inside it.
%  The states |0> and |1> lie on the north and south poles of the sphere, respectively.

[X,Y,Z] = sphere(40);

hold off;
h = surf(X,Y,Z, 2*ones(41,41));
hold on;
shading flat
alpha(0.2)
axis equal
grid on
xlabel('x');
ylabel('y');
zlabel('z');

axis_color = [0.5,0.5,0.5];

% Acerta a fonte dos eixos
ax = gca;
ax.FontSize = 13;

% line in the equator
phi = linspace(0, 2*pi, 40);
plot3(cos(phi), sin(phi), zeros(size(phi)), '-','Color', axis_color);

% axis
line([-1,1], [0,0], [0,0], 'LineWidth',1, 'Color', axis_color)
line([0,0], [-1,1], [0,0], 'LineWidth',1, 'Color', axis_color)
line([0,0], [0,0], [-1,1], 'LineWidth',1, 'Color', axis_color)

%axis labes
text(-1.2,0,0, '$\left|x_-\right>$','Interpreter','latex','FontSize',16)
text(1.2,0,0, '$\left|x_+\right>$','Interpreter','latex','FontSize',16)
text(0,-1.2,0, '$\left|y_-\right>$','Interpreter','latex','FontSize',16)
text(0,1.2,0, '$\left|y_+\right>$','Interpreter','latex','FontSize',16)
text(0,0,-1.15, '$\left|1 \right>$','Interpreter','latex','FontSize',16)
text(0,0,1.15, '$\left|0\right>$','Interpreter','latex','FontSize',16)

view(130,20); % angular position for vizualization

end

function link_bloch_vectors(rho1,rho2, color_vec, npts)
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];

blochvec1 = real([trace(sx*rho1) trace(sy*rho1) trace(sz*rho1)]);
blochvec2 = real([trace(sx*rho2) trace(sy*rho2) trace(sz*rho2)]);

[sph1(1), sph1(2),sph1(3)] = cart2sph(blochvec1(1),blochvec1(2),blochvec1(3));
[sph2(1), sph2(2), sph2(3)] = cart2sph(blochvec2(1),blochvec2(2),blochvec2(3));
lambda = linspace(0,1,npts);

for i =1:npts
    sph_pts(i,:) = lambda(i).*sph1 + (1-lambda(i)).*sph2;
    param = sph1(1)-sph2(1);
    if abs(param)>pi
        if sign(param)>0
            sph_pts(i,1) = sph2(1) + (lambda(i))*(sph1(1)-sph2(1) - 2*pi);
        elseif sign(param)<0
            sph_pts(i,1) = sph2(1) + (lambda(i))*(sph1(1)-sph2(1) + 2*pi);
        end
    end
    [Bsph(i,1), Bsph(i,2),Bsph(i,3)] = sph2cart(sph_pts(i,1), sph_pts(i,2),sph_pts(i,3));
    plot3(Bsph(i,1),Bsph(i,2),Bsph(i,3),'.','Color',color_vec,'MarkerSize',7,...
    'MarkerFaceColor',color_vec)
end



end

function plot_bloch_vector(rho,color_vec, is_ini)
%  Bloch vector plot
%  plot_bloch_vector(rho,color_vec)
%  rho is pure or mixed density operator
%  color RGB format [r,g,b]
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];

blochvec = [trace(sx*rho) trace(sy*rho) trace(sz*rho)];

% uncomment this later
%line([0,blochvec(1)], [0,blochvec(2)], [0,blochvec(3)], 'LineWidth',1.5, 'Color', color_vec)
line(blochvec(1), blochvec(2), blochvec(3), 'LineWidth',1.5, 'Color', color_vec)

%plot3(blochvec(1),blochvec(2),blochvec(3),'o','Color',color_vec,'MarkerSize',7,...
%    'MarkerFaceColor',color_vec)

if is_ini==1
    %text(blochvec(1)*1.051,blochvec(2)*1.051,blochvec(3)*1.051, '$\rho(0)$','Interpreter','latex','FontSize',16)
    plot3(blochvec(1),blochvec(2),blochvec(3),'o','Color',color_vec,'MarkerSize',7,...
    'MarkerFaceColor',color_vec)
elseif is_ini==2
    %text(blochvec(1)*1.051,blochvec(2)*1.05,blochvec(3)*1.051, '$\rho(T)$','Interpreter','latex','FontSize',16)
    plot3(blochvec(1),blochvec(2),blochvec(3),'o','Color',color_vec,'MarkerSize',7,...
    'MarkerFaceColor',color_vec)
else
    plot3(blochvec(1),blochvec(2),blochvec(3),'x','Color',color_vec,'MarkerSize',7,...
    'MarkerFaceColor',color_vec)

end

end

function plot_volution_bloch_vector(rho_0,color_traj,npontos,t_final,H_s,flag)
%  Bloch vector unitary evolution plot
%  plot_volution_bloch_vector(rho,color_trace,npontos,t_final,H_s,flag)
%  rho is pure or mixed density operator
%  color_traj in RGB format [r,g,b]
%  npontos -> number of points in the evolution
%  t_final -> time interval
%  H_s -> 1 qubit Hamiltonian 
%  Trajectory specification: flag = 1 --> dots | flag = 2 --> line

id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];

dt = (t_final)/(npontos);

if flag == 1
  tt = dt;
  n_aux = npontos - 1;
else 
  tt = 0;
  n_aux = npontos+1;
end

U_t = @(t) expm(-1i*t*H_s);

for j = 1:n_aux
  U_aux = U_t(tt);
  rho = U_aux*rho_0*ctranspose(U_aux);
  blochvec = [real(trace(sx*rho)) real(trace(sy*rho)) real(trace(sz*rho))];
  tt = tt + dt;
  if flag == 1
    plot3(blochvec(1),blochvec(2),blochvec(3),'o','Color',color_traj,'MarkerSize',4,...
      'MarkerFaceColor',color_trace)
  else
    x(j) = blochvec(1);  y(j) = blochvec(2);   z(j) = blochvec(3); 
  end
end

if flag == 2
  plot3(x,y,z,'-.','Color',color_traj,'LineWidth',1,...
    'MarkerFaceColor',color_traj)
end

end


function r = Hamiltonian_2_Pauli(H_imput)
%  Express 1 qubit Hamiltonian in terms of Pauli matrix 
%  H_imput -> hamiltonian in the computational basis
%  r -> array of din 4 such that
%  H = r(4)*id + r(1)*s_x + r(2)*s_y + r(3)*s_z
%  if H is not Hermitian a NaN vector is returned 

if H_imput == ctranspose(H_imput)  
  r(4) = (H_imput(1,1) + H_imput(2,2))/2.0;
  r(3) = (H_imput(1,1) - H_imput(2,2))/2.0;
  r(2) = imag(H_imput(2,1));
  r(1) = real(H_imput(2,1));
else
  r = [Nan Nan Nan NaN]; % H is not Hermitian
end

end

function rho_out = dephasing_ch(rho_in, p)
id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0]; sx = [0 1; 1 0]; sy = [0 -1i; 1i 0];

% Kraus operators
M{1} = sqrt(1 - p/2)*id;
M{2} = sqrt(p/2)*sz;

rho_out = 0.*id;
for k = 1:length(M)
    rho_out = rho_out + M{k}*rho_in*ctranspose(M{k});
end

end


function plot_volution_bloch_vector_deph(rho_0,color_traj,npontos,t_final,H_s,flag)
%  Bloch vector unitary evolution plot
%  plot_volution_bloch_vector(rho,color_trace,npontos,t_final,H_s,flag)
%  rho is pure or mixed density operator
%  color_traj in RGB format [r,g,b]
%  npontos -> number of points in the evolution
%  t_final -> time interval
%  H_s -> 1 qubit Hamiltonian 
%  Trajectory specification: flag = 1 --> dots | flag = 2 --> line

id = [1.0 0.0; 0.0 1.0];
sz = [1.0 0.0; 0.0 -1.0];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];

dt = (t_final)/(npontos-1);

if flag == 1
  tt = dt;
  n_aux = npontos;
else 
  tt = 0;
  n_aux = npontos;
end

U_t = @(t) expm(-1i*t*H_s);

for j = 1:n_aux
  U_aux = U_t(tt);
  rho_in = U_aux*rho_0*ctranspose(U_aux);
  p_aux = 2*(1 - exp(-tt/t_final));
  rho = dephasing_ch(rho_in, p_aux);
  blochvec = [real(trace(sx*rho)) real(trace(sy*rho)) real(trace(sz*rho))];
  tt = tt + dt;
  if flag == 1
    plot3(blochvec(1),blochvec(2),blochvec(3),'o','Color',color_traj,'MarkerSize',4,...
      'MarkerFaceColor',color_trace)
  else
    x(j) = blochvec(1);  y(j) = blochvec(2);   z(j) = blochvec(3); 
  end
end


if flag == 2
  plot3(x,y,z,'-.','Color',color_traj,'LineWidth',1,...
    'MarkerFaceColor',color_traj)
end

end
