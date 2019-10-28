addpath('input_func')
addpath('error')
addpath('matrix_rhs')
addpath('mesh')


clear;
close all;
clc

bc_type = 1;
case_flag=9; % 0 - 4 = poly examples, 5 = sin, 6 free, 7 = exp, 8 = singular, 9= L-shape, 10=cos(y)*sin(x)
% homo Dirichlet: case 4,7,8

basicName=['data/lshape_case_',num2str(case_flag),'_uniform_h_version'];
scaling=1; % scaling of domain = scaling*[-1,1]^2


draw_mesh_flag = 0; % 1 = draw mesh
draw_sol_flag  = 0; % 1 = draw solution

n     = 2; % n= nr. of Elements per edge => n^2 Elemente,
p     = 2; % p=polynomial degree
gqn_f = 5;

k=5 ;  % nr of mesh refinements

for i=1:k
    tic
    disp('----------------------------')
    disp(['      Durchlauf Nr. ', num2str(i)] )
    disp('----------------------------')
    
    name=[basicName,'_p',num2str(p),'_meshNr',num2str(i)];

    % generate mesh
    [element2vertex,vertex_coor,vertex_bc]=MeshGenerator2d_tri_simple_Lshape(n,scaling,bc_type);
    [element2face,face2vertex,face_bc,face2neighbor,face_normal]=generate_faces2(element2vertex,vertex_coor,vertex_bc);
    element2face=check_local_faceNr(element2vertex,element2face,face2vertex);
    disp('mesh done')
    
    % generate connectivity data for u = hat functions (linear) or quadratic nodal basis functions
    [element2connectivity,innerDofs,totaldof]=generate_connectivity_data(element2vertex,element2face,face_bc,vertex_bc,p);
    disp('dofs done')
    
    if draw_mesh_flag==1
        fig_handle=drawMesh(face2vertex,vertex_coor);
        axis equal
        axis off
        saveas(fig_handle,[name,'_mesh.fig'],'fig')
        saveas(fig_handle,[name,'_mesh.jpg'],'jpg')
        close(fig_handle);
        disp('draw mesh done')
    end
    
    % Interpolation of Dirichlet data u on Gamma_D
    G_u=dirichletdaten_interpolation_u(element2connectivity,element2vertex,face2vertex,face_bc,face2neighbor,vertex_coor,@(x,y)func_Diri(x,y,case_flag),innerDofs, totaldof,p);
    
     % volumen term + Neumann term  = rhs
    b=rhsVolume_2d_tri(element2connectivity,element2vertex,vertex_coor,@(x,y)volumeFunction(x,y,case_flag),gqn_f,innerDofs,p);
    b=b+rhsNeumann_2d_tri(element2connectivity,element2vertex,face2vertex,face_bc,face2neighbor,face_normal,vertex_coor,@(x,y,normal)func_NormalGrad_u(x,y,normal,case_flag),gqn_f,innerDofs,p);
    disp('rhs done')
    
    % matrix for u
    [A,b]=stiffnessMatrix_2d_tri_opt(element2connectivity,element2vertex,vertex_coor,innerDofs, b,G_u,p);
    disp('matrix done')
    
    
    % solve linear system
    % symmetric positive definite
    u=A\b;
    
    % include Dirichlet data
    u_all=G_u;
    u_all(1:innerDofs)=u;
    
    % save information for error plots
    dof(i,1)   = length(u);
    disp('solver done')
    time_galerkin(i)=toc;
    
    if draw_sol_flag==1
        % u
        fig_handle=drawSolution(element2connectivity,element2vertex,vertex_coor,u_all,0,p);
        xlabel('x')
        ylabel('y')
        axis auto
        view([-16 36])
        saveas(fig_handle,[name,'_sol_u.fig'],'fig')
        saveas(fig_handle,[name,'_sol_u.jpg'],'jpg')
        close(fig_handle);
        
        % grad u
        fig_handle=drawSolution_gradient(element2connectivity,element2vertex,vertex_coor,u_all,0,scaling,p);
        saveas(fig_handle,[name,'_sol_grad_u.fig'],'fig')
        saveas(fig_handle,[name,'_sol_grad_u.jpg'],'jpg')
        close(fig_handle); 
    end
    
    % compute error
    tic
    [L2_error_u(i,1),H1_error_u(i,1)]=error_u_func(element2connectivity,element2vertex,vertex_coor,u_all,@(x,y)func_u(x,y,case_flag),@(x,y)func_grad_u(x,y,case_flag),gqn_f,p);
    time_error(i)=toc;
        
    save(name)
    n=n*2; % refine mesh
    
    time_error+time_galerkin
end

% compute convergence rates
converL2_u=-log(L2_error_u(2:end)./L2_error_u(1:end-1))./log(dof(2:end)./dof(1:end-1));
converH1_u=-log(H1_error_u(2:end)./H1_error_u(1:end-1))./log(dof(2:end)./dof(1:end-1));
[converL2_u(:) converH1_u(:)]


% plot error curve 
options9 = optimoptions(@lsqnonlin,'TolX',10^-12,'TolFun',10^-12);
high=length(converH1_u); low=max(2,high-9);

my_func=@(x) log(L2_error_u(low:high))-x(1) +x(2)*log(dof(low:high));
xL2=lsqnonlin(my_func,[1;0.5],[],[],options9);

my_func=@(x) log(H1_error_u(low:high))-x(1) +x(2)*log(dof(low:high));
xH1=lsqnonlin(my_func,[1;0.5],[],[],options9);


figure;
loglog(dof,H1_error_u,'-S',dof,  exp(xH1(1)) *dof.^(-xH1(2)) , dof,L2_error_u,'-*',dof,  exp(xL2(1)) *dof.^(-xL2(2)) ,'Markersize',10)
xlabel('Degrees of Freedom','fontsize',16)
legend({'H^1-error',['conver. rate H1 ',num2str(xH1(2))],'L^2-error',['conver. rate L2 ',num2str(xL2(2))]},'location','best','FontSize',16)
