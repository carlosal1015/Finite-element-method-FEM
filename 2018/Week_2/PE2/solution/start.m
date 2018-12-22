function start()
% M is the number of mesh elements in every direction
% N=(M+1)*(M+1) is the total number of points 
M=5;
NEX = 5;
ERR = zeros(NEX,3);

for ex=1:NEX
    M=M*2;

    % create the mesh
    mesh = createmesh(2,M);
    
    % assemble the matrix
    A = assemblematrix2d(mesh);
    
    % set the right hand side
    f = assemblerhs2d(mesh);
    
    % solve the linear system
    u = A\f;
    
    % compute the l2 and the H1 error
    l2err = l2error(mesh,u);
    h1err = h1error(mesh,u);

    ERR(ex,:) = [1/M,h1err,l2err]';
   % disp(['h / l2 error / h1 error: ',num2str(1/M),'  ',num2str(l2err),' ',num2str(h1err)]);
end

figure('Name','error')
hold on
plot(log(ERR(:,1)),log(ERR(:,2)),'-o');
plot(log(ERR(:,1)),log(ERR(:,3)),'-s');
plot(log(ERR(:,1)),log(ERR(:,1)));
plot(log(ERR(:,1)),log(ERR(:,1).*ERR(:,1)*10));
ERR
%    % show the solution
%plotsolution('solution',mesh,u);
end


