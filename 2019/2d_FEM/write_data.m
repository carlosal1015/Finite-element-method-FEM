
fid=fopen('lshape_case_9_uniform_h_version_p1_meshNr10.txt','w');
for i=1:length(dof)
    fprintf(fid,'%9i %1.15e %1.15e %1.15e \n',dof(i),sqrt(2)/2^(i-1),H1_error_u(i),L2_error_u(i));
end
fclose(fid);

% 0 dof
% 1 h
% 2 H1_error_u
% 3 L2_error_u


