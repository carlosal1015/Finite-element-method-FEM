function val=func_NormalGrad_u(x,y,normal,flag)

grad=func_grad_u(x,y,flag);

val=grad(:,1)*normal(1)+grad(:,2)*normal(2);

end