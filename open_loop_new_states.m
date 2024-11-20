function dx = open_loop_new_states(t, states, U)

    % Get the necessary constants
    
    constants = initial_constants();
    
    g=constants('g');
    m=constants('m');
    Iz=constants('Iz');
    Cf=constants('Cf');
    Cr=constants('Cr');
    lf=constants('lf');
    lr=constants('lr');
    Ts=constants('Ts');
    
    x_dot=states(1);
    y_dot=states(2);
    psi=states(3);
    psi_dot=states(4);
    X=states(5);
    Y=states(6);
    
    % Inputs:
    delta  = U(1);
    a  = U(2);
    
    Fyf=Cf*(delta-y_dot/x_dot-lf*psi_dot/x_dot);
    Fyr=Cr*(-y_dot/x_dot+lr*psi_dot/x_dot);
   
    dx(1,1)=a+(-Fyf*sin(delta))/m+psi_dot*y_dot;
    dx(2,1)=(Fyf*cos(delta)+Fyr)/m-psi_dot*x_dot;
    dx(3,1)=psi_dot;
    dx(4,1)=(Fyf*lf*cos(delta)-Fyr*lr)/Iz;
    dx(5,1)=x_dot*cos(psi)-y_dot*sin(psi);
    dx(6,1)=x_dot*sin(psi)+y_dot*cos(psi);
    
end