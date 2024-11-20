function constants=initial_constants()
    
    % Constants    
    
    g=9.81;
    m=1500;
    Iz=3000;
    Cf=38000;
    Cr=66000;    
    rho=1.225;
    lf=2;
    lr=3;
    Ts=0.02;
    
    outputs=4;
    inputs=2;
    
    % Choose your trajectory (1,2,3)
    trajectory=1;

    delay=0;
    if trajectory==1
        hz=10;
        time_length=60;
        Q =[1 0 0 0;0 50 0 0;0 0 50 0;0 0 0 50];
        S=[1 0 0 0;0 200 0 0;0 0 50 0;0 0 0 50];
        R=[100 0;0 1];
    elseif trajectory==2
        hz=10;
        time_length=140;
        Q=[10000 0 0 0;0 10000 0 0;0 0 100 0;0 0 0 100];
        S=[10000 0 0 0;0 10000 0 0;0 0 100 0;0 0 0 100];
        R=[100 0;0 1];
    else
        hz=10;
        time_length=60;
        Q =[1 0 0 0;0 50 0 0;0 0 50 0;0 0 0 50];
        S=[1 0 0 0;0 200 0 0;0 0 50 0;0 0 0 50];
        R=[100 0;0 1];   
    
    end
    
    
    keySet = {'g', 'm', 'Iz', 'Cf', 'Cr', 'rho', 'lf', 'lr', 'Ts', 'outputs', 'inputs', 'hz', 'trajectory', 'Q', 'S', 'R','time_length','delay'};
    constants_list={g m Iz Cf Cr rho lf lr Ts  outputs inputs hz trajectory Q S R time_length delay};
    constants = containers.Map(keySet,constants_list);

end