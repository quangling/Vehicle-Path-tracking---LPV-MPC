%% Trajectory generation
function [x_dot_body,y_dot_body,psiInt,X,Y]=trajectory_generator(t)

constants = initial_constants();
Ts=constants('Ts');
trajectory=constants('trajectory');
delay=constants('delay');

if trajectory == 1
    % Trajectory 1
    
    X=15*t;
    Y=750/900.^2*X.^2+250;
elseif trajectory == 2
    % Trajectory 2
    X1=15*t(1:40/Ts);
    Y1=50*sin(2*pi*0.75/40*t(1:40/Ts))+250;

    X2=300*cos(2*pi*0.5/60*(t(40/Ts+1:100/Ts)-40)-pi/2)+600;
    Y2=300*sin(2*pi*0.5/60*(t(40/Ts+1:100/Ts)-40)-pi/2)+500;

    X3=600-15*(t(100/Ts+1:140/Ts+1)-100);
    Y3=50*cos(2*pi*0.75/40*(t(100/Ts+1:140/Ts+1)-100))+750;

    X=[X1,X2,X3];
    Y=[Y1,Y2,Y3];
    
else 
    X=15*t;
    Y=50*sin(2*pi*0.75/20*t)+250;
    
end

 
    

%%
X=round(X,8);
Y=round(Y,8);

dX=X(2:end)-X(1:end-1);
dY=Y(2:end)-Y(1:end-1);

X_dot=dX/Ts;
Y_dot=dY/Ts;

X_dot=[X_dot(1),X_dot];
Y_dot=[Y_dot(1),Y_dot];

X_dot=round(X_dot,8);
Y_dot=round(Y_dot,8);

psi=zeros(1,length(X));
psiInt=psi;
psi(1)=atan2(dY(1),dX(1));
psi(2:end)=atan2(dY(:),dX(:));
dpsi=psi(2:end)-psi(1:end-1);

psiInt(1)=psi(1);
for i = 2:length(psiInt)
    if dpsi(i-1)<-pi
        psiInt(i)=psiInt(i-1)+(dpsi(i-1)+2*pi);
    elseif dpsi(i-1)>pi
        psiInt(i)=psiInt(i-1)+(dpsi(i-1)-2*pi);
    else
        psiInt(i)=psiInt(i-1)+dpsi(i-1);
    end
end


x_dot_body=cos(psiInt).*X_dot+sin(psiInt).*Y_dot;
y_dot_body=-sin(psiInt).*X_dot+cos(psiInt).*Y_dot;

psiInt=round(psiInt,8);
x_dot_body=round(x_dot_body,8);
y_dot_body=round(y_dot_body,8);

if trajectory==2
    % Modify x_dot_body (only for trajectory 2)
    x_dot_body_temp=x_dot_body(40/Ts+1:60/Ts);
    x_dot_body_temp=linspace(x_dot_body_temp(1),16.065,1000);
    x_dot_body(40/Ts+1:60/Ts)=x_dot_body_temp;

    x_dot_body_temp=x_dot_body(80/Ts+3:100/Ts+2);
    x_dot_body_temp=linspace(16.065,x_dot_body_temp(end),1000);
    x_dot_body(80/Ts+3:100/Ts+2)=x_dot_body_temp;

    x_dot_body(60/Ts+1:80/Ts+2)=16.065;

end

end
