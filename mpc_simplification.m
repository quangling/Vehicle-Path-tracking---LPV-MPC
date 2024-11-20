function [Hdb,Fdbt,Cdb,Adc,G,ht] = mpc_simplification(Ad,Bd,Cd,Dd,hz,x_aug_t, du)

% db - double bar
% dbt - double bar transpose
% dc - double circumflex

[A_aug, B_aug, C_aug, D_aug] = augmented_matrices(Ad,Bd,Cd,Dd);

% Get the necessary constants
constants = initial_constants();

Q=constants('Q');
S=constants('S');
R=constants('R');
Cf=constants('Cf');
g=constants('g');
m=constants('m');

lf=constants('lf');
inputs=constants('inputs');

%% Constraints
d_delta_max=pi/50;
d_a_max=0.5;
d_delta_min=-pi/50;
d_a_min=-0.5;

ub_global=zeros(1,inputs*hz);
lb_global=zeros(1,inputs*hz);

% Only works for 2 inputs
for i = 1:inputs*hz
    if mod(i-1,2)==0
        ub_global(i)=d_delta_max;
        lb_global(i)=-d_delta_min;
    else
        ub_global(i)=d_a_max;
        lb_global(i)=-d_a_min;
    end
end

ub_global=ub_global(1:inputs*hz);
lb_global=lb_global(1:inputs*hz);
ublb_global=[ub_global,lb_global];

I_global=eye(inputs*hz);
I_global_negative=-I_global;
I_mega_global=[I_global;I_global_negative];


y_asterisk_max_global=[];
y_asterisk_min_global=[];

C_asterisk=[1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1];
C_asterisk_size=size(C_asterisk);
C_asterisk_global=zeros(C_asterisk_size(1)*hz,C_asterisk_size(2)*hz);

CQC=C_aug'*Q*C_aug;
CSC=C_aug'*S*C_aug;
QC=Q*C_aug;
SC=S*C_aug;

Qdb=zeros(length(CQC(:,1))*hz,length(CQC(1,:))*hz);
Tdb=zeros(length(QC(:,1))*hz,length(QC(1,:))*hz);
Rdb=zeros(length(R(:,1))*hz,length(R(1,:))*hz);
Cdb=zeros(length(B_aug(:,1))*hz,length(B_aug(1,:))*hz);
Adc=zeros(length(A_aug(:,1))*hz,length(A_aug(1,:)));

%% Advanced LPV

A_product=A_aug;
states_predicted_aug=x_aug_t;
A_aug_size=size(A_aug);
B_aug_size=size(B_aug);
A_aug_collection=zeros(A_aug_size(1),A_aug_size(2),hz);
B_aug_collection=zeros(B_aug_size(1),B_aug_size(2),hz);

%% Filling in the matrices
for i = 1:hz
    if i == hz
        Qdb(1+length(CSC(:,1))*(i-1):length(CSC(:,1))*i,1+length(CSC(1,:))*(i-1):length(CSC(1,:))*i)=CSC;
        Tdb(1+length(SC(:,1))*(i-1):length(SC(:,1))*i,1+length(SC(1,:))*(i-1):length(SC(1,:))*i)=SC;
    else
        Qdb(1+length(CQC(:,1))*(i-1):length(CQC(:,1))*i,1+length(CQC(1,:))*(i-1):length(CQC(1,:))*i)=CQC;
        Tdb(1+length(QC(:,1))*(i-1):length(QC(:,1))*i,1+length(QC(1,:))*(i-1):length(QC(1,:))*i)=QC;
    end
    
    Rdb(1+length(R(:,1))*(i-1):length(R(:,1))*i,1+length(R(1,:))*(i-1):length(R(1,:))*i)=R;
    
    Adc(A_aug_size(1)*(i-1)+1:A_aug_size(1)*(i-1)+A_aug_size(1),1:A_aug_size(2))=A_product;
    A_aug_collection(:,:,i)=A_aug;
    B_aug_collection(:,:,i)=B_aug;
    
    %% Constraints
    x_dot_max=30;  
    %if 0.17*states_predicted_aug(1) < 3
    %    y_dot_max=0.17*states_predicted_aug(1);
    %else
    y_dot_max=3;
    %end
    delta_max=pi/6;
    Fyf=Cf*(states_predicted_aug(7)-states_predicted_aug(2)/states_predicted_aug(1)-lf*states_predicted_aug(4)/states_predicted_aug(1));
    %a_max=1+(Fyf*sin(states_predicted_aug(7)))/m-states_predicted_aug(4)*states_predicted_aug(2);
    a_max = 3;
    x_dot_min=1;
    %if -0.17*states_predicted_aug(1) > -3
    %    y_dot_min=-0.17*states_predicted_aug(1);
    %else
    y_dot_min=-3;
    %end
    delta_min=-pi/6;
    %a_min=-4+(Fyf*sin(states_predicted_aug(7)))/m-states_predicted_aug(4)*states_predicted_aug(2);
    a_min = -3;

    y_asterisk_max=[x_dot_max,y_dot_max,delta_max,a_max];
    y_asterisk_min=[x_dot_min,y_dot_min,delta_min,a_min];
    
    y_asterisk_max_global=[y_asterisk_max_global,y_asterisk_max];
    y_asterisk_min_global=[y_asterisk_min_global,y_asterisk_min];
    
    C_asterisk_size=size(C_asterisk);
    C_asterisk_global(C_asterisk_size(1)*(i-1)+1:C_asterisk_size(1)*(i-1)+C_asterisk_size(1), ...
        C_asterisk_size(2)*(i-1)+1:C_asterisk_size(2)*(i-1)+C_asterisk_size(2))=C_asterisk;
    
    %% Advanced LPV
    if i<hz
        du1=du(inputs*i);
        du2=du(inputs*i+inputs-1);
        states_predicted_aug=A_aug*states_predicted_aug+B_aug*[du1,du2]';
        states_predicted=states_predicted_aug(1:6)';
        delta_predicted=states_predicted_aug(7);
        a_predicted=states_predicted_aug(8);
        [Ad, Bd, Cd, Dd]=state_space(states_predicted,delta_predicted,a_predicted);
        [A_aug, B_aug, C_aug, D_aug]=augmented_matrices(Ad, Bd, Cd, Dd);
        A_product=A_aug*A_product;
    end
end

for i = 1:hz
    for j = 1:hz
        if j <= i
            AB_product=eye(A_aug_size(1));
            for ii = i:-1:j
                if ii>j
                    AB_product=AB_product*A_aug_collection(:,:,ii);
                else
                    AB_product=AB_product*B_aug_collection(:,:,ii);
                end
            end
            Cdb(B_aug_size(1)*(i-1)+1:B_aug_size(1)*(i-1)+B_aug_size(1), ...
                B_aug_size(2)*(j-1)+1:B_aug_size(2)*(j-1)+B_aug_size(2))=AB_product;
        end
    end
end

%% Constraints
Cdb_constraints=C_asterisk_global*Cdb;
Cdb_constraints_negative=-Cdb_constraints;
Cdb_constraints_global=[Cdb_constraints;Cdb_constraints_negative];

Adc_constraints=C_asterisk_global*Adc;
Adc_constraints_x0=[Adc_constraints*x_aug_t]';

y_max_Adc_difference=y_asterisk_max_global-Adc_constraints_x0;
y_min_Adc_difference=-y_asterisk_min_global+Adc_constraints_x0;
y_Adc_difference_global=[y_max_Adc_difference,y_min_Adc_difference];

G=[I_mega_global;Cdb_constraints_global];
ht=[ublb_global,y_Adc_difference_global];

Hdb=Cdb'*Qdb*Cdb+Rdb;
Fdbt=[Adc'*Qdb*Cdb;-Tdb*Cdb];

end





%