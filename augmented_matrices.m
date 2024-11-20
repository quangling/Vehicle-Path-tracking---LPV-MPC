function [A_aug, B_aug, C_aug, D_aug] = augmented_matrices(Ad,Bd,Cd,Dd)

A_aug=[Ad,Bd;zeros(length(Bd(1,:)),length(Ad(1,:))),eye(length(Bd(1,:)))];
B_aug=[Bd;eye(length(Bd(1,:)))];
C_aug=[Cd,zeros(length(Cd(:,1)),length(Bd(1,:)))];
D_aug=Dd; % D_aug is not used because it is a zero matrix

end