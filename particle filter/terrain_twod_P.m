x0=[0 0 1 0];
dt = 0.01;
err1=0.01*(rand(1,1001));
err1=err1-mean(err1);
err2=0.01*(rand(1,1001));
err2=err2-mean(err2);
err3=0.001*(rand(1,1001));
err3=err3-mean(err3);
phi=[1 dt 0 0;0 1 0 0;0 0 1 dt;0 0 0 1];
G=[dt^2/2 dt 0 0;0 0 dt^2/2 dt]';
t = 0:dt:(1000)*dt;
z = [5.*t+err1;t.^2/2+err2];
zg=t.^3/6+25/4*t.^4+err3;
% err = 0.01*sin(100.*t);

H = [0 1 0 0;0 0 0 1];
ax = 5+zeros(1,numel(t));
az=t;
a=[ax; az];

X = zeros(numel(x0),numel(err1));
X(:,1) = x0;
R=diag([var(err1) var(err2)]);            %R= covariance of the measurement noise

P=eye(4); %P= covariance error matrix
N = 40;                %N= number of random particles

i = 1;
Pm = zeros(4,4,N);     %Pm=P(-), the best covariance error matrix
while i < N+1
    Pm(:,:,i) = P;
    i = i+1;
end
%% 
k = 2;
X(:,1) = x0;
chi=zeros(4,N);
chi(1,:) = rand(1,N);  %chi= random number
chi(1,:) = chi(1,:) - mean(chi(1,:));
chi(1,:) = chi(1,:)*sqrt(P(1,1)/var(chi(1,:)))+x0(1);
chi(2,:) = rand(1,N);
chi(2,:) = chi(2,:) - mean(chi(2,:));
chi(2,:) = chi(2,:)*sqrt(P(2,2)/var(chi(2,:)))+x0(2);
chi(3,:) = rand(1,N);
chi(3,:) = chi(3,:) - mean(chi(3,:));
chi(3,:) = chi(3,:)*sqrt(P(3,3)/var(chi(3,:)))+x0(3);
chi(4,:) = rand(1,N);
chi(4,:) = chi(4,:) - mean(chi(4,:));
chi(4,:) = chi(4,:)*sqrt(P(4,4)/var(chi(4,:)))+x0(4);
w = 1/N+zeros(1,N);    %w= weight assign to random particles
Kp = zeros(4,2,N);       %Kp= kalman gain for each random particle
zb = zeros(2,N);       %zb= updated measurement
while k < numel(t)+1
    i = 1;
    while i < N+1
        Kp(:,:,i) = Pm(:,:,i)*H'/(H*Pm(:,:,i)*H'+R);
        chi(:,i) = phi*chi(:,i)+G*a(:,k);
        Pm(:,:,i) = phi*Pm(:,:,i)*phi';
        i = i+1;
    end
    i = 1;
    %K = Kp*w';
    while i < N+1
        zb(1,i) = z(1,k) - H(1,:)*chi(:,i);
        zb(2,i) = z(2,k) - H(2,:)*chi(:,i);
        
        chi(:,i) = chi(:,i) + Kp(:,:,i)*zb(:,i);
        P(:,:,i) = Pm(:,:,i) - Kp(:,:,i)*H*Pm(:,:,i); % Pm(:,:,i) = Pm(:,:,i) - K*H*Pm(:,:,i)
        w(i) = 1/(sqrt(2*pi)*exp((chi(3,i)-(1-(chi(1,i))^2)-zg(k))^2/2));
        i = i+1;
        
    end
    w = w./sum(w);
    X(:,k) = chi*w';
    i = 1;
    P = zeros(4,4);
    while i < N+1
        P = P + w(i)*(chi(:,i)-X(:,k))*(chi(:,i)-X(:,k))';
        i = i+1;
    end
%     plot(chi(1,:),w);
%     Xu(1,k) = X(1,k)+sqrt(abs(P(1,1)));
%     Xu(2,k) = X(2,k)+sqrt(abs(P(2,2)));
%     Xu(3,k) = X(3,k)+sqrt(abs(P(3,3)));
%     Xl(1,k) = X(1,k)-sqrt(abs(P(1,1)));
%     Xl(2,k) = X(2,k)-sqrt(abs(P(2,2)));
%     Xl(3,k) = X(3,k)-sqrt(abs(P(3,3)));
    %i = 1;
    %Pm = zeros(3,3,N);
    %while i < N+1
    %    Pm(:,:,i) = P;
    %    i = i+1;
    %end
chi(1,:) = rand(1,N);  %chi= random number
chi(1,:) = chi(1,:) - mean(chi(1,:));
chi(1,:) = chi(1,:)*sqrt(P(1,1)/var(chi(1,:)))+X(1,k);
chi(2,:) = rand(1,N);
chi(2,:) = chi(2,:) - mean(chi(2,:));
chi(2,:) = chi(2,:)*sqrt(P(2,2)/var(chi(2,:)))+X(2,k);
chi(3,:) = rand(1,N);
chi(3,:) = chi(3,:) - mean(chi(3,:));
chi(3,:) = chi(3,:)*sqrt(P(3,3)/var(chi(3,:)))+X(3,k);
chi(4,:) = rand(1,N);
chi(4,:) = chi(4,:) - mean(chi(4,:));
chi(4,:) = chi(4,:)*sqrt(P(4,4)/var(chi(4,:)))+X(4,k);
    w = 1/N+zeros(1,N);
    zb = zeros(1,N);
    k = k+1;
end
figure
subplot(2,2,1)
plot (t,X(1,:)-5/2*t.^2);
title('X-position');
subplot(2,2,2)
plot(t,X(2,:)-5*t);
title('X-velocity');
subplot(2,2,3)
plot(t,X(3,:)-t.^3/6)
title('Z-position')
subplot(2,2,4)
plot(t,X(4,:)-t.^2/2)
title('Z-velocity')
% subplot(2,2,3)
% plot(t,X(3,:),t,Xu(3,:),t,Xl(3,:));
% title('acceleration');
% subplot(2,2,4)
% plot(t,z);
% title('measurement');
