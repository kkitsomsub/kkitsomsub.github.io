clear all
M_total= 1%0.00567 %in KG
g=9.8
S=1000 %stiffness 
D=0.5 %damping
mu= 0.8 %sliding friction
tmax=16*10
clockmax= 1*16*10^5
dt=tmax/clockmax
tsave=zeros(1,clockmax);
%defined a rotation matrix
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)]
angle = pi/12 %radians, works = pi/8,pi/4, doesnot = pi/12 
R = roty(angle);
%%%Coin
    num_points = 120;  %n points on circle
    radius = 0.5;
    initial_pos = [0,0,radius*cos(angle)];
    X = zeros(num_points,3); %position of each node
    U = zeros(num_points,3); %velocity of each node
    %U = [0.5,0.5,0].*ones(num_points,3);
    %U = ones(num_points,3);
    M=(M_total/num_points)*ones(num_points,1); %mass of each node
    
    % Generating points for the ball
        for num = 1:num_points
            theta = num * 2 * pi / num_points;
            X(num,:) = [0,radius*cos(theta),radius*sin(theta)];
            M(num) =  M_total / (num_points);   
        end
        %Introduce a rotation
        X = (R*X')';
        X = X+initial_pos; %Translate ball to starting position
        
        %create L
        %L = *(cross(X(1,:)-X(16,:),X(2,:)-X(17,:))/norm(cross(X(1,:)-X(16,:),X(2,:)-X(17,:))))'
        L = 1*(R*[1 0 0]');
        Xcm = (1/M_total)*sum((M.*X)); %Center of Mass
        Ucm = (1/M_total)*sum((M.*U)); %Velocity rel to COM.
        Xrel = X-Xcm;
        Urel = U-Ucm;
    
%     %Creating Angular Velocity
%     for k = 1:num_points
%         L = L + M(k)*cross(Xrel(k,:),U(k,:))';
%     end	
%     


%Output
XcmSaveY=zeros(1,clockmax); %Save the position of center of mass
XcmSaveZ=zeros(1,clockmax);
Esave=zeros(1,clockmax); %Save the total energy in the system
Ymin=zeros(1,clockmax); %Save the position of the top point
Xmin = zeros(1,clockmax); %Save the position of the bottom point
UtanSave = zeros(2,clockmax);

v = VideoWriter('subplus, angle = 15, D=0.5,mu = 0.1, spinning 1.5x Angular Momentum  ');
open(v);

for clock=1:clockmax
  % Compute I
  I = zeros(3, 3);
  for k = 1:num_points
    I = I + M(k).*( (norm(Xrel(k,:))^2).*eye(3) - Xrel(k,:)'*Xrel(k,:) );
  end
    
  % Solve System for Omega, Update Xrel,Urel
  Omega = I\L;
  L;
  norm(Omega);
  %if(norm(Omega) > 100*eps)
      unit_Omega = Omega/norm(Omega);
      Omega_cross = [0 -Omega(3) Omega(2); Omega(3) 0 -Omega(1); -Omega(2) Omega(1) 0];
      P_Omega = unit_Omega*unit_Omega';
      Xrel = ( P_Omega*(Xrel') + cos(norm(Omega)*dt).*(eye(3) - P_Omega)*(Xrel') + sin(norm(Omega)*dt).*(Omega_cross*(Xrel'))./norm(Omega) )';
      Urel = cross((Omega'.*ones(num_points,3)),(Xrel));
  %end
  
  %Update Xcm
  Xcm=Xcm +dt*Ucm;
  
  %Update Ucm
  F_normal_point = ((X(:,3)<0).*(S*(-X(:,3))-D*U(:,3)));
  F_normal_point = (F_normal_point > 0).*F_normal_point;
  
  Utan = [U(:,1),U(:,2),zeros(num_points,1)];
  Utan_norm = (vecnorm(Utan')' + 10^-10)*[1,1,1];
  Utan_unit = Utan./Utan_norm;
  F_friction = - mu*F_normal_point.*Utan_unit;
  %F_friction = (F_normal_point(:)>0).*(- mu*F_normal_point.*Utan_unit);
  %F_friction = ((X(:,3)<0).*(- mu*abs(M*g).*Utan_unit));
  %F_friction = - mu*abs(F_normal_point/S).*Utan_unit;
  %F_friction = f_friction_kit(X,mu,M_total,g,Utan_unit,num_points);
%   if ~(sum(X(:,3)<0)==0)
%       F_friction
%   end
  force_z = (-M)*g + F_normal_point;
  force_point = zeros(num_points,3);
  force_point(:,1) = F_friction(:,1);
  force_point(:,2) = F_friction(:,2);  
  force_point(:,3) = force_z;
  
  force = zeros(1,3);
   force(:,1) = sum(F_friction(:,1));
   force(:,2) = sum(F_friction(:,2));
   force(:,3) = sum(force_z);
  Ucm=Ucm+dt*force./M_total;
  
  
  %Update Angular Momentum
  
  torque = cross(Xrel,force_point);
  net_torque = sum(torque);
  L = L + dt*net_torque';
 

  tsave(clock)=clock*dt;
  %Zsave(clock)=Z(:,3);
  
  % update X
  X = Xcm + Xrel;
  U = Ucm + Urel;
  

  
  XcmSaveY(clock) = Xcm(2);
  XcmSaveZ(clock) = Xcm(3);
  [values,index] = min(X);
  Ymin(clock) = X(index(3),2);
  Xmin(clock) = X(index(3),1);
  UtanSave(:,clock) = Utan(index(3),1:2)';  
  
  energy_ground = 0;
  for i = 1:num_points
      if F_normal_point(i) > 0
          energy_ground = energy_ground + 0.5*S*X(i,3)^2;
      end
  end
  
  Esave(clock)= M_total*g*Xcm(3) + 0.5*M_total*norm(Ucm)^2 + 0.5*Omega'*I*Omega + energy_ground; 
  
  
  
  %plot
%   if mod(clock,200) == 0
%   figure(1)
%   title('Plot')
%     x = [X([1:end,1],1)];
%     y = [X([1:end,1],2)];
%     z = [X([1:end,1],3)];
%     plot3(x',y',z','linewidth',4,'Color','b')
%     hold on
%     plot3(X(1,1),X(1,2),X(1,3),'ro','MarkerSize',4);
%     plot3(X((num_points/2+1),1),X((num_points/2+1),2),X((num_points/2+1),3),'go','MarkerSize',4);
%     [value,index] = min(X);
%     
%     
%     
%     line([-10 10],[-10 -10],[-0 -0])
%     line([-10 10],[10 10],[-0 -0])
%     line([10 10],[10 -10],[-0 -0])
%     line([-10 -10],[-10 10],[-0 -0])
%     
%     for i = -9:9
%         line([-10 10],[-i -i],[-0 -0])
%         line([-i -i],[10 -10],[-0 -0])
%     end
%     view(45,15)
%     hold off 
%     
%     xlim([-10 10])
%     ylim([-10 10])
%     zlim([-10 10])
% 
%    %test for rigid body 
%    %norm(X(16,:)-X(1,:))
%    frame = getframe(gcf);
%    writeVideo(v,frame);
%   end

end

close(v)

figure(2)
%title('Top Point Z component axis vs Time')
plot(Xmin,Ymin)
figure(3)
%title('Center of Mass Z component axis vs Time')
plot(tsave,XcmSaveZ)
figure(4)
%title(Center of Mass Y component axis vs Time')
plot(tsave,XcmSaveY);

figure(5)
%title('Energy vs Time')
plot(tsave,Esave)
axis([0,tmax,0,12])
UtanVec = UtanSave(1,:).^2 + UtanSave(2,:).^2

figure(6)
plot(tsave,UtanVec)