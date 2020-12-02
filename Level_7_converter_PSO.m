clc

clear all

close all

ma=0.01:0.01:1.2;

for kk=1:length(ma)

    XX=linspace(0,pi/2,4);

    xmin=XX(1,1:3); % Upper Bound of Variables

    xmax=XX(1,2:4); % Lower Bound of Variables

    n=length(xmax);   % Number of Decision Variables

    N=10*n;              % Population Size (Swarm Size)

    maxit=9000;      % Maximum Number of Iterations

    % Velocity Limits

    VelMax=0.1*(xmax-xmin);

    VelMin=-VelMax;


    W=0.9;            % Inertia Weight


    c1=2;           % Personal Learning Coefficient


    c2=2;           % Global Learning Coefficient

    cost_t=0;

    for i=1:N

        particle(i,:)=rand(1,n).*(xmax-xmin)+xmin; % initial Particles

        x=particle(i,:);

        if x(1)<0.0175
           cost_t=1e10;
        end
        if x(1)>=x(2)
            cost_t=1e10;
        end
        if x(2)>=x(3)
            cost_t=1e10;
        end

        if x(3)>((pi/2)-0.0175)
            cost_t=1e10;
        end
        if x(2)-x(1)<0.0175
            cost_t=1e10;
        end
        if x(3)-x(2)<0.0175
            cost_t=1e10;
        end

        fitness(i,1)=cost_t+(100*((3*ma(kk))-(cos(x(1))+cos(x(2))+cos(x(3))))/(3*ma(kk))).^4 +...
                     (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2)+...
                     (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2); % Fitness Function
        cost_t=0;
        vel(i,:)=rand(1,n).*(VelMax-VelMin)+VelMin; % initial Velocity

    end

        %  Personal Best

        Pbest=particle;  

        fitness_Pbest=fitness;

        %  Global Best

        particle_mix=[particle,fitness];

        particle_sort=sortrows(particle_mix,n+1);

        Gbest=particle_sort(1,1:n);  

        fitness_Gbest=particle_sort(1,n+1:end);

        xmin=[0 0 0];
        xmax=[pi/2 pi/2 pi/2];

    for iter=1:maxit

        for i=1:N

            % Update Velocity
            vel(i,:)=W*vel(i,:)+...
                     c1*rand(1,n).*(Pbest(i,:)-particle(i,:))+...
                     c2*rand(1,n).*(Gbest(1,:)-particle(i,:));

            % Apply Velocity Limits
            vel(i,:)=min([VelMax;vel(i,:)]);
            vel(i,:)=max([VelMin;vel(i,:)]);


            % Update Position
            particle(i,:)=particle(i,:)+vel(i,:);


            % Apply Position Limits
            particle(i,:)=min([xmax;particle(i,:)]);
            particle(i,:)=max([xmin;particle(i,:)]);

            % Evaluation

            x=particle(i,:);

            if x(1)<0.0175
               cost_t=1e10;
            end
            if x(1)>=x(2)
                cost_t=1e10;
            end
            if x(2)>=x(3)
                cost_t=1e10;
            end

            if x(3)>((pi/2)-0.0175)
                cost_t=1e10;
            end
            if x(2)-x(1)<0.0175
                cost_t=1e10;
            end
            if x(3)-x(2)<0.0175
                cost_t=1e10;
            end


            fitness(i,1)=cost_t+(100*((3*ma(kk))-(cos(x(1))+cos(x(2))+cos(x(3))))/(3*ma(kk))).^4 +...
                             (1/5)*((50*(1/5)*(cos(5*x(1))+cos(5*x(2))+cos(5*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2)+...
                             (1/7)*((50*(1/7)*(cos(7*x(1))+cos(7*x(2))+cos(7*x(3)))/(cos(x(1))+cos(x(2))+cos(x(3)))).^2); % Fitness Function

            cost_t=0;


            % Update Personal Best
            if fitness(i,1)<fitness_Pbest(i,1)
               Pbest(i,:)=particle(i,:);
               fitness_Pbest(i,1)= fitness(i,1);
            end

            % Update Global Best
            if fitness(i,1)<fitness_Gbest(1,1)
               Gbest(1,:)=particle(i,:);
               fitness_Gbest(1,1)= fitness(i,1);
            end

        end

        final(iter,:)=[Gbest fitness_Gbest];

    end

    %% Results

    V1star=3*ma(kk);
    V1= cos(Gbest(1))+cos(Gbest(2))+cos(Gbest(3));
    V1nesbi=V1/V1star;
    V5=abs((1/5)*(cos(5*Gbest(1))+cos(5*Gbest(2))+cos(5*Gbest(3)))/V1);
    V7=abs((1/7)*(cos(7*Gbest(1))+cos(7*Gbest(2))+cos(7*Gbest(3)))/V1);
    V11=abs((1/11)*(cos(11*Gbest(1))+cos(11*Gbest(2))+cos(11*Gbest(3)))/V1);
    V13=abs((1/13)*(cos(13*Gbest(1))+cos(13*Gbest(2))+cos(13*Gbest(3)))/V1);
    V17=abs((1/17)*(cos(17*Gbest(1))+cos(17*Gbest(2))+cos(17*Gbest(3)))/V1);
    V19=abs((1/19)*(cos(19*Gbest(1))+cos(19*Gbest(2))+cos(19*Gbest(3)))/V1);
    V23=abs((1/23)*(cos(23*Gbest(1))+cos(23*Gbest(2))+cos(23*Gbest(3)))/V1);
    V25=abs((1/25)*(cos(25*Gbest(1))+cos(25*Gbest(2))+cos(25*Gbest(3)))/V1);
    V29=abs((1/29)*(cos(29*Gbest(1))+cos(29*Gbest(2))+cos(29*Gbest(3)))/V1);
    V31=abs((1/31)*(cos(31*Gbest(1))+cos(31*Gbest(2))+cos(31*Gbest(3)))/V1);
    V35=abs((1/35)*(cos(35*Gbest(1))+cos(35*Gbest(2))+cos(35*Gbest(3)))/V1);
    V37=abs((1/37)*(cos(37*Gbest(1))+cos(37*Gbest(2))+cos(37*Gbest(3)))/V1);
    V41=abs((1/41)*(cos(41*Gbest(1))+cos(41*Gbest(2))+cos(41*Gbest(3)))/V1);
    V43=abs((1/43)*(cos(43*Gbest(1))+cos(43*Gbest(2))+cos(43*Gbest(3)))/V1);
    V47=abs((1/47)*(cos(47*Gbest(1))+cos(47*Gbest(2))+cos(47*Gbest(3)))/V1);
    V49=abs((1/49)*(cos(49*Gbest(1))+cos(49*Gbest(2))+cos(49*Gbest(3)))/V1);

    THD=100.*sqrt((V5^2)+(V7^2)+(V11^2)+(V13^2)+(V17^2)+(V19^2)+(V23^2)+(V25^2)+(V29^2)+(V31^2)+(V35^2)+(V37^2)+(V41^2)+(V43^2)+(V47^2)+(V49^2));

    % t=linspace(0,2*pi,360);
    % y=3*sin(t);
    % 
    % th=[0  Gbest(1)  Gbest(1)  Gbest(2)   Gbest(2)   Gbest(3)     Gbest(3)      pi-Gbest(3)   pi-Gbest(3)   pi-Gbest(2)   pi-Gbest(2)   pi-Gbest(1)   pi-Gbest(1)  pi+Gbest(1)  pi+Gbest(1)     pi+Gbest(2)   pi+Gbest(2)  pi+Gbest(3)  pi+Gbest(3)   2*pi-Gbest(3)   2*pi-Gbest(3) 2*pi-Gbest(2)   2*pi-Gbest(2)   2*pi-Gbest(1)    2*pi-Gbest(1)  2*pi]/pi*180;
    % u= [0     0          1        1           2           2           3              3            2             2             1              1            0            0             -1             -1            -2           -2            -3            -3              -2           -2              -1                -1           0             0  ];
    % 
    % 
    % figure(1);
    % plot(th,u,'LineWidth',3)
    % hold on
    % plot(y,'r','LineWidth',2)
    % xlabel('0 \leq \Theta \leq 360')
    % ylabel('Output Voltage (VDC p.u.)')
    % xlim([0 360])
    % ylim([-3.1 3.1])
    % grid
    % 
    % figure(2);
    % plot(final(:,end),'LineWidth',2)
    % xlabel('Number of Iteration')
    % ylabel('Fitness value')
    % grid

    Switching_Angles=Gbest(1,1:3)/pi*180;
    fitness_Gbest;
    
    Objective_Function_Value(kk)=fitness_Gbest;
    
    Fundamental_Harmonic(kk)=V1nesbi*100;
    
    V5th(kk)=V5*100;
    V7th(kk)=V7*100;
   

end

figure(1);
    semilogy(ma,Objective_Function_Value,'k','LineWidth',2)
    legend('PSO','IPSO')
    xlabel('Modulation Index')
    ylabel('ObjectiveFunction Value')
    grid
    
    figure(2);
    plot(ma,Fundamental_Harmonic,'k','LineWidth',2)
    ylim([0 110])
    legend('PSO','IPSO')
    xlabel('Modulation Index')
    ylabel('Fundamental Harmonic')
    grid
    
    figure(3);
    plot(ma,V5th,'LineWidth',2)
    hold on
    plot(ma,V7th,'r','LineWidth',2)
    legend('5th','7th')
    xlabel('Modulation Index')
    ylabel('Harmonics')
    grid

    