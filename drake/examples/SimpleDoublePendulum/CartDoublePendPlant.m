classdef CartDoublePendPlant < Manipulator
  
  properties
    % parameters from Spong95 (except inertias are now relative to the
    % joint axes)
    mc = 5;
    l1 = 1; l2 = 1;  
    m1 = 1; m2 = 1;  
    g = 9.81;
    b1=.1;  b2=.1;
    I1=[]; I2=[]; % set in constructor
    xG;
    uG;
  end
  
  methods
    function obj = CartDoublePendPlant
      obj = obj@Manipulator(3,3);
      obj.I1 = obj.m1*obj.l1^2;
      obj.I2 = obj.m2*obj.l2^2;
      obj = setStateFrame(obj,CoordinateFrame('CartDoublePendState', 6,'x',{'theta1','theta2','cart','theta1_dot','theta2_dot','cart_dot'}));
      obj = obj.setOutputFrame(obj.getStateFrame);
      obj.xG = Point(obj.getStateFrame,[pi; 0; 0; 0; 0; 0]);
%       obj.xG = Point(obj.getStateFrame,[pi; 0; 0; 0; 0; 0]);
      obj.uG = Point(obj.getInputFrame,[0; 0; 0]);
    end
    
    function [H,C,B] = manipulatorDynamics(obj,q,qd)
      % keep it readable:
      mc = obj.mc;
      m1=obj.m1; m2=obj.m2; l1=obj.l1; l2=obj.l2; g=obj.g; b1=obj.b1; b2=obj.b2; I1=obj.I1; I2=obj.I2;
      m2l1l2 = m2*l1*l2;  % occurs often!
        
      c = cos(q(1:2,:));  s = sin(q(1:2,:));  s12 = sin(q(1,:)+q(2,:)); c12 = cos(q(1,:)+q(2,:));
      
      h12 = I2 + m2l1l2*c(2);
      H = [ I1 + I2 + m2*l1^2 + 2*m2l1l2*c(2), h12, m1*l1*c(1); h12, I2, m2*(l1*c(1)+l2*c12); m1*l1*c(1), m2*(l1*c(1)+l2*c12), mc+m1+m2 ];
      
      C = [ -2*m2l1l2*s(2)*qd(2), -m2l1l2*s(2)*qd(2), 0; m2l1l2*s(2)*qd(1), 0, 0; 0, 0, 0 ];
      G = g*[ m1*l1*s(1) + m2*(l1*s(1)+l2*s12); m2*l2*s12; 0 ];
            
      % accumate total C and add a damping term:
      C = C*qd + G + [b1;b2;0].*qd;

      B = [0 0 0; 0 0 0; 0 0 1];
      
%       [f,df,d2f,d3f]=dynamics(obj,t,x,u)
    end
    
    function x = getInitialState(obj)
%        x = 1*randn(4,1);
        x = .1*randn(6,1);
      
    end
    
%     function qdd = sodynamics(obj, t, q, qd, u)
%         % Implement second-order dynamics
%         qdd = (u - obj.m*obj.g*obj.lc*sin(q) - obj.b*qd)/obj.I;
%     end
%     
%     function [f,df,d2f,d3f]=dynamics(obj,t,x,u)
%       f=dynamics@SecondOrderSystem(obj,t,x,u);
%       if (nargout>1)
%         [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
%       end
%     end
    
    function c = balanceDoublePendLQRTree(obj, Q, R)
        options.num_branches = 15;
%         options.stabilize=true;
        options.xs=[1e-6;1e-6;1e-6;1e-6;1e-6;1e-6];
%         options.xs=[0;0;0;0;];
        options.Tslb=2;
        options.Tsub=8;
        options.degL1=4;
%         options.plot_basins=true;
        %param3=@()rand(6,1).*[2*pi;10;0;0;0;0;]-[pi;5;0;0;0;0;];
        %param3=@()rand(4,1).*[2*pi;1;0.001;0.001]-[pi;0.5;0.001;0.001];
        %param3=@()rand(4,1).*[2*pi;10;0.0;0.0]-[pi;5;0.0;0.0];
        %c=LQRTree.buildLQRTree(obj,obj.xG,obj.uG,@()rand(6,1).*[2*pi;10;0;0;0;0]-[pi;5;0;0;0;0],Q,R,options);
        c=LQRTree.buildLQRTree(obj,obj.xG,obj.uG,@()rand(6,1).*[0;0;0;0;0;0]-[0;0;0;0;0;0],Q,R,options);
    end
    
  end
  
  methods(Static)
    function run()  % runs the passive system
      pd = CartDoublePendPlant;
      pv = CartDoublePendVisualizer(pd);
      traj = simulate(pd,[0 5],randn(4,1));
%s       traj = simulate(pd,[0 5],randn(6,1));
      playback(pv,traj);
    end
  end
  
end