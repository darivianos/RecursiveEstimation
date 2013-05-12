function [posEst,oriEst,radiusEst, posVar,oriVar,radiusVar,estState] = Estimator(estState,actuate,sense,tm,knownConst,designPart)
% [posEst,oriEst,posVar,oriVar,baseEst,baseVar,estState] =
% 	Estimator(estState,actuate,sense,tm,knownConst,designPart)
%
% The estimator.
%
% The Estimator function shall be used for both estimator design parts; the
% input argument designPart is used to distinguish the two:
%   designPart==1  -> Part 1
%   designPart==2  -> Part 2
%
% The function will be called in two different modes:
% If tm==0, the estimator is initialized; otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_v(k-1), drive wheel angular velocity
%                   actuate(2): u_r(k-1), drive wheel angle
%   sense           sensor measurements z(k), [1x2]-vector, INF if no
%                   measurement
%                   sense(1): z_d(k), distance measurement
%                   sense(2): z_r(k), orientation measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   knownConst      known constants (from KnownConstants.m)
%   designPart      variable to distinguish the estimator design part
%                       designPart==1  -> Part 1
%                       designPart==2  -> Part 2
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): x position estimate
%                   posEst(2): y position estimate
%   oriEst          orientation estimate (time step k), scalar
%   radiusEst       estimate of wheel radius W (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   radiusVar       variance of wheel radius estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2013
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Markus Hehn
% hehnm@ethz.ch
%
% --
% Revision history
% [19.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    adapted version for spring 2012, added unknown wheel
%                   radius
% [06.05.13, MH]    2013 version



if (designPart == 1)
    % Mode 1: Initialization
    if (tm == 0)
        % Do the initialization of your estimator here!
        % We initialize the initial nominal values
        posEst = [0 0];
        oriEst = 0;
        posxVar = knownConst.TranslationStartBound * knownConst.TranslationStartBound / 3;
        posyVar = posxVar;
        posVar = [posxVar posyVar];
        oriVar = knownConst.RotationStartBound * knownConst.RotationStartBound / 6;
        radiusEst = knownConst.NominalWheelRadius;
        radiusVar = knownConst.WheelRadiusError * knownConst.WheelRadiusError / 3;
        Pinit = [posVar(1),0,0,0;...
            0,posVar(2),0,0;...
            0,0,oriVar,0;...
            0,0,0,radiusVar];
        % Keep the current state to be provided as an input to next iteration

        estState.tm = tm;
        estState.posEst = posEst;
        estState.oriEst = oriEst;
        estState.radiusEst = radiusEst;
        estState.P = Pinit;
        return;
    end
    
    
    % Mode 2: Estimator iteration.
    
    % Prediction step, solve the three differential equations using ode45
    Ts = tm - estState.tm;
    % previous states
    xprev = estState.posEst(1);
    yprev = estState.posEst(2);
    rprev = estState.oriEst;
    wprev = estState.radiusEst;
    
    % input Signals
    B = knownConst.WheelBase;
    uv = actuate(1);
    ur = actuate(2);
    
    % Differential equation for prediction of the states
    [~,Y] = ode45(@(t,y) propState(t,y,uv,ur,B),[0 Ts],[xprev,yprev,rprev,wprev]);
    outSize = size(Y,1);
    xpred = Y(outSize,1);
    ypred = Y(outSize,2);
    rpred = Y(outSize,3);
    wpred = Y(outSize,4);
    
    % Differential equation for prediction of Variances
    A = [0, 0, -wprev*uv*cos(ur)*sin(rprev), uv*cos(ur)*cos(rprev);...
        0, 0, wprev*uv*cos(ur)*cos(rprev), uv*cos(ur)*sin(rprev);...
        0, 0, 0,  -uv * sin(ur)/B;...
        0, 0, 0, 0];
    
    % Process Noise Characteristics
    Q = zeros(4,4);
    Q(1,1) = .5;
    Q(2,2) = Q(1,1);
    Q(3,3) = .01;
    Q(4,4) = 1e-6;
    
    Pinit = estState.P;
  
    [~,Yvar] = ode45(@(t,y) propVariance(t,y,A,Q),[0 Ts],Pinit);
    outSize = size(Yvar,1);
    Ppred = [Yvar(outSize,1),Yvar(outSize,5),Yvar(outSize,8),Yvar(outSize,13);...
            Yvar(outSize,2),Yvar(outSize,6),Yvar(outSize,10),Yvar(outSize,14);...
            Yvar(outSize,3),Yvar(outSize,7),Yvar(outSize,11),Yvar(outSize,15);...
            Yvar(outSize,4),Yvar(outSize,8),Yvar(outSize,12),Yvar(outSize,16)];
    
    % Update/Measurement Step
    H = [xpred/sqrt(xpred*xpred + ypred*ypred), ypred/sqrt(xpred*xpred + ypred*ypred), 0, 0;...
        0,0,1,0];
    
    R = [knownConst.DistNoise*knownConst.DistNoise/6, 0;...
        0, knownConst.CompassNoise];
    
    K = Ppred * H' /(H * Ppred * H' + R);
    Pmeas = (eye(4) - K*H) * Ppred;
    if sense(1) == inf && sense(2) == inf
        Pmeas = Ppred;
        sense(1) = sqrt(xpred*xpred + ypred*ypred);
        sense(2) = rpred;
    end
    if sense(1) == inf
        for i = 1:2
            for j = 1:2
                Pmeas(i,j) = Ppred(i,j);
            end
        end
        sense(1) = sqrt(xpred*xpred + ypred*ypred);
    end
    if sense(2) == inf
        Pmeas(3,3) = Ppred(3,3);
        sense(2) = rpred;
    end
    
    e = [sense(1) - sqrt(xpred*xpred + ypred*ypred);...
        sense(2) - rpred];
    temp = K * e;
    
    posEst(1) = xpred + temp(1);
    posEst(2) = ypred + temp(2);
    oriEst = rpred + temp(3);
    radiusEst = wpred + temp(4);
   
    posVar(1) = Pmeas(1,1);
    posVar(2) = Pmeas(2,2);
    oriVar = Pmeas(3,3);
    radiusVar = Pmeas(4,4);
    
    % Keep the current state to be provided as an input to next iteration
    estState.tm = tm;
    estState.posEst = posEst;
    estState.oriEst = oriEst;
    estState.radiusEst = radiusEst;
    estState.P = Pmeas;
end

return