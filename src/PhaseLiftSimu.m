function [fieldvec,Obf2F,Phase] = PhaseLiftSimu(fieldvec,mask,typePropagation,Obj1,wvl, delta1, deltan, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = mask.*fieldvec;
%% Forward propagation
OP = propagation(typePropagation,J,wvl, delta1, deltan, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In = abs(OP).^2;
In = Obj1;
Phasei = angle(OP);
Obj2 = sqrt(complex(In)).*exp(1i.*Phasei); %propagador schmidt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Backward propagation
backwardP = propagation(typePropagation,Obj2,wvl, delta1, deltan, -z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InverMn = ~mask;
Obf2F = backwardP.*mask+ fieldvec.*InverMn;
Phase = angle(Obf2F);
fieldvec = Obf2F;
end