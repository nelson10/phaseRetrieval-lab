function [fieldvec,Obf2F,Phase] = PhaseLift(fieldvec,mask,typePropagation,Obj1,wvl, delta1, deltan, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forward propagation
OP = propagation(typePropagation,Obj1,wvl, delta1, deltan, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Obj1 = abs(OP).^2;
Phasei = angle(OP);
In = Obj1;
Obj2 = sqrt(complex(In)).*exp(1i.*Phasei); %propagador schmidt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Backward propagation
backwardP = propagation(typePropagation,Obj2,wvl, delta1, deltan, -z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InverMn = ~mask;
Obf2F = backwardP.*mask+ fieldvec.*InverMn;
fieldvec = Obf2F;
Phase = angle(Obf2F);
end