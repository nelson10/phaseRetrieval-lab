function [OP] = propagation(typePropagation,Obj1,wvl, delta1, deltan, z)
if(typePropagation==0)
    [OP] = ang_spec_pro(Obj1, wvl, delta1, deltan, z);
elseif(typePropagation==1)
    temp = delta1;
    delta1 = deltan;
    deltan = temp;
    [OP] = ang_spec_multi_prop(Obj1, wvl, delta1, deltan, z);
elseif(typePropagation==2)
    OP = fresnel(Obj1,wvl,delta1,z);
elseif(typePropagation==3)
    %[OP] = ang_spec_multi_prop_vac2(Uin, wvl, delta1, deltan, z,opc);
end
end
