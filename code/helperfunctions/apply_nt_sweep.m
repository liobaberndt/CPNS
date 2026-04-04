function Mout = apply_nt_sweep(Morig, param_list, adjustment, mode)
    Mout = Morig;
    for p = 1:numel(param_list)
        param = param_list{p};
        i = param.i;
        j = param.j;
        if strcmp(mode, 'abs')
            Mout(i, j) = Morig(i, j) + adjustment;
        else
            Mout(i, j) = Morig(i, j) + abs(Morig(i, j)) * adjustment;
        end
    end
end
