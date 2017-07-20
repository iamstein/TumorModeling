function[pout] = micro2macro_oral1_1cmt(p)
%converts pk parameters (CL,V2,ka) to equation parameters (C,ke,ka)
%equation with absorption:

    %deal with disambiguation
        f  = fieldnames(p);
        ii = find(ismember(f,{'ka','KA','Ka'}),1,'first');
        ka  = p.(f{ii});

        ii = find(ismember(f,{'alag','tlag','Tlag'}),1,'first');
        if isempty(ii)
            tlag = 0;
        else
            tlag  = p.(f{ii});
        end         

        ii= find(ismember(f,{'V1','V2','VF','V','Vc'}),1,'first');
        V2= p.(f{ii});   
        
        ii= find(ismember(f,{'CL','Cl','CLF'}),1,'first');
        if ~isempty(ii)
            CL= p.(f{ii});        
        else
            CL = V2*p.ke;
        end
                
        if isfield(p,'FD')
            FD  = p.FD;
        else
            FD = 1;
        end

    %convert
        pout.ka = ka;
        pout.ke = CL/V2;
        pout.C  = ka*FD/(V2*(ka-pout.ke));
        pout.FD = FD;
        pout.tlag = tlag;