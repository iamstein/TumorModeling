function[] = AddR_AFI_TFI()

a = get(gcf,'Children');

str_old  = {'AFI' ,'TFI','Tss'};
str_new = {'AFIR','TFIR','Ttotss'}; 

for i_str = 1:length(str_old)
for i_axis=1:length(a)
    s_old = str_old{i_str};
    s_new = str_new{i_str};
    
    ai = a(i_axis);
    if isa(ai,'matlab.graphics.axis.Axes')
        tax = get(ai,'Title');
        tstr= get(tax,'String');
        tstr= strrep(tstr,s_old,s_new);
        set(tax,'String',tstr);
        
        yax = get(ai,'YLabel');
        ystr= get(yax,'String');
        ystr= strrep(ystr,s_old,s_new);
        ystr= strrep(ystr,' ratio','');
        set(yax,'String',ystr);
    end    
    1;
end
end
