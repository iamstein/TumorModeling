function[h] = add_title_label(Position)

    a = get(gcf,'Children');
    set(gcf,'CurrentAxes',a(end));
    h = text(0,0,'hi');
    
    otxt.Position = Position;
    otxt.String = ['Adjusted Var:' char(10) ... 
                   'Param. Range:' char(10) ... 
                   'Derived Var:'];
    otxt.HorizontalAlignment = 'Left';
    otxt.FontSize = 11;
    %otxt.BackgroundColor = [1 1 1];
    %otxt.FontWeight = 'Bold';
    set(h,otxt);