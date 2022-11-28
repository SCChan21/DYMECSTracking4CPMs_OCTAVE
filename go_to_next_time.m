function [month,day,hour,minute]=go_to_next_time(month,day,hour,minute,timestep) 
    minute=minute+timestep;
    if minute>=60
      hour=hour+1;
      minute=minute-60;
      if hour==24
        hour=0;
        day=day+1;
        if (month==2 && day>28) || ((month==4 || month==6 || month==9 || month == 11) && day>30) || day>31
          day=1;
          month=month+1;
        end
      end
    end;
return;
