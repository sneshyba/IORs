function f=getfirstnum(line)
% Interprets the first non-blank field in a character string as a number

i0=find(line == ' ');
i1 = find(diff(i0)>1);
i2=i1+1;
i3=i0(i2);
f=str2num(line(i2:i3));
return
