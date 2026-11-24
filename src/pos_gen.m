function pos=pos_gen(grid,dir)
%POS_GEN Extracted unchanged from stat_script.m and suddeninput.m.
pos(:,1)=grid'*dir';
pos(:,1)=pos(:,1)-min(pos(:,1))+0.01;
pos(:,2)=1;
pos(2:2:end,2)=-1;
end
