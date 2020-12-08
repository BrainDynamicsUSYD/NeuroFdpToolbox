% visualize LFP in required order
function OrderedLFPVisualize(R,Vorder)
LFP = R.LFP.LFP_broad(Vorder,:);
seg_size = 1e4; % 1s
step_tot = R.step_tot;
seg_num = ceil(step_tot/seg_size);
for seg = 1:seg_num
    seg_ind = get_seg(step_tot, seg_size, seg);
    ShowLFP = LFP(:,seg_ind);
    imagesc(ShowLFP)
    next = input('\t Next figure?');
    delete(gcf);
end
end