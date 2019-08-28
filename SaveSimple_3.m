
%% Save stim parameters

if Param.StimProtocols.DriftingGrating
%     DG.tempo=reshape(DG.tempo,DG.n_reps,DG.directions);
    save(stim_fname, '-append','Param','DG');
end
if Param.StimProtocols.Looming
%     LO.tempo=reshape(LO.tempo,LO.n_reps,LO.directions);
    save(stim_fname, '-append','Param','LO');
end
if ismember(1,Param.StimProtocols.RFMapping)
%     RFM.tempo=reshape(RFM.tempo,DGP.directions*length(DGP.phase),DGP.n_reps);
    save(stim_fname, '-append','Param','RFM');
    for i=1:length(screen)
        mon_pos = get(0,'MonitorPositions');
        figure; set(gcf,'Position',[mon_pos(2,1),mon_pos(2,2),500,500]);
        imshow(RFM.ScreenImg_merge(:,:,:,i),'InitialMagnification','fit'); %#ok<*NODEF>
    end
end
% disp(['* StimInfo file saved in ' save_dir ]);

