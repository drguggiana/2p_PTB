
%% Save stim parameters

if Param.StimProtocols.DriftingGrating
    DG.tempo=reshape(DG.tempo,DG.n_reps,DG.directions);
    save(stim_fname, '-append','Param','DG');
end
if Param.StimProtocols.DriftingGrating2
    DG2.tempo=reshape(DG2.tempo,DG2.n_reps,DG2.directions);
    save(stim_fname, '-append','Param','DG2');
end
if Param.StimProtocols.DriftingGrating3
    DG3.tempo=reshape(DG3.tempo,DG3.n_reps,DG3.directions);
    save(stim_fname, '-append','Param','DG3');
end
if Param.StimProtocols.DriftingGrating4
    DG4.tempo=reshape(DG4.tempo,DG4.n_reps,DG4.directions);
    save(stim_fname, '-append','Param','DG4');
end
if Param.StimProtocols.DriftingGrating5
    DG5.tempo=reshape(DG5.tempo,DG5.n_reps,DG5.directions);
    save(stim_fname, '-append','Param','DG5');
end
if Param.StimProtocols.DriftingGrating6
    DG6.tempo=reshape(DG6.tempo,DG6.n_reps,DG6.directions);
    save(stim_fname, '-append','Param','DG6');
end
if Param.StimProtocols.DriftingGratingDisparity
    DGD.tempo=reshape(DGD.tempo,DGD.n_reps,DGD.directions*DGD.nPhases);
    save(stim_fname, '-append','Param','DGD');
end
if Param.StimProtocols.DriftingGratingDisparity2
    DGD2.tempo=reshape(DGD2.tempo,DGD2.n_reps,DGD2.directions*DGD2.nPhases);
    save(stim_fname, '-append','Param','DGD2');
end
if Param.StimProtocols.DriftingGratingDisparity3
    DGD3.tempo=reshape(DGD3.tempo,DGD3.n_reps,DGD3.directions*DGD3.nPhases);
    save(stim_fname, '-append','Param','DGD3');
end
if Param.StimProtocols.DriftingGratingFlow
    DGF.tempo=reshape(DGF.tempo,DGF.n_reps,DGF.directions*DGF.nPhases);
    save(stim_fname, '-append','Param','DGF');
end
if Param.StimProtocols.DriftingGratingFlow2
    DGF2.tempo=reshape(DGF2.tempo,DGF2.n_reps,DGF2.directions*DGF2.nPhases);
    save(stim_fname, '-append','Param','DGF2');
end
if Param.StimProtocols.DriftingGratingFlow3
    DGF3.tempo=reshape(DGF3.tempo,DGF3.n_reps,DGF3.directions*DGF3.nPhases);
    save(stim_fname, '-append','Param','DGF3');
end
if Param.StimProtocols.DriftingGratingInDepth
    DGID.tempo=reshape(DGID.tempo,DGID.n_reps,DGID.nConditions_tot);
    save(stim_fname, '-append','Param','DGID');
end
if Param.StimProtocols.DriftingGratingPhi
    DGP.tempo=reshape(DGP.tempo,DGP.directions*length(DGP.phase),DGP.n_reps);
    save(stim_fname, '-append','Param','DGP');
end
if Param.StimProtocols.DriftingGratingPhi
    SG.tempo=reshape(SG.tempo,SG.directions*length(SG.phase),SG.n_reps);
    save(stim_fname, '-append','Param','SG');
end
if Param.StimProtocols.RetinotopicPosition
    RP.tempo=reshape(RP.tempo,RP.n_reps,prod(RP.positions));
    save(stim_fname, '-append','Param','RP');
    if Param.Dichoptic
        RP2.tempo=reshape(RP2.tempo,RP2.n_reps,prod(RP2.positions));
        save(stim_fname, '-append','Param','RP2');
    end
end
if Param.StimProtocols.RetinotopicPositionStripes
    RPS.tempo=reshape(RPS.tempo,RPS.n_reps,sum(RPS.positions));
    save(stim_fname, '-append','Param','RPS');
    if Param.Dichoptic && exist('RPS2','var')
        RPS2.tempo=reshape(RPS2.tempo,RPS2.n_reps,sum(RPS2.positions));
        save(stim_fname, '-append','Param','RPS2');
    end
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
if ismember(2,Param.StimProtocols.RFMapping)
    save(stim_fname, '-append','Param','RFM2');
end
if Param.StimProtocols.PhiMotionGratings
    PMG.tempo=reshape(PMG.tempo,PMG.orientations*PMG.nCombos,PMG.n_reps);
    save(stim_fname, '-append','Param','PMG');
end
if Param.StimProtocols.DriftingGratingDisparity_phi
    DGDphi.tempo=reshape(DGDphi.tempo,DGDphi.n_reps,DGDphi.directions*DGDphi.nPhases*2);
    save(stim_fname, '-append','Param','DGDphi');
end
if Param.StimProtocols.PhiMotionGratingsControl
    PMGc1.tempo=reshape(PMGc1.tempo,PMGc1.directions*length(PMGc1.phase),PMGc1.n_reps);
    save(stim_fname, '-append','Param','PMGc1');
    PMGc2.tempo=reshape(PMGc2.tempo,PMGc2.directions*length(PMGc2.phase),PMGc2.n_reps);
    save(stim_fname, '-append','Param','PMGc2');
end
if Param.StimProtocols.PhiGratings
    PG.tempo=reshape(PG.tempo,PG.nbr_stimuli,PG.n_reps);
    save(stim_fname, '-append','Param','PG');
end
if Param.StimProtocols.PhiMotionBars
    PMB.tempo=reshape(PMB.tempo,PMB.nbr_stimuli,PMB.n_reps);
    save(stim_fname, '-append','Param','PMB');
%     imshow(ScreenImg_merge);
end
if Param.StimProtocols.DotMotion
    DM.tempo=reshape(DM.tempo,DM.n_reps,DM.directions);
    save(stim_fname, '-append','Param','DM');
end
if Param.StimProtocols.DotMotion2
    DM2.tempo=reshape(DM2.tempo,DM2.n_reps,DM2.directions);
    save(stim_fname, '-append','Param','DM2');
end
if Param.StimProtocols.DotMotion3
    DM3.tempo=reshape(DM3.tempo,DM3.n_reps,DM3.directions);
    save(stim_fname, '-append','Param','DM3');
end
if Param.StimProtocols.RandomDotStereogram
    RDS.tempo=reshape(RDS.tempo,RDS.n_reps,RDS.orientations*RDS.n_disparity_intervals);
    save(stim_fname, '-append','Param','RDS');
end
if Param.StimProtocols.RandomDotStereogram2
    RDS2.tempo=reshape(RDS2.tempo,RDS2.n_reps,RDS2.orientations*RDS2.n_disparity_intervals);
    save(stim_fname, '-append','Param','RDS2');
end
if Param.StimProtocols.RandomBarStereogram
    try
    RBS.tempo=reshape(RBS.tempo,RBS.n_reps,RBS.orientations*RBS.n_disparity_intervals);
    end
    save(stim_fname, '-append','Param','RBS');
end
if Param.StimProtocols.RandomBarStereogram2
    try
    RBS2.tempo=reshape(RBS2.tempo,RBS2.n_reps,RBS2.orientations*RBS2.n_disparity_intervals);
    end
    save(stim_fname, '-append','Param','RBS2');
end
if Param.StimProtocols.NaturalImageStereogram
    NIS.tempo=reshape(NIS.tempo,NIS.n_reps,NIS.orientations*NIS.n_disparity_intervals);
    NIS.seqImages=reshape(NIS.seqImages,NIS.n_reps,NIS.orientations*NIS.n_disparity_intervals);
    save(stim_fname, '-append','Param','NIS');
end
if Param.StimProtocols.RandomDotStereoBar
    RDSB.mask_bar = [];
    try
    RDSB.tempo=reshape(RDSB.tempo,RDSB.n_reps,RDSB.directions*RDSB.n_disparities);
    end
    save(stim_fname, '-append','Param','RDSB');
end
if Param.StimProtocols.RandomDotStereoEdge
    RDSE.mask_edge = []; RDSE.mask_gap = [];
    try
    RDSE.tempo=reshape(RDSE.tempo,RDSE.n_reps,RDSE.directions*RDSE.n_disparities);
    end
    save(stim_fname, '-append','Param','RDSE');
end

% disp(['* StimInfo file saved in ' save_dir ]);

