
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
if Param.StimProtocols.DriftingGrating7
    DG7.tempo=reshape(DG7.tempo,DG7.n_reps,DG7.directions);
    save(stim_fname, '-append','Param','DG7');
end
if Param.StimProtocols.DriftingGrating8
    DG8.tempo=reshape(DG8.tempo,DG8.n_reps,DG8.directions);
    save(stim_fname, '-append','Param','DG8');
end
if Param.StimProtocols.DriftingGrating9
    DG9.tempo=reshape(DG9.tempo,DG9.n_reps,DG9.directions);
    save(stim_fname, '-append','Param','DG9');
end



% disp(['* StimInfo file saved in ' save_dir ]);

