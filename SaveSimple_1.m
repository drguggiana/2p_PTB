
%% Save stim parameters

if Param.StimProtocols.DriftingGrating
    DG.tempo=reshape(DG.tempo,DG.n_reps,DG.directions);
    save(stim_fname, '-append','Param','DG');
end

% disp(['* StimInfo file saved in ' save_dir ]);

