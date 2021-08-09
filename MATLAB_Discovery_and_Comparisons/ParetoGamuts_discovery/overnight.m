
% figure
% [buffer1, des_figh1, perf_figh1, pts1, labels1, filledPtsArray1, idxToVis1] = appDiscover(1);
% lowerEnvelopeCompare(buffer1, pts1, [0, 1], [0, 0.5]);
% 
% figure
% [buffer2, des_figh2, perf_figh2, pts2, labels2, filledPtsArray2, idxToVis2] = appDiscover(2);
% buffer2.visualizeAppBuff(idxToVis2, true, pts2, labels2);
% lowerEnvelopeCompare(buffer2, pts2, [0, 1], [0, 0.5]);

% % figure
% [buffer3, des_figh3, perf_figh3, pts3, labels3, filledPtsArray3, idxToVis3] = appDiscover(3);
% lowerEnvelopeCompare(buffer3, pts3, [0,1], [0.25, 0.75]);

% [buffer4, des_figh4, perf_figh4, pts4, labels4, filledPtsArray4, idxToVis4] = appDiscover(4);
% lowerEnvelopeCompare(buffer4, pts4, [0,1], [0.0, 0.5]);

% [buffer6, des_figh6, perf_figh6, pts6, labels6, filledPtsArray6, idxToVis6] = appDiscover(5); % ZDT6 is mapping function ID 5
% lowerEnvelopeCompare(buffer6, pts6, [0,1], [0.0, 0.5]);


% [buffer7, des_figh7, perf_figh7, pts7, labels7, filledPtsArray7, idxToVis7] = appDiscover(10); % Fourier test
% lowerEnvelopeCompare(buffer7, pts7, [0,1], [0.0, 0.5]);

% 
% figure
% [t_buffer, t_des_figh, t_perf_figh, t_pts, t_labels, t_filledPtsArray, t_idxToVis] = appDiscover(31); % turbine
% singleDesignVsGamut(t_buffer, t_pts, t_labels, 20, size(t_pts, 1));
% DESIGN POINTS FOR FIGURE
% despts = [0.517244622093831,4.13561279074443e-05,0.999917301678813;
%           0.884962532816126,0.110135100912814,0.398716488916166;
%           0.920827329900260,0.651861821012161,0.909234056856701;
%           0.661304101472483,0.111901028688005,0.000187967582329334;
%           0.946458470015279,4.46806145549155e-05,0.745354890021532;
%           0.283444622409020,4.13561279074523e-05,0.999917301678813
%           ];
% 
% singleDesignVsGamut(t_buffer, t_pts, t_labels, 100, size(despts, 1), despts);



% [l_buffer, l_des_figh, l_perf_figh, l_pts, l_labels, l_filledPtsArray, l_idxToVis] = appDiscover(29); % lamp
% singleDesignVsGamut(t_buffer, t_pts, t_labels, 100, 5)
% GENERATING THE LAMP GAMUT IMAGE
% ------ because I forot to resclae metric
% l_buffer = buffer_l;
% l_pts = pts_l;
% 
% lpoints_rescaled = l_pts;
% lpoints_rescaled(:, l_buffer.rD+l_buffer.rA+1) = lpoints_rescaled(:, l_buffer.rD+l_buffer.rA+1) / 3;
% l_perf1idx = l_buffer.rD + l_buffer.rA + 1;
% scatter3(lpoints_rescaled(:, l_perf1idx), lpoints_rescaled(:, l_perf1idx+1), lpoints_rescaled(:, l_perf1idx+2), 40, lpoints_rescaled(:, l_perf1idx-1));
% xlabel("Stability"); ylabel("Mass"); zlabel("Focal Illumination"); 
% hcb = colorbar;
% colormap jet;
% colorTitleHandle = get(hcb,'Title');
% titleString = 'Focal Height';
% set(colorTitleHandle ,'String',titleString);


% [r_buffer, r_des_figh, r_perf_figh, r_pts, r_labels, r_filledPtsArray, r_idxToVis] = appDiscover(32); % rocker
% 
% [s_buffer, s_des_figh, s_perf_figh, s_pts, s_labels, s_filledPtsArray, s_idxToVis] = appDiscover(33); % stay
% s_buffer.visualizeAppBuff(s_idxToVis, false, s_pts, s_labels);

