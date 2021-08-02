% Generate gamuts for each of the engineering paper examples

% Lamp -- 3 performance metrics, 4d gamut -- takes about 8 hours on avg
% [stats_l, buffer_l, des_figh_l, perf_figh_l, pts_l, labels_l, filledPtsArray_l, idxToVis_l] = appDiscover(29);
% view_4D_gamut(buffer_l, pts_l);

% % Turbine
% [stats_t, buffer_t, des_figh_t, perf_figh_t, pts_t, labels_t, filledPtsArray_t, idxToVis_t] = appDiscover(31);
% 
% % Bike Rocker 
% [stats_r, buffer_r, des_figh_r, perf_figh_r, pts_r, labels_r, filledPtsArray_r, idxToVis_r] = appDiscover(32);
% 
% % Bike Seat Stay
[stats_s, buffer_s, des_figh_s, perf_figh_s, pts_s, labels_s, filledPtsArray_s, idxToVis_s] = appDiscover(33);

% Bicopter
% [stats_b, buffer_b, des_figh_b, perf_figh_b, pts_b, labels_b, filledPtsArray_b, idxToVis_b] = appDiscover(11);

% Bicopter 2 contexts -- 4D gamut -- takes about 18-20 hours on avg
% [stats_b2, buffer_b2, des_figh_b2, perf_figh_b2, pts_b2, labels_b2, filledPtsArray_b2, idxToVis_b2] = appDiscover(12);

% Solar roofing -- about 20 minutes
% [stats_g100, buffer_g100, des_figh_g100, perf_figh_g100, pts_g100, labels_g100, filledPtsArray_g100, idxToVis_g100] = appDiscover(37);

% Lamp -- 3 performance metrics, 4d gamut
% [stats_l, buffer_l, des_figh_l, perf_figh_l, pts_l, labels_l, filledPtsArray_l, idxToVis_l] = appDiscover(29);
% view_4D_gamut(buffer_l, pts_l);
