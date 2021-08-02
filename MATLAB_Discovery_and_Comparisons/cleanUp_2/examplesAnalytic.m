% we should pull out the parameters here so they're all easily accessible

% ZDT1
[stats1, buffer1, des_figh1, perf_figh1, pts1, labels1, filledPtsArray1, idxToVis1] = appDiscover(1);

% ZDT2
[stats2, buffer2, des_figh2, perf_figh2, pts2, labels2, filledPtsArray2, idxToVis2] = appDiscover(2);

% ZDT3
[stats3, buffer3, des_figh3, perf_figh3, pts3, labels3, filledPtsArray3, idxToVis3] = appDiscover(3);

% ZDT4
[stats4, buffer4, des_figh4, perf_figh4, pts4, labels4, filledPtsArray4, idxToVis4] = appDiscover(4);

% ZDT6
% Note: ZDT6 uses problem ID 5! (since we do not use the discrete valued ZDT5 problem)
[stats6, buffer6, des_figh6, perf_figh6, pts6, labels6, filledPtsArray6, idxToVis6] = appDiscover(5);

% Fourier test 
[statsf, bufferf, des_fighf, perf_fighf, ptsf, labelsf, filledPtsArrayf, idxToVisf] = appDiscover(10);
