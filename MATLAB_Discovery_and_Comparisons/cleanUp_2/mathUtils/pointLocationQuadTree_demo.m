%% Quadtree point search in a triangulation
% Given a triangulation, we want to find the indices of the triangles
% containing some query points, along with their respective barycentric
% coordinates. We compare the computation times for the pointLocation and
% pointLocationQuadTree algorithms depending on the size of the search
% problem, setting the number query points equal to the number of
% triangulation points.

%% Setup

close all
clear

nPoints = 4.^(1:7);

%% Start computation loop

for i=1:length(nPoints)

  %% Test data
  %
  % Delaunay triangulation of a random point set.

  delaunayTri = delaunayTriangulation(rand(nPoints(i),2));
  nTriangles(i) = size(delaunayTri,1);

  %%
  % Standard triangulation of the same set, neglecting the Delaunay
  % property.

  standardTri = triangulation(delaunayTri.ConnectivityList, ...
                              delaunayTri.Points);

  %%
  % Random set of query points.

  queryPoints = rand(nPoints(i),2);

  %% Matlab search in a Delaunay triangulation
  %
  % The Matlab search in a Delaunay triangulation is very fast. We will not
  % be able to beat it.

  tic
  [indMatlabDelaunay,barMatlabDelaunay] = pointLocation(delaunayTri,queryPoints);
  tictocMatlabDelaunay(i) = toc; %#ok<*SAGROW>

  %% Matlab search in a standard triangulation
  %
  % As of Matlab R2018a, the Matlab search in a standard triangulation is
  % quite slow for large triangulations and a large numbers of query
  % points. (Think of interpolating finite element solutions.)

  tic
  [indMatlabStandard,barMatlabStandard] = pointLocation(standardTri,queryPoints);
  tictocMatlabStandard(i) = toc;
  
  %% Quadtree search in a Delaunay triangulation
  % Quadtree search works with Delaunay triangulations, but directly calls
  % the standard Matlab search in a Delaunay triangulation.

  tic
  [indQuadtreeDelaunay,barQuadtreeDelaunay] = pointLocationQuadTree(delaunayTri,queryPoints);
  tictocQuadtreeDelaunay(i) = toc;

  %% Quadtree search in a standard triangulation
  % Quadtree search in a standard triangulation is much faster than Matlab
  % search in a standard triangulation if the triangulation and number of
  % query points are large.

  tic
  [indQuadtreeStandard,barQuadtreeStandard] = pointLocationQuadTree(standardTri,queryPoints);
  tictocQuadtreeStandard(i) = toc;
  
  %% Test for correct results
  % Check that all methods find exactly the same triangle indices. The case
  % that a random query point is exactly on the edge between two triangles
  % is rare and, therefore, not taken care of.

  assert(isequaln(indMatlabDelaunay,indMatlabStandard,indQuadtreeDelaunay,indQuadtreeStandard))

  %% 
  % The barycentric coordinates are subject to round-off, so we can not
  % expect exact equality of the results. Checking against an absolute
  % tolerance is sufficient for this test.

  tolerance = 1e-12;
  pick = ~isnan(indMatlabDelaunay);
  if any(pick)  % assert() does not work for empty inputs
    approximationError = @(x)max(max(abs(barMatlabDelaunay(pick,:)-x(pick,:))));
    assert(approximationError(barMatlabStandard)<tolerance);
    assert(approximationError(barQuadtreeDelaunay)<tolerance);
    assert(approximationError(barQuadtreeStandard)<tolerance);
  end
  
%% Terminate computation loop

end

%% Display results
% The slopes suggest that the cost of the Matlab search in a standard
% triangulation depends about quadratically on the problem size, while the
% Quadtree search in a standard triangulation depends about linearly on the
% problem size. The vertical dotted line indicates the point at which
% pointLocationQuadTree actually builds a quadtree. For smaller sizes,
% pointLocationQuadTree directly calls pointLocation.

clf
loglog(nPoints,tictocMatlabStandard,'-s','DisplayName','Matlab (standard)')
hold on
loglog(nPoints,tictocQuadtreeStandard,'-x','DisplayName','Quadtree (standard)')
loglog(nPoints,tictocQuadtreeDelaunay,'-o','DisplayName','Quadtree (Delaunay)')
loglog(nPoints,tictocMatlabDelaunay,'-^','DisplayName','Matlab (Delaunay)')
legend('Location','NorthWest')
xlabel('n = number of query points = number of mesh points')
ylabel('computation time in seconds')
title('quadtree vs. standard search')

threshold = find(nPoints.*nTriangles >= 1e6,1,'first');
pick = threshold:length(nPoints);

if length(pick)>1

  quadraticSlope = 2*tictocMatlabStandard(end)*(nPoints(pick)/nPoints(end)).^2;
  loglog(nPoints(pick),quadraticSlope,'k--','DisplayName','O(n^2)')

  linearSlope = 0.5*tictocQuadtreeStandard(end)*(nPoints(pick)/nPoints(end)).^1;
  loglog(nPoints(pick),linearSlope,'k-','DisplayName','O(n)')

  loglog(nPoints(threshold)*[1 1],get(gca,'YLim'),'k:','DisplayName','Quadtree threshold');

end
