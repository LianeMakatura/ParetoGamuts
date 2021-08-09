filename = "norcoART.jpg";

if false
    imshow(filename);
    title("Please select points in the order shown in the console.");

    disp("Order of points:")
    disp("1: Rear Contact patch (where rear tire meets the ground)")
    disp("2: Front contact patch (where front tire meets the ground)")
    disp("3: Rear Axle")
    disp("4: Front Axle")
    disp("5: Top of Rear Gear (that the chain goes around)")
    disp("6: Top of Front Gear (that the chain goes around, near the pedals)")
    disp("7: Linkage pt -- main pivot (where rear suspension linkage is attached to main frame near the pedals)")
    disp("8: Linkage pt -- chainstay pivot (pivot between linkage bars, near the rear axle)")
    disp("9: Linkage pt -- seat tube pivot (where the rocker is fixed to the near-vertical seat tube on the main frame)")
    disp("10: Linkage pt -- pivot between rocker and the seat stay")
    disp("11: Linkage pt -- fixed shock pivot (where the shock is fixed to the main frame)")
    disp("12: Linkage pt -- pivot between rocker and shock")
    disp("13: Guess for average rider's center of mass above seat")
    disp("14-18: 5 points along the upper seat tube")
    disp("19-23: 5 points along the top bar")
    disp("24-28: 5 points along the lower diagonal bar")
    disp("29-33: 5 points along the front steering tube")
    disp("34-36: 3 points at seat corners")
    disp("37-39: 3 points along handlebars")
    [x, y] = getpts;
    
    orig_pixel_pts = struct('x', x, 'y', y);
    save('rawBikeMeasurements.mat', 'orig_pixel_pts');

else
    raw = load('rawBikeMeasurements.mat');
    orig_pixel_pts = raw.orig_pixel_pts;
    x = orig_pixel_pts.x;
    y = orig_pixel_pts.y;
end

% correct the values that should line up exactly
% rear contact/axle x should match
newRearX = (x(1) + x(3)) / 2;
x(1) = newRearX;
x(3) = newRearX;

% front contact/axle x should match
newFrontX = (x(2) + x(4)) / 2;
x(2) = newFrontX;
x(4) = newFrontX;

% rear/front contact y should match 
newGroundY = (y(1) + y(2)) / 2;
y(1) = newGroundY;
y(2) = newGroundY;

% rear/front axle y should match 
newAxleY = (y(3) + y(4)) / 2;
y(3) = newAxleY;
y(4) = newAxleY;

% offset so coordinate system has 0,0 at rear contact patch
x = x - x(1);
y = -(y - y(1));

% get scale factor 
knownWheelBase = 1211; % mm
knownWheelHeight = 698.5; %mm, 27.5"
knownAxleHeight = knownWheelHeight / 2;

imgWheelBase = norm([x(2) - x(1), y(2)-y(1)]);
scaleFactorToReal = knownWheelBase / imgWheelBase;

scaledAxleHeight = y(3) * scaleFactorToReal;
fprintf("Scaled height is \t%0.6f\n Actual height is \t%0.6f\n", scaledAxleHeight, knownAxleHeight);


x = x*scaleFactorToReal;
y = y*scaleFactorToReal;

% place into struct
scaled_bike = struct(...
    'rearContact', [x(1), y(1)], ...
    'frontContact', [x(2), y(2)], ...
    'rearAxle', [x(3), y(3)], ...
    'frontAxle', [x(4), y(4)], ...
    'rearGearTop', [x(5), y(5)], ...
    'frontGearTop', [x(6), y(6)], ...
    'pedalPivot', [x(7), y(7)], ...
    'rearAxlePivot', [x(8), y(8)], ...
    'seatTubePivot', [x(9), y(9)], ...
    'rockerSeatStayPivot', [x(10), y(10)],...
    'fixedShockPt', [x(11), y(11)],...
    'rockerShockPt', [x(12), y(12)],...
    'riderCOM', [x(13), y(13)],...
    'frame_upperSeatTube', [[x(14), y(14)]; [x(15), y(15)]; [x(16), y(16)]; [x(17), y(17)]; [x(18), y(18)]],...
    'frame_topBar', [[x(19), y(19)]; [x(20), y(20)]; [x(21), y(21)]; [x(22), y(22)]; [x(23), y(23)]],...
    'frame_lowerBar', [[x(24), y(24)]; [x(25), y(25)]; [x(26), y(26)]; [x(27), y(27)]; [x(28), y(28)]],...
    'frame_frontSteerBar', [[x(29), y(29)]; [x(30), y(30)]; [x(31), y(31)]; [x(32), y(32)]; [x(33), y(33)]],...
    'frame_seat', [[x(34), y(34)]; [x(35), y(35)]; [x(36), y(36)]],...
    'frame_handlebars', [[x(37), y(37)]; [x(38), y(38)]; [x(39), y(39)]]...
);

bike = struct('original_pixel_measures', orig_pixel_pts, 'scaleFactorToReal', scaleFactorToReal, 'scaled_bike', scaled_bike);

save('fullBikeMeasurements.mat', 'bike');
