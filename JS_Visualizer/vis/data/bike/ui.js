// this is the function that determines where most of the joints are.
// It finds the two intersection points of two bars, A and B, with joints
// jointA and jointB respectively. It returns one of two solutions.
// The solutions are determined by whether the intersectionNum parameter
// is 0, 1, 2, or 3.  Each number corresponds to high, low, left, or right solution,
// as shown by the variables below.
const high = 0
const low = 1
const left = 2
const right = 3
function circIntersection(jointA, jointB, lengthA, lengthB, intersectionNum) {
    const xPos_a = jointA[0];
    const yPos_a = jointA[1];
    const xPos_b = jointB[0];
    const yPos_b = jointB[1];
    const Lc = Math.sqrt((Math.pow(xPos_a - xPos_b, 2)) + (Math.pow(yPos_a - yPos_b, 2)));
    const bb = ((lengthB * lengthB) - (lengthA * lengthA) + (Lc * Lc)) / (Lc * 2);
    const h = Math.sqrt((((lengthB) * lengthB) - (bb * bb)));
    const Xp = xPos_b + ((bb * (xPos_a - xPos_b)) / Lc);
    const Yp = yPos_b + ((bb * (yPos_a - yPos_b)) / Lc);
    const Xsolution1 = Xp + ((h * (yPos_b - yPos_a)) / Lc);
    const Ysolution1 = Yp - ((h * (xPos_b - xPos_a)) / Lc);
    const Xsolution2 = Xp - ((h * (yPos_b - yPos_a)) / Lc);
    const Ysolution2 = Yp + ((h * (xPos_b - xPos_a)) / Lc);
    const solution1 = [Xsolution1, Ysolution1];
    const solution2 = [Xsolution2, Ysolution2];
    if (intersectionNum == 0)
        if (Ysolution1 > Ysolution2)
            return solution1;
        else return solution2;
    else if(intersectionNum == 1)
        if (Ysolution1 < Ysolution2)
            return solution1;
        else return solution2;
    else if(intersectionNum == 2)
        if (Xsolution1 < Xsolution2)
            return solution1;
        else return solution2;
    else {
        if (Xsolution1 > Xsolution2)
            return solution1;
        else return solution2;
    }
}

// This function finds the joint that is in a straight line from jointA to jointB,
// with the length of the bar Lb.
function lineExtension(jointA, jointB, Lb) {
    const xPos_a = jointA[0];
    const yPos_a = jointA[1];
    const xPos_b = jointB[0];
    const yPos_b = jointB[1];
    const theta = Math.atan2((yPos_b - yPos_a) , (xPos_b - xPos_a));
    const X3 = (xPos_b + (Lb * Math.cos(theta)));
    const Y3 = (yPos_b + (Lb * Math.sin(theta)));
    const solution = [X3, Y3];
    return solution;
}

// This function is used in special cases to find the joint at a
// certain angle away from the line jointA to jointB, with bar
// length Lb
function lineExtensionAngle(jointA, jointB, Lb, angle) {
    const xPos_a = jointA[0];
    const yPos_a = jointA[1];
    const xPos_b = jointB[0];
    const yPos_b = jointB[1];
    const theta = Math.atan2((yPos_b - yPos_a) , (xPos_b - xPos_a)) + angle;
    const X3 = (xPos_b + (Lb * Math.cos(theta)));
    const Y3 = (yPos_b + (Lb * Math.sin(theta)));
    const solution = [X3, Y3];
    return solution;
}

const json_data = '{"rearContact":[0,-0],"frontContact":[1211,-0],"rearAxle":[0,340.91176470588238],"frontAxle":[1211,340.91176470588238],"rearGearTop":[2.5441176470588296,369.74509803921563],"frontGearTop":[440.13235294117646,390.0980392156863],"pedalPivot":[440.13235294117646,374.83333333333331],"rearAxlePivot":[109.39705882352945,323.95098039215696],"seatTubePivot":[482.53431372549034,568.186274509804],"rockerSeatStayPivot":[353.63235294117658,566.4901960784315],"fixedShockPt":[480.83823529411774,390.0980392156863],"rockerShockPt":[560.55392156862752,610.58823529411779],"riderCOM":[338.36764705882354,1048.1764705882356],"frame_upperSeatTube":[[482.53431372549034,571.57843137254906],[441.82843137254906,649.59803921568641],[392.64215686274525,734.40196078431381],[353.63235294117658,809.02941176470608],[311.23039215686276,881.96078431372564]],"frame_topBar":[[455.39705882352956,630.9411764705884],[567.33823529411779,678.43137254901967],[658.92647058823547,742.88235294117669],[752.21078431372575,825.99019607843161],[859.06372549019625,917.57843137254929]],"frame_lowerBar":[[440.13235294117646,327.343137254902],[574.122549019608,361.26470588235304],[667.40686274509824,508.82352941176481],[745.42647058823536,651.2941176470589],[864.151960784314,888.74509803921592]],"frame_frontSteerBar":[[840.40686274509824,978.637254901961],[930.29901960784321,807.3333333333336],[1026.9754901960787,617.372549019608],[1111.7794117647061,451.15686274509807],[1179.622549019608,320.55882352941188]],"frame_seat":[[226.42647058823522,863.30392156862763],[212.85784313725495,887.04901960784332],[458.78921568627453,939.62745098039238]],"frame_handlebars":[[800.35784313725514,1000.3039215686276],[900.98529411764719,975.2549019607844]]}';
const linkage_data = JSON.parse(json_data);

function magnitude(p1, p2) {
    const xd = p1[0] - p2[0];
    const yd = p1[1] - p2[1];
    return Math.sqrt(xd * xd + yd * yd);
}

function computeLinkage(angleBar1, rockerStayPivot) {
    const f0 = linkage_data.pedalPivot;
    const bar1 = magnitude(linkage_data.rearAxlePivot, linkage_data.pedalPivot);
    const bar2 = magnitude(linkage_data.rearAxlePivot, rockerStayPivot);
    const bar3 = magnitude(rockerStayPivot, linkage_data.seatTubePivot);
    const f1 = linkage_data.seatTubePivot;
    const f2 = linkage_data.fixedShockPt;
    const bar4 = magnitude(linkage_data.rearAxlePivot, linkage_data.rearAxle);
    const bar5 = magnitude(rockerStayPivot, linkage_data.rearAxle);
    const bar6 = magnitude(rockerStayPivot, linkage_data.rockerShockPt);
    const bar7 = magnitude(linkage_data.seatTubePivot, linkage_data.rockerShockPt);

    const v1 = [Math.cos(angleBar1) * bar1 + f0[0], Math.sin(angleBar1) * bar1 + f0[1]];
    const v2 = circIntersection(v1, f1, bar2, bar3, high);
    const t1 = circIntersection(v1, v2, bar4, bar5, left);
    const t2 = circIntersection(v2, f1, bar6, bar7, high);

    const rgt = circIntersection(t1, v1, magnitude(linkage_data.rearGearTop, linkage_data.rearAxle), magnitude(linkage_data.rearGearTop, linkage_data.rearAxlePivot), high);
    const fgt = linkage_data.frontGearTop;
    const rc = circIntersection(t1, v1, magnitude(linkage_data.rearAxle, linkage_data.rearContact), magnitude(linkage_data.rearAxlePivot, linkage_data.rearContact), low);
    const fc = linkage_data.frontContact;
    const fa = linkage_data.frontAxle;

    return {
        pedalPivot: f0,
        seatTubePivot: f1,
        fixedShockPt: f2,
        rearAxlePivot: v1,
        rockerSeatStayPivot: v2,
        rearAxle: t1,
        rockerShockPt: t2,
        rearContact: rc,
        frontContact: fc,
        frontAxle: fa,
        rearGearTop: rgt,
        frontGearTop: fgt,
        riderCOM: linkage_data.riderCOM
    };
}

// line intercept math by Paul Bourke http://paulbourke.net/geometry/pointlineplane/
// Determine the intersection point of two line segments
// Return FALSE if the lines don't intersect
function intersect(x1, y1, x2, y2, x3, y3, x4, y4) {
    // Check if none of the lines are of length 0
    if ((x1 === x2 && y1 === y2) || (x3 === x4 && y3 === y4)) {
        return false
    }
    denominator = ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1))
    if (denominator === 0) {
        return false     // Lines are parallel
    }
    let ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator
    let ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

    // is the intersection along the segments
    //   if (ua < 0 || ua > 1 || ub < 0 || ub > 1) {
    //       return false
    //   }

    // Return a object with the x and y coordinates of the intersection
      let x = x1 + ua * (x2 - x1)
      let y = y1 + ua * (y2 - y1)

      return [x, y]
}

function intersectLines(p1_l1, p2_l1, p1_l2, p2_l2) {
    return intersect(...p1_l1, ...p2_l1, ...p1_l2, ...p2_l2);
}

function antiSquat(system_CoG, rockerStayPivot, linkage_results) {
    let bike_data = linkage_results;

    // identify instant center / virual pivot point
    let vpp = intersectLines(bike_data.rearAxlePivot, bike_data.pedalPivot, bike_data.seatTubePivot, rockerStayPivot);

    // identify instantaneous force center
    let ifc = intersectLines(bike_data.rearAxle, vpp, bike_data.rearGearTop, bike_data.frontGearTop);

    // identify height at which anti-squat line crosses x-coordinate of front axle
    let ASL_intPt = intersectLines(bike_data.rearContact, ifc, bike_data.frontContact, bike_data.frontAxle);
    let ASLh = ASL_intPt[1];

    // calculate system COM -- assume a 14cm variance
    let antiSquatPercentage = ASLh / system_CoG * 100;
    return {antiSquatPercentage, system_CoG, ASL_intPt};
}

function generate_circle_points(radius, center) {
    return range(361).map(i => {
        const theta = i / 360 * 2 * Math.PI;
        return {x: Math.cos(theta) * radius + center[0], y: Math.sin(theta) * radius + center[1]};
    });
}

// the name of the class must be exactly this
class CustomCenterDiv {
    last_selection_data = null;

    constructor(parent_div) {
        this.angle_range = [Math.PI + 0.15, Math.PI - 0.35];
        this.max_angle_index = 200;
        this.current_angle_index = 0;
        this.animation_paused = false;

        parent_div.innerHTML = `
<div style="margin-top:10px;height:120px">
    <label for="rider_height_input">Rider Height</label>
    <input type="number" name="rider_height_input" id="rider_height_input" value="180" />
    <hr>
    <div>
        <button type="button" id="left-arrow-button">&lt;</button>
        <button type="button" id="pause-button">||</button>
        <button type="button" id="right-arrow-button">&gt;</button>
        <label id="anti_squad_output" style="margin-left: 10px;">Select a point to show something here!</label>
    </div>
    <div style="width: 50%; display: inline-block;">
        <canvas id="linkage-canvas">
        </canvas>
    </div>
    <div style="width: 40%; display: inline-block;">
        <canvas id="antisquat-canvas">
        </canvas>
    </div>
</div>
`;

        const ctx_linkage = document.getElementById("linkage-canvas").getContext('2d');
        this.linkage_chart = Chart.Scatter(ctx_linkage, {
            data: {
                datasets: []
            },
            options: {
                responsive: true,
                maintainAspectRatio: true,
                aspectRatio: 1,
                legend: {
                    display: false
                },
                title: {
                    display: true,
                    text: "Suspension Linkage Analysis"
                },
                plugins: gfopt,
                animation: chartjs_animation
            }
        });
        this.linkage_chart.options.scales.xAxes[0].ticks.min = -400;
        this.linkage_chart.options.scales.xAxes[0].ticks.max = 1600;
        this.linkage_chart.options.scales.yAxes[0].ticks.min = -400;
        this.linkage_chart.options.scales.yAxes[0].ticks.max = 1600;

        const ctx_antisquat = document.getElementById("antisquat-canvas").getContext('2d');
        this.antisquat_chart = Chart.Scatter(ctx_antisquat, {
            data: {
                datasets: []
            },
            options: {
                responsive: true,
                maintainAspectRatio: true,
                aspectRatio: 0.75,
                legend: {
                    display: false
                },
                title: {
                    display: true,
                    text: "Anti-Squat Curve"
                },
                plugins: gfopt,
                animation: chartjs_animation
            }
        });
        this.antisquat_chart.options.scales.xAxes[0].ticks.min = 0;
        this.antisquat_chart.options.scales.xAxes[0].ticks.max = 220;
        this.antisquat_chart.options.scales.xAxes[0].scaleLabel = { display: true, labelString: "Rear Wheel Travel (mm)" };
        this.antisquat_chart.options.scales.yAxes[0].ticks.min =-100;
        this.antisquat_chart.options.scales.yAxes[0].ticks.max = 150;
        this.antisquat_chart.options.scales.yAxes[0].scaleLabel = { display: true, labelString: "Anti-Squat (%)" };

        setInterval(() => this.animation_tick(), 30);

        document.getElementById("left-arrow-button").onclick = () => {this.move_animation(-1); this.update_animation()};
        document.getElementById("right-arrow-button").onclick = () => {this.move_animation(+1); this.update_animation()};
        document.getElementById("pause-button").onclick = () => this.pause_animation();

        document.getElementById('rider_height_input').onchange = () => {this.update_animation(); this.update_antisquat()};
    }

    get_computations(angle) {
        let rider_height = parseFloat(document.getElementById('rider_height_input').value);
        let nRider = (rider_height - 140) / (200-140); // height normalized within permissible range
        let offset = -20 + (20- (-20))*nRider; // lerp(-20,20)
        let system_CoG = linkage_data.riderCOM[1] + offset;

        let app_point = this.last_selection_data[0].selected_points_info[0].application_point;
        let pos_rockerSeatStayPivot = [linkage_data.rockerSeatStayPivot[0], linkage_data.rockerSeatStayPivot[1] + (-9 + app_point*(40 - -9))];

        const linkage_results = computeLinkage(angle, pos_rockerSeatStayPivot);
        const antiSquatData = antiSquat(system_CoG, pos_rockerSeatStayPivot, linkage_results);

        return [linkage_results, antiSquatData];
    }

    update(angle) {
        if (this.last_selection_data == null)
            return;

        const [linkage_results, antiSquatData] = this.get_computations(angle);

        const result_label = document.getElementById('anti_squad_output');
        result_label.textContent = `Current Anti-Squat: ` + Math.round(antiSquatData.antiSquatPercentage * 100) / 100 + "%";

        const rear_wheel_radius = magnitude(linkage_results.rearAxle, linkage_results.rearContact);
        const rear_wheel_points = generate_circle_points(rear_wheel_radius, linkage_results.rearAxle);
        const front_wheel_points = generate_circle_points(rear_wheel_radius, linkage_results.frontAxle);
        const rear_gear_radius = magnitude(linkage_results.rearAxle, linkage_results.rearGearTop)*2;
        const rear_gear_points = generate_circle_points(rear_gear_radius, [linkage_results.rearAxle[0], linkage_results.rearAxle[1]]);
        const pedal_gear_radius = magnitude(linkage_results.frontGearTop, linkage_data.frame_lowerBar[0])
        const pedal_gear_points = generate_circle_points(pedal_gear_radius, linkage_data.frame_lowerBar[0]);
        const upper_chain_points = [rear_gear_points[Math.floor(rear_gear_points.length / 4)], pedal_gear_points[Math.floor(pedal_gear_points.length / 4)]];
        const lower_chain_points = [rear_gear_points[Math.floor(rear_gear_points.length * 3/4)], pedal_gear_points[Math.floor(pedal_gear_points.length *3/4)]];
        const lowerSeatTube = [linkage_data.frame_lowerBar[0], linkage_data.frame_upperSeatTube[0]];

        const generate_vis_dataset = (points, borderWidth = 5) => (
            {
                data: points,
                showLine: true,
                fill: false,
                lineTension: 0,
                borderColor: "#d6d6d6",
                pointRadius: 0,
                borderWidth: borderWidth
            }
        );

        const asl_color = "#55a6f2";
        const asl_dashes = [10,5];

        const HZ_MIN = -1000;
        const HZ_MAX = 10000;
        const VT_MIN = -1000;
        const VT_MAX = 10000;
        const pp = (p) => ({x: p[0], y: p[1]});
        this.linkage_chart.data.datasets = [
            {
                data: [pp(linkage_results.pedalPivot), pp(linkage_results.rearAxlePivot)],
                showLine: true,
                fill: false,
                lineTension: 0,
                pointBackgroundColor: "red",
                borderColor: "green",
                lineWidth: 10
            },
            {
                data: [pp(linkage_results.rearAxlePivot), pp(linkage_results.rearAxle), pp(linkage_results.rockerSeatStayPivot)],
                showLine: true,
                fill: false,
                lineTension: 0.1,
                pointBackgroundColor: "red",
                borderColor: "#edd605",
                lineWidth: 10
            },
            {
                data: [pp(linkage_results.rockerSeatStayPivot), pp(linkage_results.seatTubePivot), pp(linkage_results.rockerShockPt), pp(linkage_results.rockerSeatStayPivot)],
                showLine: true,
                fill: false,
                lineTension: 0.1,
                pointBackgroundColor: "red",
                borderColor: "red",
                lineWidth: 10
            },
            {
                data: [pp(linkage_results.rockerShockPt), pp(linkage_results.fixedShockPt)],
                showLine: true,
                fill: false,
                lineTension: 0,
                pointBackgroundColor: "red",
                borderColor: "gray",
                lineWidth: 10
            },
            {
                data: [pp([HZ_MIN, antiSquatData.system_CoG]), pp([HZ_MAX, antiSquatData.system_CoG])],
                showLine: true,
                fill: false,
                lineTension: 0,
                pointBackgroundColor: asl_color,
                borderColor: asl_color,
                borderDash: asl_dashes
            },
            {
                data: [pp([linkage_results.frontContact[0], VT_MIN]), pp([linkage_results.frontContact[0], VT_MAX])],
                showLine: true,
                fill: false,
                lineTension: 0,
                pointBackgroundColor: asl_color,
                borderColor: asl_color,
                borderDash: asl_dashes
            },
            {
                data: [pp(linkage_results.rearContact), pp(antiSquatData.ASL_intPt)],
                showLine: true,
                fill: false,
                lineTension: 0,
                pointBackgroundColor: asl_color,
                borderColor: asl_color,
                borderDash: asl_dashes
            },
            generate_vis_dataset(rear_wheel_points),
            generate_vis_dataset(front_wheel_points),
            generate_vis_dataset(linkage_data.frame_upperSeatTube.map(pp), 10),
            generate_vis_dataset(linkage_data.frame_topBar.map(pp), 10),
            generate_vis_dataset(linkage_data.frame_lowerBar.map(pp), 10),
            generate_vis_dataset(linkage_data.frame_frontSteerBar.map(pp), 10),
            generate_vis_dataset(linkage_data.frame_seat.map(pp), 10),
            generate_vis_dataset(linkage_data.frame_handlebars.map(pp), 10),
            generate_vis_dataset(lowerSeatTube.map(pp), 10),
            generate_vis_dataset(rear_gear_points, 3),
            generate_vis_dataset(pedal_gear_points, 3),
            generate_vis_dataset(upper_chain_points, 3),
            generate_vis_dataset(lower_chain_points, 3)
        ];
        this.linkage_chart.update();
    }

    animation_tick() {
        if(this.animation_paused)
            return;

        this.move_animation(+1);
        this.update_animation();
    }

    update_animation() {
        if(this.current_angle_index == null || this.last_selection_data == null)
            return;

        const n = this.max_angle_index / 2;
        const angle = this.current_angle_index < n ?
        lerp(this.current_angle_index / n, this.angle_range[0], this.angle_range[1]) :
        lerp((this.current_angle_index - n) / n, this.angle_range[1], this.angle_range[0]);
        this.update(angle);
    }

    pause_animation() {
        this.animation_paused = !this.animation_paused;
    }

    move_animation(delta) {
        this.current_angle_index = this.current_angle_index + delta;
        if(this.current_angle_index > this.max_angle_index)
            this.current_angle_index = 0;
        if(this.current_angle_index < 0)
            this.current_angle_index = 0;
    }

    update_antisquat() {
        const points = range(100).map(i => {
            const theta = lerp(i / 100, this.angle_range[0], this.angle_range[1]);
            const [linkage_results, antiSquatData] = this.get_computations(theta);
            return {x: linkage_results.rearContact[1], y: antiSquatData.antiSquatPercentage};
        });

        this.antisquat_chart.data.datasets = [{
            data: points,
            showLine: true,
            fill: false,
            lineTension: 0,
            pointBackgroundColor: "red",
            borderColor: "black",
            pointRadius: 0,
            borderWidth: 3
        }];
        this.antisquat_chart.options.scales.xAxes[0].ticks.min = 0;
        this.antisquat_chart.options.scales.xAxes[0].ticks.max = 220;
        this.antisquat_chart.options.scales.xAxes[0].scaleLabel = { display: true, labelString: "Rear Wheel Travel (mm)" };
        this.antisquat_chart.options.scales.yAxes[0].ticks.min =-100;
        this.antisquat_chart.options.scales.yAxes[0].ticks.max = 150;
        this.antisquat_chart.options.scales.yAxes[0].scaleLabel = { display: true, labelString: "Anti-Squat (%)" };

        this.antisquat_chart.update();
    }

    //selection_data = [{component_name, selected_points_info}], selected_points_info = [{design_point, performance_point, application_point, mesh}]. mesh = [triangle]. triangle = [v1, v2, v3], vi = [x, y, z]
    set_points(selection_data) {
        this.last_selection_data = selection_data;
        this.update_animation();
        this.update_antisquat();
    }
}