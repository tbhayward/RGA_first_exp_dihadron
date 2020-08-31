/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

//import org.apache.commons.math3.*;

public class analysis_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public analysis_fitter(double beam) {
        super(beam);
        mybeam = beam;
    } 
    
    ////////////////////////////////////////////////////////////////////////////////////////////////  
    // General cuts to test the validity of the event and particle
    public boolean banks_test(DataEvent event) {
        return event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Track") && 
            event.hasBank("REC::Traj") && event.hasBank("RUN::config");
    }
    
    public boolean forward_detector_cut(int current_Part, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", current_Part);
//        System.out.println(status);
        return (Math.abs(status)<4000 || Math.abs(status)>4999) && Math.abs(status)>1999;
    }
    
    public boolean highest_e_in_fd_cut(int p_max_index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", p_max_index);
        return status<-1999;
    }
    
    public boolean theta_cut(int current_Part, HipoDataBank rec_Bank) {
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double r = Math.pow(px*px + py*py + pz*pz, 0.5);
        double theta = (180/Math.PI)*Math.acos(pz/r);
        
        return 5<theta && theta<35;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Fiducial volume cuts on the pcal and drift chambers
    
    public boolean pcal_fiducial_cut(int current_Part, HipoDataBank cal_Bank) {
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part && 
                    cal_Bank.getInt("layer", current_Row)==1) {
                float lv = cal_Bank.getFloat("lv", current_Row);
                float lw = cal_Bank.getFloat("lw", current_Row);
                return lv > 9 && lw > 9; // "loose cut, change to 14 and 19 for medium and tight"
                
            }
        }
        return false;
    }

    public boolean dc_fiducial_cut(int current_Part, HipoDataBank rec_Bank, 
    HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        // trajectory crosses  (layers: 6 = DC region 1 start,  18 = DC region 2 center,  
        // 36 = DC region 3 end)
        
        
        int event_number = run_Bank.getInt("event", 0);
        //fitted values
    double[][][][] maxparams_had = {
        {{{-35.1716, 25.102, -0.750281, 5.34679e-05},{-39.1633, 28.5551, -1.13429, 0.00419047},{-33.7705, 24.8068, -0.811239, 0.00138345}},
      {{-36.2389, 26.7979, -1.08147, 0.0050898},{-43.643, 31.6783, -1.49203, 0.00872922},{-54.4042, 40.6516, -2.52393, 0.0205649}},
      {{-38.3238, 26.1667, -0.777077, 0.000264835},{-34.2011, 24.2843, -0.696392, 3.75866e-12},{-36.4636, 25.8712, -0.786592, 2.24421e-10}},
      {{-31.8019, 23.154, -0.653992, 2.69968e-05},{-34.6637, 24.6043, -0.714901, 2.02675e-10},{-36.7209, 26.2469, -0.828638, 0.000340435}},
      {{-33.4016, 24.6901, -0.779889, 0.000430557},{-35.4583, 24.7491, -0.707953, 2.18559e-10},{-37.7335, 28.1547, -1.1986, 0.00582395}},
      {{-34.7808, 24.6988, -0.719936, 5.73299e-10},{-54.5797, 40.9138, -2.57493, 0.0213354},{-38.4972, 28.3142, -1.21741, 0.00640373}}},
    {{{-2.25358e-08, 12.631, -0.767619, 0.00739811},{-8.09501, 15.9098, -0.844083, 0.00667995},{-1.48113e-06, 12.2061, -0.73167, 0.0074309}},
      {{-2.10872e-07, 12.6689, -0.765156, 0.00720044},{-4.88862, 14.0376, -0.687202, 0.00506307},{-4.59793e-06, 11.5553, -0.591469, 0.00536957}},
      {{-1.13504e-08, 12.6011, -0.746025, 0.00687498},{-6.97884, 15.1788, -0.765889, 0.00570532},{-1.29468, 12.3844, -0.667561, 0.00619226}},
      {{-2.91953e-09, 13.883, -0.999624, 0.0104257},{-4.9855, 13.8864, -0.661348, 0.0048371},{-2.29438e-08, 11.8341, -0.668486, 0.00669247}},
      {{-2.02824e-08, 13.3855, -0.91158, 0.00926769},{-3.29092e-08, 10.8294, -0.382323, 0.00178367},{-4.59027e-06, 11.9414, -0.663872, 0.00625769}},
      {{-3.73322e-09, 12.6126, -0.723548, 0.0062217},{-4.56248, 14.1574, -0.727805, 0.00560108},{-2.39381e-08, 12.0663, -0.6651, 0.00602544}}},
    {{{-1.45923e-08, 13.0297, -0.828302, 0.00795271},{-5.41905, 13.2753, -0.503236, 0.00255607},{-3.67719, 12.1358, -0.462905, 0.00308219}},
      {{-9.953e-10, 11.549, -0.52816, 0.00378771},{-8.47154, 15.9863, -0.826166, 0.0062936},{-6.43715, 13.9081, -0.618535, 0.0046102}},
      {{-4.68458e-08, 12.9481, -0.781613, 0.00689754},{-3.46617, 12.2786, -0.440121, 0.00205448},{-4.43519, 10.9372, -0.210059, 3.69283e-10}},
      {{-4.18414e-07, 13.1542, -0.811251, 0.00714402},{-4.63166, 13.7769, -0.657207, 0.0047586},{-1.99278e-05, 11.3993, -0.575232, 0.00532141}},
      {{-7.07189e-10, 13.2814, -0.88476, 0.00874389},{-5.08373, 14.4384, -0.750795, 0.00586116},{-6.9642e-05, 9.50651, -0.189316, 3.07274e-06}},
      {{-5.85515e-08, 12.5116, -0.688741, 0.00557297},{-1.86306, 11.985, -0.482567, 0.00279836},{-4.94295e-07, 10.1342, -0.316715, 0.00176254}}},
    {{{-0.0157256, 11.1508, -0.415185, 0.00186904},{-13.6561, 19.4418, -1.15773, 0.00989432},{-6.24969e-07, 10.5776, -0.329325, 0.00103488}},
      {{-2.5686e-08, 11.4797, -0.476772, 0.00264288},{-0.0475099, 10.1207, -0.244786, 3.13032e-06},{-4.6875e-07, 12.019, -0.63598, 0.00543214}},
      {{-0.00702545, 11.1294, -0.407207, 0.00171263},{-7.27687, 15.5, -0.807858, 0.0062086},{-5.15078, 12.6368, -0.348584, 9.2687e-12}},
      {{-8.14106e-08, 13.28, -0.818164, 0.00703758},{-7.60722, 14.4871, -0.588662, 0.00326244},{-1.70764e-06, 12.0413, -0.63961, 0.00541784}},
      {{-1.09281, 11.5573, -0.41311, 0.00155228},{-3.71599, 12.8335, -0.521472, 0.00296792},{-0.000410815, 12.4833, -0.72999, 0.0066601}},
      {{-0.652641, 12.2766, -0.554202, 0.00314615},{-8.42824, 15.5087, -0.710609, 0.00447051},{-14.9692, 21.5885, -1.47528, 0.0136615}}},
    {{{-5.58945, 17.4004, -1.34516, 0.0142099},{-14.9585, 20.4538, -1.25118, 0.0106617},{-12.0069, 16.4545, -0.727162, 0.00495418}},
      {{-7.03048, 17.3519, -1.1831, 0.0111308},{-7.30641, 15.8503, -0.850952, 0.00648446},{-10.2549, 15.6139, -0.648352, 0.00380506}},
      {{-9.73111e-09, 13.498, -0.932479, 0.00939708},{-8.38053, 15.5588, -0.711323, 0.00433827},{-12.3097, 16.6403, -0.741362, 0.0050708}},
      {{-7.38905, 17.2652, -1.15517, 0.0109165},{-1.11835e-07, 10.4637, -0.301972, 0.000612754},{-12.2182, 17.4958, -0.919555, 0.00747512}},
      {{-0.492676, 14.4148, -1.0959, 0.0116708},{-5.34309, 14.3258, -0.691954, 0.00480109},{-12.5443, 16.1047, -0.59594, 0.00280171}},
      {{-4.08375e-07, 12.2846, -0.655278, 0.00525956},{-8.93101, 16.4266, -0.861853, 0.00644623},{-11.8406, 17.0417, -0.826301, 0.00596028}}},
    {{{-9.29415, 16.5566, -0.831923, 0.00562504},{-0.954483, 10.5813, -0.265766, 3.24615e-05},{-6.87423, 14.892, -0.76495, 0.00639603}},
      {{-18.8913, 19.3123, -0.711917, 0.00227889},{-13.9788, 18.5678, -0.940183, 0.00664397},{-11.7696, 18.3415, -1.04368, 0.0083506}},
      {{-3.82873, 12.7727, -0.425968, 0.000789835},{-9.81221, 14.6531, -0.471092, 0.00131406},{-14.2392, 15.9895, -0.430525, 2.20712e-12}},
      {{-1.76975e-07, 11.4006, -0.420134, 0.00141302},{-3.11764, 10.9707, -0.245823, 2.23044e-12},{-17.6005, 22.2881, -1.39992, 0.0117791}},
      {{-0.767518, 11.6824, -0.456275, 0.00214005},{-5.28047, 12.65, -0.350658, 9.80081e-05},{-0.0888832, 11.508, -0.49197, 0.00301269}},
      {{-4.72388, 15.8507, -1.00574, 0.00876768},{-2.80649, 11.4056, -0.301812, 0.000190262},{-13.0484, 18.665, -1.08614, 0.00960977}}}};

    double[][][][] minparams_had = {
        {{{37.289, -27.5201,1.12866, -0.00526111},{45.3103, -33.5226,1.72923, -0.0114495},{61.5709, -47.6158,3.4295, -0.0316429}},
      {{36.6259, -27.4064,1.16617, -0.00604629},{50.3751, -37.5848,2.19621, -0.0169241},{35.1563, -26.514,1.09795, -0.00545864}},
      {{27.2367, -20.3068,0.517752, -0.000335432},{39.0489, -28.6903,1.24306, -0.0065226},{41.0208, -30.0339,1.30776, -0.00626721}},
      {{29.261, -21.7041,0.613556, -0.000774652},{39.5304, -29.1388,1.34116, -0.00823818},{44.5313, -33.4056,1.77581, -0.0123965}},
      {{36.5659, -25.119,0.714074, -2.65397e-11},{31.6524, -22.6934,0.613977, -5.46634e-10},{34.7312, -24.9901,0.749061, -1.22922e-09}},
      {{33.154, -23.8803,0.685794, -1.13236e-10},{42.6731, -31.0799,1.40425, -0.00730816},{46.4732, -35.6988,2.10144, -0.0164771}}},
    {{{2.40451, -15.0848,1.05504, -0.0103356},{8.93825, -16.5995,0.925874, -0.00767902},{7.23439e-08, -12.5963,0.814574, -0.00864749}},
      {{6.2953e-07, -12.6365,0.732206, -0.00639165},{12.6866, -18.7831,1.0952, -0.00923029},{3.12805e-07, -12.5395,0.795535, -0.00828991}},
      {{2.69495, -14.8778,1.00751, -0.00975373},{6.05446, -14.6778,0.767457, -0.00636729},{3.94741e-07, -11.1038,0.524109, -0.00471514}},
      {{2.31558e-07, -11.5073,0.494316, -0.00303611},{5.66995, -14.5948,0.740956, -0.00561851},{4.40475e-06, -9.57062,0.20354, -0.000213213}},
      {{2.74277e-08, -13.3573,0.886651, -0.00857992},{9.98849e-05, -11.524,0.531486, -0.00391441},{8.50811e-07, -9.72224,0.240264, -0.000781498}},
      {{6.9021e-08, -11.8859,0.53864, -0.00325092},{10.0169, -16.9153,0.921593, -0.00752414},{9.90518e-07, -11.9578,0.697029, -0.00717645}}},
    {{{6.87966e-10, -12.8497,0.757379, -0.00651612},{16.2087, -19.3776,0.951508, -0.00645029},{14.513, -18.8625,1.05542, -0.00918985}},
      {{1.07197e-07, -12.5469,0.703086, -0.00585238},{0.0871522, -9.22628,0.159628, -0.000343326},{12.1181, -17.5575,0.940249, -0.00788125}},
      {{2.10191e-09, -12.2782,0.661926, -0.00555279},{12.5105, -17.9998,0.951807, -0.00732845},{12.8043, -17.8322,0.972401, -0.00841528}},
      {{8.11926e-10, -12.7225,0.737941, -0.00647355},{7.50649, -15.987,0.889398, -0.00729282},{0.174541, -10.0266,0.306882, -0.00186093}},
      {{3.81202e-09, -12.0926,0.598943, -0.00430458},{8.72368, -17.2511,1.06348, -0.00953327},{1.5205, -9.86713,0.183806, -6.40377e-12}},
      {{1.37378e-07, -12.9247,0.769722, -0.00664936},{8.53877, -16.6167,0.946138, -0.00788745},{8.47417, -14.3897,0.581492, -0.00387111}}},
    {{{2.50079e-07, -12.5209,0.678491, -0.00528954},{12.6171, -18.4205,1.01802, -0.00807702},{10.4903, -18.0981,1.10546, -0.00971519}},
      {{5.87069e-07, -12.0075,0.585538, -0.00416654},{11.1348, -17.5468,0.943652, -0.00729083},{0.949201, -10.5869,0.267536, -6.04802e-05}},
      {{1.14857, -11.1478,0.345528, -0.000841836},{10.9482, -17.1647,0.909605, -0.00722404},{8.7569e-08, -10.4446,0.316302, -0.00101964}},
      {{1.09759e-06, -11.5019,0.48435, -0.00277852},{0.637937, -10.7065,0.316211, -0.000801127},{5.67144e-07, -12.88,0.831252, -0.00835441}},
      {{1.68853, -11.2582,0.308152, -7.81686e-12},{9.44238, -17.1892,1.00561, -0.00864837},{1.20713e-07, -12.2246,0.669321, -0.0057622}},
      {{0.00217558, -10.8858,0.347928, -0.000790679},{11.8583, -17.6423,0.923581, -0.00703041},{3.24078, -13.4024,0.668777, -0.00504175}}},
    {{{6.04158, -16.8155,1.13335, -0.0105359},{8.24786, -17.0204,1.05097, -0.00941875},{11.7617, -17.202,0.864472, -0.00649032}},
      {{3.70947, -13.0663,0.513818, -0.00222627},{16.7022, -21.9618,1.42869, -0.012705},{6.8993, -14.8192,0.740813, -0.00585407}},
      {{2.18472e-06, -11.9461,0.583354, -0.00423414},{6.51489e-07, -10.5669,0.353028, -0.00166977},{12.5113, -16.5038,0.709888, -0.00471964}},
      {{0.812719, -11.3245,0.390183, -0.00134086},{2.97251, -11.9374,0.338592, -4.36096e-13},{13.8844, -17.5707,0.818446, -0.00581811}},
      {{1.55496, -14.4569,0.949497, -0.00857237},{0.34359, -10.5041,0.286497, -0.000346977},{14.4141, -18.7457,1.01652, -0.00845189}},
      {{1.26317e-08, -11.1424,0.434251, -0.00236267},{6.58119, -15.8546,0.930324, -0.00801288},{4.41865, -11.1991,0.234652, -7.43723e-10}}},
    {{{6.87926, -12.8949,0.334733, -6.38494e-06},{35.2336, -32.2007,2.21489, -0.020555},{6.80949, -16.8945,1.19056, -0.0127558}},
      {{0.95782, -12.4625,0.599979, -0.00405342},{20.4051, -23.1936,1.42408, -0.0120792},{10.277, -16.1457,0.785186, -0.00612069}},
      {{0.236196, -11.6165,0.458613, -0.002018},{12.8771, -19.6785,1.26163, -0.0115917},{5.21194e-08, -12.551,0.78718, -0.00794713}},
      {{8.40778, -14.9001,0.534967, -0.00147246},{15.9376, -20.9945,1.2908, -0.0110556},{10.4773, -16.2238,0.783386, -0.00593478}},
      {{3.21187, -12.1221,0.348938, -8.70415e-14},{13.8983, -19.1128,1.04727, -0.00797426},{11.6342, -18.8428,1.18853, -0.0107619}},
      {{3.7311, -12.4292,0.419345, -0.00134704},{6.92884, -13.2494,0.391862, -0.000767396},{5.5939, -14.4175,0.729195, -0.00568477}}}};

   
   double[][][][] maxparams_elec =
        {{{{-7.49907, 0.583375},{-18.8174, 0.599219},{-23.9353, 0.574699}},
     {{-14.0547, 0.631533},{-14.4223, 0.597079},{-14.838, 0.547436}},
     {{-7.72508, 0.578501},{-18.7928, 0.56725},{-29.9003, 0.612354}},
     {{-6.12844, 0.566777},{-13.6772, 0.573262},{-26.1895, 0.591816}},
     {{-20.0718, 0.670941},{-9.4775, 0.511748},{-28.0869, 0.590488}},
     {{-9.52924, 0.591687},{-17.8564, 0.596417},{-23.5661, 0.576317}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}}};

    double[][][][] minparams_elec =
        {{{{7.62814, -0.56319},{18.2833, -0.587275},{20.2027, -0.54605}},
     {{9.20907, -0.586977},{10.493, -0.544243},{23.0759, -0.581959}},
     {{12.5459, -0.631322},{20.5635, -0.618555},{26.3621, -0.576806}},
     {{8.36343, -0.552394},{14.7596, -0.554798},{29.5554, -0.60545}},
     {{16.3732, -0.663303},{10.0255, -0.533019},{31.6086, -0.617053}},
     {{8.20222, -0.567211},{20.0181, -0.605458},{22.2098, -0.567599}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
  {{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
    {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}}};


        boolean track_success = true; 
        for (int current_Row = 0; current_Row < traj_Bank.rows(); current_Row++) {
            if (!track_success) { continue; }
            // loop over all entries in the trajectory bank
            if (current_Part == traj_Bank.getInt("pindex", current_Row)) {
                // require that the particle examined corresponds to this track
                int region = -1;
                if (traj_Bank.getInt("layer", current_Row) == 6) {
                    region = 0;
                } else if (traj_Bank.getInt("layer", current_Row) == 18) {
                    region = 1;
                } else if (traj_Bank.getInt("layer", current_Row) == 36) {
                    region = 2;
                } // note these are CLAS region - 1 because of java indexing rules

                double theta_DCr = 5000;
                double phi_DCr_raw = 5000;
                if (region == 0 || region == 1 || region == 2) {
                    float cx = traj_Bank.getFloat("cx", current_Row);
                    float cy = traj_Bank.getFloat("cy", current_Row);
                    float x = traj_Bank.getFloat("x", current_Row);
                    float y = traj_Bank.getFloat("y", current_Row);
                    float z = traj_Bank.getFloat("z", current_Row);
                    float cz = traj_Bank.getFloat("cz", current_Row);
                    double r = Math.sqrt(x*x + y*y + z*z);
                    theta_DCr = (180/Math.PI) * Math.acos(z/r);
                    phi_DCr_raw = (180/Math.PI) * Math.atan2(y/r,x/r);
                    
                    int sector = 0;
                    if(phi_DCr_raw < 30 && phi_DCr_raw >= -30){        sector = 0;} // note this is CLAS sector - 1
                    else if(phi_DCr_raw < 90 && phi_DCr_raw >= 30){    sector = 1;} // because of Java indexing reasons
                    else if(phi_DCr_raw < 150 && phi_DCr_raw >= 90){   sector = 2;}
                    else if(phi_DCr_raw >= 150 || phi_DCr_raw < -150){ sector = 3;}
                    else if(phi_DCr_raw < -90 && phi_DCr_raw >= -150){ sector = 4;}
                    else if(phi_DCr_raw < -30 && phi_DCr_raw >= -90){  sector = 5;}
                    double phi_DCr = 5000;
                    if (sector == 0) {phi_DCr = phi_DCr_raw;}
                    if (sector == 1) {phi_DCr = phi_DCr_raw - 60;}
                    if (sector == 2) {phi_DCr = phi_DCr_raw - 120;}
                    if (sector == 3 && phi_DCr_raw > 0) {phi_DCr = phi_DCr_raw - 180;}
                    if (sector == 3 && phi_DCr_raw < 0) {phi_DCr = phi_DCr_raw + 180;}
                    if (sector == 4) {phi_DCr = phi_DCr_raw + 120;}
                    if (sector == 5) {phi_DCr = phi_DCr_raw + 60;}
                    
                    
//                    System.out.println(phi_DCr_raw+" "+phi_DCr);

                    int pid = rec_Bank.getInt("pid", current_Part);
                    switch (pid) { // these are to set the index in the paramter lists
                    case 11:
                        pid = 0;
                        break;

                    case 2212:
                        pid = 1;
                        break;

                    case 211:
                        pid = 2;
                        break;

                    case -211:
                        pid = 3;
                        break;

                    case 321:
                        pid = 4;
                        break;

                    case -321:
                        pid = 5;
                        break;
                    }
                    
//                    if (pid==0 && sector!=0) {
//                        return false;
//                    }

                    int pdg_pid = rec_Bank.getInt("pid", current_Part);
                    if (rec_Bank.getInt("pid", current_Part) != 11) { // hadrons
                        double calc_phi_min = minparams_had[pid][sector][region][0] + 
                            minparams_had[pid][sector][region][1]*Math.log(theta_DCr) + 
                            minparams_had[pid][sector][region][2]*theta_DCr +
                            minparams_had[pid][sector][region][3]*theta_DCr*theta_DCr;

                        double calc_phi_max = maxparams_had[pid][sector][region][0] + 
                            maxparams_had[pid][sector][region][1]*Math.log(theta_DCr) + 
                            maxparams_had[pid][sector][region][2]*theta_DCr +
                            maxparams_had[pid][sector][region][3]*theta_DCr*theta_DCr; 

                        
                        track_success = phi_DCr > calc_phi_min && phi_DCr < calc_phi_max;
//                        if (track_success == false) {
//                        System.out.println("event: "+event_number+", pid: "+pdg_pid+", sector: "+(sector+1)+", region: "+(region+1)+", theta_DCr: "+theta_DCr+", calc_phi_min: "+calc_phi_min+", phi_DCr: "+phi_DCr+", calc_phi_max: "+calc_phi_max+", track success: "+track_success);
//                        }
                    } else { // electrons
                        double x_New = 10000;
                        double y_New = 10000;
                        switch (sector) {
                            case 0:
                                x_New = x;
                                y_New = y;
                                break;
                            case 1:
                                x_New = x*Math.cos(-60*Math.PI/180) - y*Math.sin(-60*Math.PI/180);
                                y_New = x*Math.sin(-60*Math.PI/180) + y*Math.cos(-60*Math.PI/180);      
                                break;
                            case 2:
                                x_New = x*Math.cos(-120*Math.PI/180) - y*Math.sin(-120*Math.PI/180);
                                y_New = x*Math.sin(-120*Math.PI/180) + y*Math.cos(-120*Math.PI/180);      
                                break;
                            case 3:
                                x_New = x*Math.cos(-180*Math.PI/180) - y*Math.sin(-180*Math.PI/180);
                                y_New = x*Math.sin(-180*Math.PI/180) + y*Math.cos(-180*Math.PI/180);      
                                break;
                            case 4:
                                x_New = x*Math.cos(120*Math.PI/180) - y*Math.sin(120*Math.PI/180);
                                y_New = x*Math.sin(120*Math.PI/180) + y*Math.cos(120*Math.PI/180);      
                                break;
                            case 5:
                                x_New = x*Math.cos(60*Math.PI/180) - y*Math.sin(60*Math.PI/180);
                                y_New = x*Math.sin(60*Math.PI/180) + y*Math.cos(60*Math.PI/180);      
                                break;
                        }
                        double calc_min = minparams_elec[pid][sector][region][0]+minparams_elec[pid][sector][region][1]*x_New;
                        double calc_max = maxparams_elec[pid][sector][region][0]+maxparams_elec[pid][sector][region][1]*x_New;
                        
                        track_success = y_New > calc_min && y_New < calc_max;
//                        System.out.println("event: "+event_number+", pid: "+pdg_pid+", sector: "+(sector+1)+", region: "+(region+1)+", x: "+x+", y: "+y+", calc_min: "+calc_min+", y_New: "+y_New+", calc_max: "+calc_max+", track success: "+track_success);
                    }
                    
                
                }
            }

        }
        return track_success;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // PID enhancements
    
    public boolean calorimeter_energy_cut(int current_Part, HipoDataBank cal_Bank) {
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part && 
                    cal_Bank.getInt("layer", current_Row)==1) {
                return cal_Bank.getFloat("energy", current_Row) > 0.07;
            }
        }
        return false;
    }
    
//  TWO DIFFERENT SAMPLING FRACTION FUNCTIONS: ONE > 0.17 AND ONE BASED ON STD AWAY FROM MEAN
//    public boolean calorimeter_sampling_fraction_cut(int current_Part, double p, HipoDataBank cal_Bank) {
//        double cal_energy = 0;
//        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
//            if (cal_Bank.getInt("pindex", current_Row)==current_Part)  {
//                cal_energy+= cal_Bank.getFloat("energy", current_Row);
//            }
//        }
//        return cal_energy/p > 0.17;
//    }
    
    public boolean calorimeter_sampling_fraction_cut(int current_Part, double p, HipoDataBank cal_Bank) {
        double scale = 3.5; // how many std away from mean to cut on
        int sector = -1;
        double cal_energy = 0;
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part)  {
                sector = cal_Bank.getInt("sector", current_Row) - 1; // subtract one to start at index of 0
                cal_energy+= cal_Bank.getFloat("energy", current_Row);
            }
        }
        
        double e_cal_sampl_mu[][] = {{0.2531, 0.2550, 0.2514, 0.2494, 0.2528, 0.2521},
            {-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703},
            {4.939, 5.350, 5.102, 6.440, 6.149, 4.957}};
        
        double e_cal_sampl_sigma[][] = {{0.002726, 0.004157, 0.00522, 0.005398, 0.008453, 0.006553}, 
            {1.062, 0.859, 0.5564, 0.6576, 0.3242, 0.4423}, 
            {-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274}};
        
        double mean = e_cal_sampl_mu[0][sector]+(e_cal_sampl_mu[1][sector]/1000)*(p-e_cal_sampl_mu[2][sector])*
                (p-e_cal_sampl_mu[2][sector]);
        
        double std = e_cal_sampl_sigma[0][sector] + e_cal_sampl_sigma[1][sector] / 
                (10 * (p-e_cal_sampl_sigma[2][sector]));
        
        // sampling fraction is cal_energy/p
        return ((cal_energy/p) > (mean-scale*std)) && ((cal_energy/p) < (mean+scale*std));
    }
    
    public boolean calorimeter_diagonal_cut(int current_Part, double p, HipoDataBank cal_Bank) {
        if (p < 4.5) {
            return true; // only apply diagonal cut above 4.5 GeV
        }
        double pcal_plus_ecal_inner = 0;
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part 
                && (cal_Bank.getInt("layer", current_Row)==1 || cal_Bank.getInt("layer", current_Row)==4 ) )  {
                pcal_plus_ecal_inner+=cal_Bank.getFloat("energy", current_Row);
            }
        }
        return 0.2 < pcal_plus_ecal_inner/p;
    }
    
    public boolean electron_z_vertex_cut(float vz) {
        return vz>-13 && vz<12;
    }
    
    public boolean pion_z_vertex_cut(float vz, double trigger_electron_vz) {
        return Math.abs(trigger_electron_vz - vz) < 20;
    }
    
    public boolean pion_chi2pid_cut(int current_Part, HipoDataBank rec_Bank) {
        float chi2pid = rec_Bank.getFloat("chi2pid", current_Part);
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double p = Math.sqrt(px*px+py*py+pz*pz);
        
        int pid = rec_Bank.getInt("pid", current_Part); // slightly different cuts for pi+ and pi-
        double C;
        if (pid==211) { 
            C = 0.88;  
        } else {
            C = 0.93;
        }
        
        if (p<2.44) {
            return chi2pid < 3*C && chi2pid > -3*C;
        } else {
            return chi2pid < C * (0.00869 + 14.98587*Math.exp(-p/1.18236) + 1.81751*Math.exp(-p/4.86394))
                && chi2pid > -3*C;
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    public boolean particle_test(int current_Part, HipoDataBank rec_Bank) {
        return true
            && forward_detector_cut(current_Part, rec_Bank) 
            && theta_cut(current_Part, rec_Bank)
            ;
    }
    
    public boolean electron_test(int current_Part, double p, float vz, HipoDataBank rec_Bank, HipoDataBank cal_Bank, 
            HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        return true
            && p > 2.0
            && calorimeter_energy_cut(current_Part, cal_Bank) 
            && calorimeter_sampling_fraction_cut(current_Part, p, cal_Bank)
            && calorimeter_diagonal_cut(current_Part, p, cal_Bank)
            && electron_z_vertex_cut(vz)
            && pcal_fiducial_cut(current_Part, cal_Bank)
            && dc_fiducial_cut(current_Part, rec_Bank, track_Bank, traj_Bank, run_Bank)
                ;
    }
    
    public boolean pion_test(int current_Part, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
            && p > 1.25
            && pion_z_vertex_cut(vz, trigger_electron_vz)
            && pion_chi2pid_cut(current_Part, rec_Bank)
            && dc_fiducial_cut(current_Part, rec_Bank, track_Bank, traj_Bank, run_Bank)
              ;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) {
        
        if (banks_test(event)) {
            PhysicsEvent physEvent = new PhysicsEvent();
            HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle"); // load particle bank
            HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
            HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
            HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
            HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
            
            double trigger_electron_vz = -99;
            double p_max = 0;
            int p_max_index = -99; // find the index of the highest energy electorn before FD cut
            
            
            // cycle over the particles in recBank and investigate electron and pion IDs
            for (int current_Part = 0; current_Part < rec_Bank.rows(); current_Part++) {
                int pid = rec_Bank.getInt("pid", current_Part);
                if (pid!=11) { continue; }
//                System.out.println(rec_Bank.getInt("status", current_Part));
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", current_Part);
                float py = rec_Bank.getFloat("py", current_Part);
                float pz = rec_Bank.getFloat("pz", current_Part);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",current_Part);
                float vy = rec_Bank.getFloat("vy",current_Part);
                float vz = rec_Bank.getFloat("vz",current_Part);
                if (p>p_max) {
                    p_max = p;
                    p_max_index = current_Part;
                    trigger_electron_vz = vz;
                }
            }
            if (p_max_index >= 0 && highest_e_in_fd_cut(p_max_index, rec_Bank)) {
//                System.out.println(p_max);
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", p_max_index);
                float py = rec_Bank.getFloat("py", p_max_index);
                float pz = rec_Bank.getFloat("pz", p_max_index);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",p_max_index);
                float vy = rec_Bank.getFloat("vy",p_max_index);
                float vz = rec_Bank.getFloat("vz",p_max_index);
                if (particle_test(p_max_index, rec_Bank)
                        && electron_test(p_max_index, p, vz, rec_Bank, cal_Bank, track_Bank, traj_Bank, 
                                run_Bank)) {

                    Particle part = new Particle(11,px,py,pz,vx,vy,vz);
                    physEvent.addParticle(part);

                }
            }
            
            
            for (int current_Part = 0; current_Part < rec_Bank.rows(); current_Part++) {
                int pid = rec_Bank.getInt("pid", current_Part);
                
                if ((Math.abs(pid)!=211) || trigger_electron_vz == -99) { continue; }
                // requires the particle to be pion by EventBuilder and for an electron to have been assigned to event
                
//                System.out.println(rec_Bank.getInt("status", current_Part));
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", current_Part);
                float py = rec_Bank.getFloat("py", current_Part);
                float pz = rec_Bank.getFloat("pz", current_Part);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",current_Part);
                float vy = rec_Bank.getFloat("vy",current_Part);
                float vz = rec_Bank.getFloat("vz",current_Part);
               if (particle_test(current_Part, rec_Bank) 
                    && pion_test(current_Part, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
//                if (pion_test(current_Part, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
//                    track_Bank, traj_Bank, run_Bank)) {
                   
                   Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                   physEvent.addParticle(part);   
               }
            }
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}