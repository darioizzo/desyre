// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_MISC_DATA_H
#define DSYRE_MISC_DATA_H

#include <vector>

namespace dsyre
{
namespace gym
{
// Points from some temperature profile used in industrial applications (from Luca Gui (dow corporation))
inline std::vector<std::vector<double>> luca1_points = {
    {0},     {60},    {120},   {180},   {240},   {300},   {360},   {420},   {480},   {540},   {600},   {660},   {720},
    {780},   {840},   {900},   {960},   {1020},  {1080},  {1140},  {1200},  {1260},  {1320},  {1380},  {1440},  {1500},
    {1560},  {1620},  {1680},  {1740},  {1800},  {1860},  {1920},  {1980},  {2040},  {2100},  {2160},  {2220},  {2340},
    {2400},  {2460},  {2520},  {2580},  {2640},  {2700},  {2760},  {2820},  {2880},  {2940},  {3000},  {3060},  {3120},
    {3180},  {3240},  {3300},  {3360},  {3420},  {3480},  {3540},  {3600},  {3660},  {3720},  {3780},  {3840},  {3900},
    {3960},  {4020},  {4080},  {4140},  {4200},  {4260},  {4320},  {4380},  {4440},  {4500},  {4560},  {4620},  {4680},
    {4740},  {4800},  {4860},  {4980},  {5040},  {5100},  {5160},  {5220},  {5280},  {5340},  {5400},  {5460},  {5520},
    {5580},  {5640},  {5700},  {5760},  {5820},  {5880},  {5940},  {6000},  {6060},  {6180},  {6240},  {6300},  {6360},
    {6420},  {6480},  {6540},  {6600},  {6660},  {6720},  {6780},  {6840},  {6900},  {6960},  {7020},  {7080},  {7140},
    {7200},  {7260},  {7320},  {7380},  {7440},  {7500},  {7560},  {7620},  {7680},  {7740},  {7800},  {7860},  {7920},
    {7980},  {8040},  {8100},  {8160},  {8220},  {8280},  {8340},  {8400},  {8520},  {8580},  {8640},  {8700},  {8760},
    {8820},  {8880},  {8940},  {9000},  {9060},  {9120},  {9180},  {9240},  {9300},  {9360},  {9420},  {9480},  {9540},
    {9600},  {9660},  {9720},  {9780},  {9840},  {9900},  {9960},  {10020}, {10080}, {10140}, {10200}, {10260}, {10320},
    {10380}, {10440}, {10500}, {10560}, {10620}, {10680}, {10740}, {10800}, {10860}, {10920}, {10980}, {11040}, {11100},
    {11160}, {11220}, {11280}, {11340}, {11400}, {11460}, {11520}, {11580}, {11640}, {11700}, {11760}, {11820}, {11880},
    {11940}, {12000}, {12060}, {12120}, {12180}, {12240}, {12300}, {12360}, {12420}, {12480}, {12540}, {12600}, {12660},
    {12720}, {12780}, {12840}, {12900}, {12960}, {13020}, {13080}, {13140}, {13200}, {13260}, {13320}, {13380}, {13440},
    {13500}, {13560}, {13620}, {13680}, {13740}, {13800}, {13860}, {13920}, {13980}, {14040}, {14100}, {14160}, {14220},
    {14340}, {14400}, {14460}, {14520}, {14580}, {14700}, {14760}, {14820}, {14880}, {14940}, {15000}, {15060}, {15120},
    {15180}, {15240}, {15300}, {15360}, {15420}, {15480}, {15540}, {15600}, {15660}, {15720}, {15780}, {15840}, {15900},
    {15960}, {16020}, {16080}, {16140}, {16200}, {16260}, {16320}, {16380}, {16440}, {16500}, {16560}, {16620}, {16680},
    {16740}, {16800}, {16860}, {16920}, {16980}, {17040}, {17100}, {17160}, {17220}, {17280}, {17340}, {17400}, {17460},
    {17520}, {17580}, {17640}, {17700}, {17760}, {17820}, {17880}, {17940}, {18000}, {18060}, {18120}, {18180}, {18240},
    {18300}, {18360}, {18420}, {18480}, {18540}, {18600}, {18660}, {18720}, {18780}, {18840}, {18900}, {18960}, {19020},
    {19080}, {19140}, {19200}, {19260}, {19320}, {19380}, {19440}, {19500}, {19560}, {19620}, {19680}, {19740}, {19800},
    {19860}, {19920}, {19980}, {20040}, {20100}, {20160}, {20220}, {20280}, {20340}, {20400}, {20460}, {20520}, {20580},
    {20640}, {20700}, {20760}, {20820}, {20880}, {20940}, {21000}, {21060}, {21120}, {21180}, {21240}, {21300}, {21360},
    {21420}, {21480}, {21540}, {21600}, {21660}, {21720}, {21780}, {21840}, {21900}, {21960}, {22020}, {22080}, {22140},
    {22200}, {22260}, {22320}, {22380}, {22440}, {22500}, {22560}, {22620}, {22680}, {22740}, {22800}, {22860}, {22920},
    {23040}, {23100}, {23160}, {23220}, {23280}, {23340}, {23400}, {23460}, {23520}, {23580}, {23640}, {23700}, {23760},
    {23820}, {23880}, {23940}, {24000}, {24060}, {24120}, {24180}, {24240}, {24300}, {24360}, {24420}, {24480}, {24540},
    {24600}, {24660}, {24720}, {24780}, {24840}, {24900}, {24960}, {25020}, {25080}, {25140}, {25200}, {25260}, {25320},
    {25380}, {25440}, {25500}, {25560}, {25620}, {25680}, {25740}, {25800}, {25860}, {25920}, {25980}, {26040}, {26100},
    {26160}, {26220}, {26280}, {26340}, {26400}, {26460}, {26520}, {26580}, {26640}, {26700}, {26760}, {26820}, {26880},
    {26940}, {27000}, {27060}, {27120}, {27180}, {27240}, {27300}, {27360}, {27420}, {27480}, {27540}, {27600}, {27660},
    {27720}, {27780}, {27840}, {27900}, {27960}, {28020}, {28080}, {28140}, {28200}, {28260}, {28320}, {28380}, {28440},
    {28500}, {28560}, {28620}, {28680}, {28740}, {28800}, {28860}, {28920}, {28980}, {29040}, {29100}, {29160}, {29220},
    {29280}, {29340}, {29400}, {29460}, {29520}, {29580}, {29640}, {29700}, {29760}, {29820}, {29880}, {29940}, {30000},
    {30060}, {30120}, {30180}, {30240}, {30300}, {30360}, {30420}, {30480}, {30540}, {30600}, {30660}, {30720}, {30780},
    {30840}, {30900}, {30960}, {31020}, {31080}, {31140}, {31200}, {31260}, {31320}, {31380}, {31440}, {31500}, {31560},
    {31620}, {31680}, {31740}, {31800}, {31860}, {31920}, {31980}, {32040}, {32100}, {32160}, {32220}, {32280}, {32340},
    {32400}, {32460}, {32520}, {32580}, {32640}, {32700}, {32760}, {32820}, {32880}, {32940}, {33000}, {33060}, {33120},
    {33180}, {33240}, {33300}, {33360}, {33420}, {33480}, {33540}, {33600}, {33660}, {33720}, {33780}, {33840}, {33900},
    {33960}, {34020}, {34080}, {34140}, {34200}, {34260}, {34320}, {34380}, {34440}, {34500}, {34560}, {34620}, {34680},
    {34740}, {34800}, {34860}, {34920}, {34980}, {35040}, {35100}, {35160}, {35220}, {35280}, {35340}, {35400}, {35460},
    {35520}, {35580}, {35640}, {35700}, {35760}, {35820}, {35880}, {35940}, {36000}, {36060}, {36120}, {36180}, {36240},
    {36300}, {36360}, {36420}, {36480}, {36540}, {36600}, {36660}, {36720}, {36780}, {36840}, {36900}, {36960}, {37020},
    {37080}, {37140}, {37200}, {37260}, {37320}, {37380}, {37440}, {37500}, {37560}, {37620}, {37680}, {37740}, {37800},
    {37860}, {37920}, {37980}, {38040}, {38100}, {38160}, {38220}};

inline std::vector<double> luca1_labels = {
    337.1442908, 336.9744037, 336.9824174, 336.5960922, 335.255398,  335.5252764, 334.3072227, 334.9410065, 334.4445171,
    332.8497379, 333.9870947, 332.2910747, 334.1350198, 332.0067778, 333.8810337, 331.8969311, 332.7588098, 331.5910934,
    330.9132806, 330.5956971, 332.4400475, 332.3796281, 331.814891,  329.9387639, 329.5636731, 330.0219347, 329.7540461,
    330.6068489, 329.597415,  329.3769882, 330.0688911, 328.3634273, 328.5975643, 328.434937,  328.2632795, 328.1946805,
    326.8960196, 328.6338314, 327.1479591, 326.5556453, 327.9179143, 327.9942847, 328.4017316, 327.2529231, 327.3911318,
    326.6970198, 326.6950885, 325.4782815, 324.8187344, 326.9839012, 325.1536111, 326.8089144, 325.8877202, 324.7849329,
    325.6072575, 324.9948122, 325.4887436, 326.0431928, 323.6113226, 325.9838914, 324.1924508, 324.5970406, 324.8230764,
    324.9237099, 324.2666625, 323.6965632, 321.9786247, 322.8748345, 322.8899082, 323.9822995, 321.1392998, 323.6433758,
    323.2831224, 323.0677065, 322.035195,  320.2140103, 322.0256243, 321.7293684, 320.3009872, 319.6279384, 319.4813218,
    319.9187445, 319.6932013, 319.1582722, 320.7196917, 319.4144958, 318.8384286, 319.9096955, 318.1323327, 319.0491528,
    318.2715763, 319.2434569, 319.3683525, 318.7220745, 317.3729954, 316.1924805, 317.6607081, 317.4401266, 316.8007575,
    316.3831204, 315.3087716, 316.4643641, 315.545849,  315.0653973, 315.4820921, 316.0636043, 315.8050276, 317.15843,
    315.0013158, 316.3060765, 314.7660777, 316.5689015, 314.9915883, 314.8268058, 313.6818342, 313.9843479, 314.85659,
    313.6457786, 313.9740967, 313.055651,  314.9535052, 314.8975508, 314.655686,  313.8751524, 312.7129912, 314.6638683,
    313.5441028, 312.707546,  312.5758631, 313.6651482, 312.0994893, 311.4212523, 313.4662091, 314.06175,   313.9324801,
    311.2821648, 314.2439145, 313.9891592, 312.6230376, 312.6748777, 311.7325204, 311.9418737, 312.6439556, 313.2380496,
    313.1371983, 313.1932226, 312.6721982, 311.6082809, 312.6828219, 311.4888788, 310.3933265, 312.8019786, 311.7547666,
    310.553449,  312.2085369, 312.7470519, 310.6600828, 311.9492027, 311.3306919, 310.1125428, 310.6858603, 310.5374625,
    309.616072,  311.7214697, 310.7412576, 309.255049,  310.9959456, 309.5496629, 309.8922506, 309.2589913, 309.6225286,
    311.4041382, 308.9472363, 309.4356973, 308.6938563, 308.8859748, 311.5421611, 308.8003615, 310.5179815, 310.9961059,
    309.5454753, 308.8832005, 309.0879641, 309.1384432, 308.5626005, 310.3581126, 309.722625,  309.5528934, 309.8827313,
    308.7161182, 309.6421889, 310.9313043, 310.30758,   308.8335661, 309.2443618, 307.8976394, 308.1901883, 309.0817311,
    310.5508054, 308.9074166, 307.7294688, 310.0647507, 308.2555717, 309.0918719, 308.4474127, 308.9598839, 310.5361995,
    309.8934223, 309.2699092, 308.8803725, 308.8688793, 307.6696575, 309.4188047, 309.0433909, 309.3836838, 307.4465841,
    309.3869415, 308.4820491, 309.6412089, 307.6391334, 310.2361889, 307.692014,  307.6061464, 309.9988997, 309.6635524,
    307.2992441, 309.5478944, 309.7384236, 307.6155164, 309.7216317, 307.0322128, 309.5891705, 308.6442173, 308.2093349,
    307.0311527, 307.4621258, 308.9132837, 307.6651043, 308.9875313, 306.7033955, 307.7732745, 307.8833902, 307.8650017,
    309.1836698, 307.1342834, 308.7995078, 308.7250358, 307.0102584, 307.4185696, 306.7896699, 306.3319703, 307.0610922,
    307.4909859, 306.189724,  308.5292718, 306.9459035, 308.0883136, 306.6247203, 308.6184201, 308.6830487, 306.0353728,
    308.4789734, 308.097129,  307.4498101, 306.4571618, 306.3348016, 307.449898,  306.8918064, 306.4379689, 308.3057104,
    305.8244024, 307.7832466, 305.9537608, 308.4256222, 307.6612404, 305.9978287, 306.1627269, 307.4821939, 307.0098737,
    305.9258368, 308.0316925, 308.404896,  306.6772489, 307.1441669, 307.769201,  307.0541718, 306.3905305, 307.5552371,
    307.0951414, 307.311043,  308.6336743, 307.0834437, 306.7422114, 306.7223216, 308.2104013, 307.5611921, 306.5807061,
    307.7215757, 308.1577335, 307.6429405, 307.8988021, 308.1819402, 308.2631038, 307.4730703, 305.9205743, 308.0726969,
    306.2175604, 307.7846751, 305.6169576, 305.2916826, 306.4418186, 307.9873754, 306.3515944, 307.8798629, 306.8339679,
    307.6220998, 305.4361687, 306.584915,  306.5069071, 308.011441,  306.3209287, 305.6135625, 306.9757229, 306.7890305,
    306.916095,  306.5813547, 306.744411,  306.8335505, 306.5867868, 308.1741714, 306.0897144, 305.4680507, 307.8125832,
    307.3571521, 306.0402141, 307.184291,  307.1835549, 307.0372405, 307.0840199, 305.6962037, 305.2478807, 306.3991056,
    307.9055191, 305.8134413, 305.8399075, 305.5569301, 306.6018084, 304.9970244, 307.2292339, 307.7746954, 307.6803197,
    307.3820439, 307.9751343, 307.4492357, 305.1299735, 307.5105224, 306.8530443, 305.7680432, 307.3193381, 306.8631509,
    306.8185668, 307.0390632, 305.4182244, 307.700534,  306.1784029, 307.5712078, 306.1657809, 305.189449,  305.8711327,
    307.4817421, 305.3787481, 306.0472523, 306.2071925, 305.9580448, 307.9557773, 305.754064,  306.1769804, 307.3841333,
    307.2805372, 305.7400142, 305.5328498, 307.6712751, 307.7788903, 306.8111616, 306.1689617, 305.6462088, 305.859433,
    306.1031148, 306.9101377, 305.4886536, 307.6575571, 304.9825334, 305.8299527, 305.4030654, 306.3010376, 306.357146,
    307.057029,  306.4965379, 305.5599069, 305.5645026, 306.2943541, 307.1203741, 307.5704016, 307.8517386, 306.0768634,
    307.2010159, 305.5804168, 306.5364965, 306.3469362, 306.8309923, 306.6511972, 306.9782457, 305.2733407, 305.4728793,
    305.2912188, 307.8001396, 305.627772,  307.1668523, 306.4876935, 307.3750693, 305.7695079, 306.8056921, 307.0369189,
    306.5110383, 306.0620944, 307.4957863, 306.9752775, 306.5710425, 305.4477758, 305.4528515, 305.7233082, 306.3863738,
    305.7117551, 306.2012565, 305.927307,  307.4829772, 305.0606872, 307.7354487, 306.1107522, 306.1136096, 304.9820865,
    307.0581987, 307.050172,  305.0797664, 305.622729,  304.9810833, 306.5723329, 307.2196545, 307.2654846, 307.7160309,
    307.4282317, 306.9246549, 307.1567933, 304.9786722, 305.1201814, 306.8880356, 305.4620093, 305.0589237, 307.4431505,
    304.7908602, 305.8305925, 304.8857783, 307.2814716, 307.3583644, 305.4137979, 306.2091732, 304.8884508, 307.556866,
    307.3039035, 306.9252168, 305.1672539, 305.2869234, 306.2486871, 305.4460185, 305.2625659, 306.6338623, 307.4689719,
    305.0662185, 305.8718175, 306.7560277, 307.0427027, 305.0387176, 305.4424005, 305.6777791, 305.1410258, 306.6966423,
    304.8519419, 306.1619919, 305.3890482, 305.7128601, 306.9203805, 306.0250744, 305.2912094, 306.6868383, 304.846985,
    304.7523684, 307.6827671, 304.8480043, 306.0399757, 307.5736516, 306.2654723, 305.4405891, 304.9586484, 307.2184755,
    306.2588872, 306.7142365, 305.5502765, 305.8823466, 307.7645003, 306.8094696, 307.2397953, 307.0704621, 307.0551947,
    307.2942226, 307.4862859, 306.9189325, 305.3964363, 305.019566,  306.7367947, 307.0990798, 307.3439834, 305.4991304,
    305.0246621, 305.6974766, 307.3646664, 305.5134816, 307.1263851, 305.3737373, 306.3000549, 306.8270476, 307.1980537,
    305.1637352, 305.0388343, 306.6245695, 307.4167664, 306.1423806, 307.6053206, 306.1912851, 304.9290219, 307.7331187,
    306.0684097, 307.8110142, 305.9730081, 305.3432728, 307.4122775, 306.4808796, 306.9049017, 306.4334538, 307.0952438,
    306.3211153, 306.511188,  306.4861401, 306.7099309, 305.2563924, 306.5516842, 306.8426088, 307.2826085, 306.9847809,
    307.4888327, 307.5039328, 307.4879898, 305.4206951, 307.7067811, 307.4928402, 305.6940031, 306.9157713, 307.7426473,
    307.2470943, 305.7712094, 305.5839218, 307.1594446, 306.4227357, 307.5168061, 307.71615,   305.9652921, 307.1660465,
    307.6245405, 305.0945084, 306.1365997, 305.3019045, 306.5104383, 304.8732433, 306.5101915, 307.0313022, 306.9501966,
    304.8972697, 306.0793411, 305.1392824, 305.8013785, 306.9592517, 305.8477576, 305.6857825, 305.0814681, 307.5763435,
    307.0051189, 305.0225708, 305.2497508, 306.5764733, 306.3426292, 305.5495827, 306.8298108, 306.8085155, 306.2967501,
    307.6348655, 305.0412606, 306.7269981, 307.5228473, 306.6590694, 305.8149035, 306.1341824, 307.5910359, 307.311384,
    305.4218405, 306.4511726, 307.118042,  305.3233664, 307.7590667, 305.9716328, 306.8776062, 305.3084645, 306.1507973,
    305.0561959, 307.3528002, 306.7036441, 305.1209326, 306.9318865, 306.9471204, 305.1141925, 306.5349586, 306.1560218,
    306.8049769};
} // namespace gym
} // namespace dsyre
#endif // MISC_DATA_H