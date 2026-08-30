// PROBE: is a species' cohort set really a ONE-PARAMETER FAMILY IN HEIGHT that a
// boundary could answer with a sweep, instead of the per-cohort recording plant
// pays for now?
//
// docs/design/one-order.md, "The one-parameter family in height", claims twelve of
// the leaf's sixteen scalar inputs are traits and only three channels move with the
// cohort: kmax (proportional to 1/height), ppfd (the canopy read at the crown), and
// the 2L root resistances (through rooting_depth = min(height, rooting_depth_max)).
// It leaves the deciding question unmeasured: how does the ANSWER move with height?
//
// This reproduces plant's height -> leaf-input map for the century fixture
// (scripts/refined-century.rds) and asks, at a held soil and a held canopy:
//
//   1. the operating point (opt_root_psi_, opt_psi_stem_, ci_, profit_) and its
//      KIND, over a realistic height range, densely;
//   2. how many operating-point KIND switches happen over that range, and how far
//      the stand's own soil sits from the nearest switch;
//   3. what a table of K heights costs in accuracy -- linear and cubic
//      interpolation error against exact solves, so a node count can be read off;
//   4. the same for the two DERIVATIVE rows a boundary would have to carry,
//      dprofit*/dkmax and dprofit*/dppfd, taken at FReal<double>. A boundary that
//      interpolates the value and not the row is no boundary at all.
//
//   make -C tests/cpp CXX=g++ ODELIA_INC=../../../odelia/inst/include
//        ODELIA_SRC=../../../odelia/src probe_height && ./probe_height

#include <phylloptim.hpp>

#include "leaf_inputs.hpp"
#include "root_network.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace pl = phylloptim;
using T = pl::tangent;  // FReal<double>

namespace {

// ---------------------------------------------------------------------------
// THE CENTURY FIXTURE'S OWN NUMBERS.
// Dumped from scripts/refined-century.rds (one TF24 species, lma 0.0825,
// hmat 5.13, k_I 0.5, a_l1 5.44, a_l2 0.306) and from Environment("TF24").
// Written out rather than re-derived: a probe that recomputes the hyperpar chain
// is measuring its own copy of it.
// ---------------------------------------------------------------------------
constexpr double kEta = 12.0;                       // canopy shape
constexpr double kK_s = 1.0;                        // max hydraulic conductivity
constexpr double kThetaHuber = 0.00021417862497322766;  // sapwood per leaf area
constexpr double kPPFD = 1800.0;                    // above-canopy PPFD
constexpr double kK_I = 0.5;                        // self-shading coefficient
constexpr double kA_r1 = 0.07;                      // root mass per leaf area
constexpr double kRootMassCarbonScale = 83.26 * 0.5;
constexpr double kRootDepthShapeEta = 0.2;
constexpr double kRootingDepthMax = 1.5;            // metres

constexpr double kVpd = 1.0, kCa = 40.0, kLeafTemp = 25.0, kO2 = 21.0,
                 kAtm = 101.3;

// The twelve traits + R_d_25, byte-identical across every cohort and every stage.
constexpr double kVcmax25 = 96.0;
constexpr double kStemC = 2.04;
constexpr double kStemB = 3.4572679073200265;
constexpr double kPsiCrit = 5.9198804427473473;
constexpr double kRootC = 2.6801469999999998;
constexpr double kRootB = 3.8982450000000002;
constexpr double kRootPsiCrit = 5.8702825428827037;
constexpr double kBeta2 = 1.5;
constexpr double kJmax25 = 157.44;
constexpr double kQuantumYield = 0.29999999999999999;
constexpr double kCurvElec = 0.69999999999999996;
constexpr double kCurvColim = 0.98999999999999999;
constexpr double kCostScale = 7.5;
constexpr double kR_d_25 = 1.4399999999999999;

// plant's Control for this run.
constexpr double kGSSTolAbs = 0.1;
constexpr double kNControl = 1600.0;
constexpr double kCiAbsTol = 0.001;
constexpr double kCiNiter = 1000.0;

// The soil column: five layers over 1.5 m.
const std::vector<double> kSoilDepth = {0.3, 0.6, 0.9, 1.2, 1.5};
constexpr double kSoilMoistSat = 0.428, kA_psi = 1780.0, kN_psi = 6.57;
constexpr double kSoilMoistResidual = 1e-2;

// The soil the century stand ACTUALLY reached at t = 105.32 (its last step):
// layer moistures 0.1770, 0.1495, 0.1495, 0.1495, 0.1495.
const std::vector<double> kStandMoist = {0.1770139, 0.14949996, 0.14950129,
                                         0.14952376, 0.1495291};

double psi_from_moist(double m) {
  const double t = std::max(m, kSoilMoistResidual);
  return kA_psi * std::pow(t / kSoilMoistSat, -kN_psi) / 1e6;  // Pa -> MPa
}

std::vector<double> psi_soil_from_moist(const std::vector<double>& m) {
  std::vector<double> out;
  for (double x : m) out.push_back(psi_from_moist(x));
  return out;
}

// A whole soil column at one moisture, for the dryness sweep.
std::vector<double> uniform_psi(double moist) {
  return std::vector<double>(kSoilDepth.size(), psi_from_moist(moist));
}

double eta_c() { return 1.0 - 2.0 / (1.0 + kEta) + 1.0 / (1.0 + 2.0 * kEta); }

// ---------------------------------------------------------------------------
// THE CANOPY the stand actually cast at t = 105.32, as plant's own light
// interpolator holds it: 65 knots of (height, openness, slope), read back as a
// cubic Hermite so the ppfd channel this sweep sees is the one the stand saw --
// including its knot structure, which a smooth analytic stand-in would hide.
// ---------------------------------------------------------------------------
struct Knot { double z, y, dy; };
const Knot kCanopy[] = {
  {0, 0.18830394517088384, -0},
  {0.092198484603845901, 0.188303945171371, 6.340160860566147e-11},
  {0.1843969692076918, 0.18830394716605187, 1.2983256779235365e-07},
  {0.27659545381153772, 0.18830420224880282, 1.1075071577270442e-05},
  {0.3687939384153836, 0.18831032212757579, 0.00014906555017977509},
  {0.46099242301922949, 0.18831983961286217, 0.00014255459704089104},
  {0.55319090762307543, 0.18832686974973331, 0.000107732405544714},
  {0.64538939222692127, 0.18833863203481063, 0.00010217279553978952},
  {0.73758787683076721, 0.18835562455808791, 0.00024849597853614642},
  {0.82978636143461315, 0.18838228152344985, 0.00042907181555137925},
  {0.92198484603845898, 0.18843343884076469, 0.00074523091181930598},
  {1.0141833306423049, 0.18852132062995455, 0.0010380891902707196},
  {1.1063818152461509, 0.18864714356608933, 0.001508207510075995},
  {1.1985802998499968, 0.18879881506536117, 0.0019622597719328974},
  {1.2907787844538425, 0.18900937620206248, 0.0027783290226103137},
  {1.3829772690576885, 0.18932160572345988, 0.0040140459431507129},
  {1.4751757536615344, 0.18976061171564604, 0.0054875395543967269},
  {1.5673742382653804, 0.19032685127526394, 0.0069164575181494643},
  {1.6595727228692263, 0.19103686512174936, 0.0085036452955813361},
  {1.751771207473072, 0.19189376331546545, 0.010046922752953102},
  {1.843969692076918, 0.19288149470972382, 0.011308079890686499},
  {1.9361681766807639, 0.19398097103190132, 0.012523260716135852},
  {2.0283666612846099, 0.19516707980666936, 0.013331935451400706},
  {2.1205651458884556, 0.19642509271974837, 0.013734124531521586},
  {2.2127636304923017, 0.19769980121376143, 0.014179132789249267},
  {2.3049621150961475, 0.19902127550519363, 0.014048528244194702},
  {2.3971605996999936, 0.20024963815120375, 0.012048184439382232},
  {2.4893590843038393, 0.20133145745285377, 0.011262830990611646},
  {2.5815575689076851, 0.20236938169907065, 0.011065336839229655},
  {2.6737560535115312, 0.20341164243326648, 0.012064959182265382},
  {2.7659545381153769, 0.20465819774907387, 0.014909299401239891},
  {2.8581530227192231, 0.20612692866574525, 0.016719049616477171},
  {2.9503515073230688, 0.20768577694473978, 0.017038798053785807},
  {3.0425499919269146, 0.20925485457957971, 0.017096854140222564},
  {3.1347484765307607, 0.21084835185051767, 0.018142468624078895},
  {3.2269469611346064, 0.21263257381168077, 0.02019047970410881},
  {3.3191454457384526, 0.21446333102486267, 0.019490518474628415},
  {3.4113439303422983, 0.21634116151692345, 0.020777580852542197},
  {3.503542414946144, 0.21824243854654565, 0.02131036694981538},
  {3.5957408995499902, 0.2202524557241409, 0.021706841869565927},
  {3.6879393841538359, 0.22222335428224047, 0.022221177045765456},
  {3.7801378687576821, 0.22433665882699474, 0.023135244170011584},
  {3.8723363533615278, 0.22653565636438527, 0.025791785750459199},
  {3.9645348379653735, 0.2291080041451021, 0.029904850860109557},
  {4.0567333225692197, 0.23217465741173487, 0.036728302337223365},
  {4.1489318071730654, 0.2359079710619085, 0.044720731268894295},
  {4.2411302917769111, 0.24047232204824814, 0.054430053872202899},
  {4.3333287763807578, 0.24604860457687758, 0.067109019606375608},
  {4.4255272609846035, 0.25290826917096731, 0.08230765582753119},
  {4.5177257455884492, 0.26143801373099307, 0.10342773354238372},
  {4.6099242301922949, 0.27211930469837076, 0.12928624190278429},
  {4.7021227147961406, 0.2856032403263728, 0.16460929201807878},
  {4.7943211993999872, 0.30276224085107845, 0.20936350180330304},
  {4.886519684003833, 0.32461858514165581, 0.26718566736386401},
  {4.9787181686076787, 0.35251347909026354, 0.34081507496223312},
  {5.0709166532115244, 0.38804663846489312, 0.43337657897571003},
  {5.1631151378153701, 0.43307967621626248, 0.54736958549678161},
  {5.2553136224192167, 0.4897325828312375, 0.68485512660706116},
  {5.3475121070230625, 0.55977680726494561, 0.83511688278875196},
  {5.4397105916269082, 0.6432669447346977, 0.96922260687439032},
  {5.5319090762307539, 0.73702075577442727, 1.053015863205873},
  {5.6241075608345996, 0.83416843373396821, 1.0305634270918751},
  {5.7163060454384462, 0.92098305424695825, 0.82256391798452011},
  {5.8085045300422919, 0.98032715311692609, 0.4410309926488078},
  {5.9007030146461377, 1, -0},
};
const int kNCanopy = static_cast<int>(sizeof(kCanopy) / sizeof(kCanopy[0]));

double canopy_openness(double z) {
  if (z <= kCanopy[0].z) return kCanopy[0].y;
  if (z >= kCanopy[kNCanopy - 1].z) return 1.0;
  int lo = 0, hi = kNCanopy - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (kCanopy[mid].z <= z) lo = mid; else hi = mid;
  }
  const double h = kCanopy[hi].z - kCanopy[lo].z;
  const double t = (z - kCanopy[lo].z) / h, t2 = t * t, t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * kCanopy[lo].y +
         (t3 - 2 * t2 + t) * h * kCanopy[lo].dy +
         (-2 * t3 + 3 * t2) * kCanopy[hi].y +
         (t3 - t2) * h * kCanopy[hi].dy;
}

// The crown's leaf-area density, plant's CanopyShape::q_from_height at eta = 12.
double q_shape(double z, double height) {
  if (z <= 0.0) return 0.0;
  const double u = z / height;
  const double ue = std::pow(u, kEta);
  return 2.0 * kEta * (1.0 - ue) * ue / z;
}

// plant's default TF24 shading model is mean-light: one optimisation at the
// leaf-area-weighted MEAN openness over the crown. Its adaptive rule is replaced
// here by a converged 400-point composite Simpson, because what is being measured
// is the smoothness of the model and not the smoothness of a quadrature.
double mean_openness(double height) {
  const int n = 400;  // even
  const double h = height / n;
  double s = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double z = i * h;
    const double f = std::max(canopy_openness(z), 1e-4) * q_shape(z, height);
    const double w = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
    s += w * f;
  }
  return s * h / 3.0;
}


// ---------------------------------------------------------------------------
// THE STAND'S OWN COHORTS at its last step (t = 105.32): 170 of them, sorted by
// height, each with the operating point PLANT ITSELF recorded in the aux slots.
// Two jobs: it says what the height spacing actually is, and it referees this
// probe's reconstruction of the height -> leaf-input map against the model's own
// answer. A probe that cannot reproduce the numbers plant wrote down is measuring
// its own arithmetic.
// ---------------------------------------------------------------------------
struct Cohort { double height, root_psi, psi_stem, profit; };
const Cohort kCohorts[] = {
  {0.39204578130703827, 1.5248509014651659, 1.5562185920015232, 5.3675035032780247},
  {0.48626199804756715, 1.5244314741531917, 1.5630645053286913, 5.3530886716027624},
  {0.62634711009065214, 1.5253408633892582, 1.5745725785277886, 5.3294197736256779},
  {0.75989145328444208, 1.5265637118146911, 1.5856873705900587, 5.3065482477714099},
  {0.88544102161376981, 1.5274769490975775, 1.5957190609971925, 5.2857268571157769},
  {1.0068560654940442, 1.5281027611240245, 1.6050101528095426, 5.2668482077120231},
  {1.1232636662862541, 1.5285420075276515, 1.6136328902875539, 5.2503508304885704},
  {1.2222514498343569, 1.5287618403295029, 1.6207207714702592, 5.2376154750541426},
  {1.3030170419568885, 1.5288875079063473, 1.6263955200411515, 5.2281257225132123},
  {1.3754779766229377, 1.5289987895340671, 1.6314509102969044, 5.2205506995750053},
  {1.4396562419147558, 1.5291177660057977, 1.63592762062602, 5.2148408493731679},
  {1.4978023238941236, 1.5292592710072341, 1.6400047531289967, 5.210643993164032},
  {1.550718269587569, 1.52747284320805, 1.6418784975377849, 5.2111303380551446},
  {1.5975346217007123, 1.525903427258219, 1.6435509748536119, 5.2123651191474405},
  {1.643617380226694, 1.5244416582314519, 1.64528247694281, 5.214230706360957},
  {1.6860849847038886, 1.5231709462660183, 1.6469580970028777, 5.2165566541307218},
  {1.7264138787119967, 1.5220323429402798, 1.6486215959896675, 5.2193111531468865},
  {1.7666259444573376, 1.5209645236720173, 1.6503529180287793, 5.2226037782260413},
  {1.8032593819506935, 1.5200496702059847, 1.6519931755544712, 5.2260749785964817},
  {1.8401071651855729, 1.5191845407142779, 1.6537035401902873, 5.230017411561426},
  {1.8783850942237783, 1.5183433724406856, 1.6555441066413614, 5.2345864031771541},
  {1.9146665754240959, 1.5175976423934268, 1.6573462745087317, 5.2393423716824925},
  {1.9528424059892928, 1.516864552572863, 1.6593005706497272, 5.2447729387855455},
  {1.9944420081166696, 1.5161235921377361, 1.6614958076186426, 5.2511704706934594},
  {2.0359251610994691, 1.5154418000144569, 1.6637502256334216, 5.2580247467066785},
  {2.0801239073160427, 1.5147725723631602, 1.6662181043229916, 5.2658026775995195},
  {2.1334102860871802, 1.5140374233355742, 1.669276635129195, 5.2757747174259713},
  {2.1913236177518187, 1.5133169554761068, 1.6726926779912796, 5.2872620503533874},
  {2.2495435199048202, 1.5126559328419729, 1.6762010708826414, 5.2993220967139214},
  {2.3191424678887662, 1.511942861333412, 1.6804870055347272, 5.3143644767485503},
  {2.3992373999976819, 1.511190873030102, 1.68550098821139, 5.3321998629291096},
  {2.4750716918099092, 1.5104855078538506, 1.6902532283122307, 5.3490405538321726},
  {2.5351443610414242, 1.5098965036442769, 1.6939769729369536, 5.3620230816379166},
  {2.5883419620554426, 1.5093278787452848, 1.6972123767426659, 5.3730209358964265},
  {2.6375090921673907, 1.5087781511278076, 1.7001693934311746, 5.3829082334899487},
  {2.6773908685260079, 1.5083129713320989, 1.7025413817787298, 5.3907105824231003},
  {2.7066831803960474, 1.5079574198106953, 1.7042645007712625, 5.3962890020355587},
  {2.7343132963955687, 1.5076217695486822, 1.7058885393984449, 5.4015318744559293},
  {2.7490800470509344, 1.5074427247397897, 1.7067565624668093, 5.4043306153506867},
  {2.7644435868116704, 1.5072570983258557, 1.7076602379199524, 5.407243869482885},
  {2.7805050030047544, 1.5070639823985053, 1.7086058868464815, 5.4102931828731569},
  {2.7974912223247199, 1.5068615290478475, 1.7096079754809046, 5.4135292329864111},
  {2.8158929790131562, 1.5066462083162295, 1.7106984198116548, 5.4170664140066069},
  {2.8364470596613693, 1.5064120158988867, 1.7119242130620638, 5.4210693894758553},
  {2.888625300294394, 1.5058526134027281, 1.7150801220972898, 5.4315291565612629},
  {2.9240099372787793, 1.5055085797835539, 1.7172654476407745, 5.4389309940295103},
  {2.9684291155734082, 1.5051130132735446, 1.7200555698523381, 5.4485389578762993},
  {3.0251454627653334, 1.5046726678295483, 1.7237027651714731, 5.461380268108857},
  {3.0980759611043402, 1.5041984686074383, 1.72851468382485, 5.4787055170525791},
  {3.3048522653670118, 1.5033370258681442, 1.7428176230070518, 5.5321144422579254},
  {3.4551501369574735, 1.5031778067326735, 1.7538749503794073, 5.5751968800875034},
  {3.6320058011715024, 1.5034132189590808, 1.7675184553452277, 5.629822496013638},
  {3.8081329680953946, 1.5039675434102358, 1.7816226351188016, 5.6872148975472445},
  {3.9633550239111925, 1.504729586163601, 1.7945105586242318, 5.7404284677594628},
  {4.1085770572090192, 1.5059142054342858, 1.8073210380469042, 5.7949474994480559},
  {4.2543838200114283, 1.5080294715020628, 1.8216331537997235, 5.8591983464455115},
  {4.4148984981195882, 1.512097477807238, 1.840134858935297, 5.9481319901843781},
  {4.6029527902478975, 1.520591925124591, 1.8677890986524239, 6.0920603017942012},
  {4.8038060608802846, 1.5363892795457756, 1.9083936986338002, 6.3190429952092826},
  {5.150295894206522, 1.5913425399709298, 2.0261861992605295, 7.0226113714502176},
  {5.4556668229819776, 1.684055947791772, 2.2133058456495842, 8.1627051105106592},
  {5.5148598762633432, 1.7071296979606565, 2.2605204213301024, 8.4449930315121016},
  {5.5637029414290602, 1.7270646428276755, 2.3018208085973662, 8.688312267750252},
  {5.6461922067652424, 1.7618518833586083, 2.3753387599313198, 9.1104229312799045},
  {5.6681603074387983, 1.771202802745145, 2.3954756131639887, 9.2230212596757646},
  {5.6857231397370587, 1.7786650582122461, 2.4116758495209227, 9.3125085986309788},
  {5.7135517886820573, 1.7904232171130181, 2.4374561356758511, 9.4527171485522938},
  {5.7376409400019561, 1.8004938837076412, 2.4598039115515582, 9.5718820919698082},
  {5.7590975713642161, 1.809348241695115, 2.479673791913835, 9.6758179402737863},
  {5.7784597711793939, 1.8172151508228489, 2.4975172939539605, 9.767410369615078},
  {5.796263678461127, 1.8243245648878619, 2.5138089732800251, 9.8494996232496437},
  {5.8136623217520205, 1.8311405384106183, 2.5295880772084263, 9.9275251281571162},
  {5.8310090145399887, 1.8377937902120709, 2.5451524734433923, 10.002986617852034},
  {5.8417187797296926, 1.8418231186062475, 2.5546622429968378, 10.04832568638227},
  {5.8448034803936562, 1.8429719120707349, 2.5573857218965559, 10.061199833730335},
  {5.8468172667533267, 1.8437189725076442, 2.5591597709674221, 10.069559125384995},
  {5.8483870083046625, 1.8442996959414464, 2.5605404532156832, 10.076050158103296},
  {5.8499294464245191, 1.8448689358194268, 2.5618952323583613, 10.082406831684763},
  {5.8553911474774498, 1.8468734353562257, 2.5666770821507758, 10.104743101511492},
  {5.8562760233133764, 1.8471965388304703, 2.5674495153680494, 10.108336408128604},
  {5.8573440365118001, 1.8475858923526436, 2.5683809469315038, 10.112663866851962},
  {5.8602563216134236, 1.8486441197312002, 2.5709159246079203, 10.124410887125743},
  {5.87119950281856, 1.8525742464186297, 2.5803755382065758, 10.167847043518943},
  {5.8836367589286276, 1.8569494847643537, 2.5909937299263603, 10.21583679421583},
  {5.8871711395110582, 1.8581745314681728, 2.5939840276096415, 10.229202337423381},
  {5.8890678171195212, 1.8588284979558087, 2.5955835428456893, 10.236324247209858},
  {5.8904323305480375, 1.859297482565148, 2.596732006989086, 10.241426021392495},
  {5.8915433936671286, 1.8596784307792573, 2.5976657453996737, 10.245566630648154},
  {5.8932977132033733, 1.8602782368499764, 2.5991374981967952, 10.25207970758551},
  {5.8947283860900521, 1.8607658484848175, 2.6003353864700554, 10.257368756064832},
  {5.8959291741983142, 1.8611740399967593, 2.6013391598974867, 10.261792360874058},
  {5.8969580428846493, 1.8615230123258171, 2.6021980301011429, 10.265571309931493},
  {5.8974280842376388, 1.8616822021675967, 2.6025900405458771, 10.267294255368761},
  {5.8978660288973064, 1.8618303867070414, 2.6029550753092274, 10.268897586745954},
  {5.8982782181026057, 1.861969737459203, 2.6032984595220205, 10.270404896358791},
  {5.8986687360373526, 1.8621016551211642, 2.6036236261383956, 10.271831409721148},
  {5.8990407981317929, 1.8622272418498824, 2.603933276975595, 10.273189105006537},
  {5.8992210778970824, 1.8622880600256571, 2.6040832637164084, 10.273846473787984},
  {5.8993976238528898, 1.862347597139113, 2.6042301110442883, 10.274489916517165},
  {5.8995706949734306, 1.8624059418205596, 2.6043740363547259, 10.275120395908587},
  {5.8997405799025513, 1.8624631925497372, 2.6045152814824095, 10.27573898053978},
  {5.899907544183721, 1.8625194398848555, 2.6046540688533595, 10.276346652819054},
  {5.8999900128034293, 1.8625472151162443, 2.6047226092817235, 10.276646698687474},
  {5.9000718162105121, 1.8625747617328967, 2.6047905898025978, 10.276944257983942},
  {5.9001529769091823, 1.8626020874216842, 2.6048580292841925, 10.277239414150406},
  {5.900233499615287, 1.8626291938760324, 2.6049249318070395, 10.277532185874609},
  {5.9003134011812683, 1.8626560868720827, 2.6049913115310561, 10.277822635930281},
  {5.9003531191839134, 1.8626694534043509, 2.6050243055141764, 10.277966991326736},
  {5.9003926830777731, 1.8626827670056261, 2.6050571698325529, 10.278110771130141},
  {5.9004320929102914, 1.8626960277042715, 2.605089904544609, 10.278253975693291},
  {5.9004713481493143, 1.8627092353335979, 2.6051225092272952, 10.278396603262491},
  {5.9005104475721151, 1.8627223894946072, 2.6051549828840468, 10.278538649575404},
  {5.900529938536617, 1.8627289464221264, 2.6051711703213574, 10.278609453669912},
  {5.9005493905781599, 1.8627354899976818, 2.6051873250351494, 10.278680112630983},
  {5.9005688031239112, 1.8627420200301226, 2.6052034465520664, 10.278750624400388},
  {5.9005881758947298, 1.8627485364271215, 2.6052195346427025, 10.278820987987007},
  {5.9006075077873348, 1.8627550388191529, 2.6052355883932652, 10.278891199406692},
  {5.9006171574597079, 1.8627582844461028, 2.6052436016073948, 10.278926244630748},
  {5.9006267952718101, 1.8627615260205928, 2.6052516048749994, 10.278961245863584},
  {5.9006364193226597, 1.8627647629035029, 2.6052595966178931, 10.278996196205121},
  {5.9006460284580919, 1.8627679947069287, 2.6052675758780919, 10.279031091468049},
  {5.9006556181377663, 1.8627712199041739, 2.6052755388858984, 10.279065915169147},
  {5.9006604046468754, 1.8627728296775732, 2.6052795134360602, 10.279083296428027},
  {5.9006651845153817, 1.8627744372020663, 2.6052834824480837, 10.279100653347138},
  {5.9006699566052818, 1.8627760420950024, 2.605287444977133, 10.279117981795121},
  {5.9006747201075589, 1.8627776440843804, 2.605291400351553, 10.279135278835424},
  {5.9006794732993519, 1.8627792425908849, 2.6052953471409088, 10.279152538213392},
  {5.9006818454783589, 1.8627800403529742, 2.6052973168591707, 10.279161151780713},
  {5.9006842143748033, 1.8627808370073131, 2.6052992838458837, 10.279169753373338},
  {5.900686579984848, 1.862781632552629, 2.6053012480978901, 10.279178342977593},
  {5.9006889417574699, 1.8627824268036326, 2.6053032091576833, 10.279186918592913},
  {5.9006912992518723, 1.8627832196121055, 2.6053051666592837, 10.279195478619084},
  {5.9006924762656361, 1.8627836154319835, 2.6053061439685741, 10.279199752317536},
  {5.9006936520481199, 1.862784010836849, 2.6053071202540403, 10.279204021531603},
  {5.9006948265379897, 1.8627844058060803, 2.6053080954647645, 10.279208286038591},
  {5.9006959996694972, 1.8627848003175698, 2.6053090695461565, 10.279212545599883},
  {5.9006971713728831, 1.8627851943478628, 2.6053100424403022, 10.279216799962198},
  {5.90069775666593, 1.8627853911747918, 2.6053105284229749, 10.279218925109847},
  {5.900698341573472, 1.8627855878718469, 2.6053110140851934, 10.279221048854406},
  {5.9006989260881575, 1.8627857844365563, 2.6053114994208562, 10.279223171169185},
  {5.9006995102158708, 1.8627859808709031, 2.6053119844348553, 10.279225292075505},
  {5.9007000939266723, 1.8627861771648107, 2.6053124691023157, 10.279227411464724},
  {5.9007003856223426, 1.8627862752579651, 2.6053127113032901, 10.279228470578094},
  {5.9007006772095476, 1.8627863733145844, 2.60531295341411, 10.279229529296822},
  {5.9007009686867971, 1.8627864713341693, 2.6053131954335416, 10.279230587615459},
  {5.9007012600526645, 1.8627865693162409, 2.6053134373604037, 10.279231645528832},
  {5.900701551305775, 1.862786667260335, 2.6053136791935492, 10.279232703031971},
  {5.9007016968895503, 1.8627867162179756, 2.6053138000745713, 10.2792332316279},
  {5.9007018424445112, 1.8627867651659105, 2.6053139209316418, 10.279233760119006},
  {5.9007019879704981, 1.8627868141040871, 2.6053140417646325, 10.2792342885047},
  {5.9007021334672665, 1.8627868630324242, 2.6053141625733418, 10.279234816784076},
  {5.9007022789346468, 1.8627869119508644, 2.6053142833576275, 10.279235344956545},
  {5.9007023516572588, 1.8627869364063545, 2.6053143437405657, 10.279235609002495},
  {5.9007024243724029, 1.8627869608593282, 2.6053144041172938, 10.279235873021253},
  {5.9007024970801139, 1.8627869853097991, 2.6053144644878459, 10.279236137012997},
  {5.9007025697803561, 1.862787009757755, 2.6053145248521918, 10.279236400977556},
  {5.9007026424731128, 1.8627870342031898, 2.6053145852103157, 10.279236664914894},
  {5.9007026788166748, 1.8627870464249581, 2.6053146153870355, 10.279236796873304},
  {5.9007027151583422, 1.8627870586460902, 2.6053146455621854, 10.279236928824828},
  {5.9007027514981241, 1.8627870708665859, 2.605314675735765, 10.279237060769496},
  {5.9007027878359866, 1.8627870830864355, 2.6053147059077499, 10.279237192707182},
  {5.9007028241720043, 1.8627870953056638, 2.6053147360782019, 10.279237324638149},
  {5.900702847984129, 1.8627871033133052, 2.6053147558498395, 10.279237411096627},
  {5.9007028717954491, 1.8627871113206749, 2.6053147756208075, 10.279237497552153},
  {5.9007028956059404, 1.8627871193277632, 2.6053147953910822, 10.279237584004704},
  {5.9007029194156209, 1.8627871273345833, 2.6053148151606926, 10.2792376704543},
  {5.9007029432244744, 1.8627871353411214, 2.6053148349296085, 10.279237756900846},
  {5.9007029670325224, 1.8627871433473893, 2.6053148546978577, 10.279237843344502},
  {5.9007029908397293, 1.8627871513533731, 2.6053148744654062, 10.279237929785074},
  {5.9007030146461377, 1.8627871593590895, 2.605314894232293, 10.27923801622274},
};
const int kNCohorts = static_cast<int>(sizeof(kCohorts) / sizeof(kCohorts[0]));

enum class Light { Held, CrownCentre, MeanLight };

const char* light_name(Light m) {
  if (m == Light::Held) return "held";
  if (m == Light::CrownCentre) return "crown-centre";
  return "mean-light";
}

// ---------------------------------------------------------------------------
// plant's height -> leaf-input map, the three channels the design names.
// ---------------------------------------------------------------------------
double kmax_at(double height) {
  return kK_s * kThetaHuber / (height * eta_c());
}

double ppfd_at(double height, Light mode, double held_openness) {
  double light;
  if (mode == Light::Held) {
    light = held_openness;
  } else if (mode == Light::CrownCentre) {
    light = canopy_openness(height * eta_c());
  } else {
    light = mean_openness(height);
  }
  return kK_I * std::max(light, 1e-4) * kPPFD;
}

// The root profile, exactly as tf24_strategy.h builds it: carbon per leaf area
// over the five layers, capped at rooting_depth_max.
double Q_root(double z, double rooting_depth) {
  if (z > rooting_depth) return 0.0;
  const double u = std::pow(z / rooting_depth, kRootDepthShapeEta);
  const double t = 1.0 - u;
  return t * t;
}

std::vector<double> root_carbon_at(double height) {
  const double rd = std::min(height, kRootingDepthMax);
  const double scale = kRootMassCarbonScale * kA_r1;
  std::vector<double> carbon(kSoilDepth.size(), 0.0);
  double prev = 1.0;
  for (std::size_t a = 0; a < kSoilDepth.size(); ++a) {
    if (prev == 0.0) break;
    const double q = Q_root(kSoilDepth[a], rd);
    carbon[a] = scale * (prev - q);
    prev = q;
  }
  return carbon;
}

// ---------------------------------------------------------------------------
// One placement.
// ---------------------------------------------------------------------------
struct Point {
  double height = 0, kmax = 0, ppfd = 0;
  double root_psi = 0, psi_stem = 0, ci = 0, profit = 0;
  double dprofit_dkmax = 0, dprofit_dppfd = 0;
  const char* kind = "unsolved";
  bool ok = false;
};

pl::Leaf make_leaf() {
  pl::Leaf l(kVcmax25, kStemC, kStemB, kPsiCrit, kRootC, kRootB, kRootPsiCrit,
             kBeta2, kJmax25, kQuantumYield, kCurvElec, kCurvColim, kGSSTolAbs,
             kNControl, kCiAbsTol, kCiNiter, kCostScale);
  l.R_d_25 = kR_d_25;
  return l;
}

void drive_at(pl::Leaf& l, double height, const std::vector<double>& psi_soil,
              Light mode, double held_openness) {
  const std::vector<double> carbon = root_carbon_at(height);
  l.set_physiology(fixture::root_network(carbon, kSoilDepth),
                   ppfd_at(height, mode, held_openness), psi_soil, kSoilDepth,
                   kmax_at(height), kVpd, kCa, kLeafTemp, kO2, kAtm);
}

// The two derivative rows a boundary would carry, at the TOTAL (collar free)
// profit -- the composition a graft actually needs, not the frozen-collar
// partial. Seeded one at a time at FReal<double>, exactly as probe_rank does.
void rows_at(pl::Leaf& l, double height, Point& p) {
  fixture::Soil soil;
  soil.depth = kSoilDepth;
  soil.carbon = root_carbon_at(height);
  const pl::LeafInputs<double> base =
      fixture::leaf_inputs<double>(l, fixture::Input::None, 0, soil);
  const fixture::Input which[2] = {fixture::Input::kmax, fixture::Input::ppfd};
  double out[2] = {0.0, 0.0};
  for (int k = 0; k < 2; ++k) {
    pl::LeafInputs<T> in = base.template rebind_from<T>();
    // Seed by rebuilding at T with the one direction on: leaf_inputs does the
    // seeding, and it is the same routine every other probe here uses.
    in = fixture::leaf_inputs<T>(l, which[k], 0, soil);
    const auto draw = l.supply_draw_at<T>(T(l.opt_root_psi_), in.supply);
    const T moved = l.collar_at<T>(in, draw);
    const pl::Leaf::LeafOutputs<T> o = l.outputs_at<T>(moved, in, draw);
    out[k] = pl::derivative_along(o.profit);
  }
  p.dprofit_dkmax = out[0];
  p.dprofit_dppfd = out[1];
}

Point place(pl::Leaf& l, double height, const std::vector<double>& psi_soil,
            Light mode, double held_openness, bool with_rows) {
  Point p;
  p.height = height;
  p.kmax = kmax_at(height);
  p.ppfd = ppfd_at(height, mode, held_openness);
  try {
    drive_at(l, height, psi_soil, mode, held_openness);
    l.find_root_collar_psi();
  } catch (const std::exception& e) {
    p.kind = "threw";
    return p;
  }
  p.kind = pl::Leaf::operating_point_kind_name(l.operating_point_kind());
  p.root_psi = l.opt_root_psi_;
  p.psi_stem = l.opt_psi_stem_;
  p.ci = l.ci_;
  p.profit = l.profit_;
  p.ok = true;
  if (with_rows) {
    try {
      rows_at(l, height, p);
    } catch (const std::exception& e) {
      p.dprofit_dkmax = std::nan("");
      p.dprofit_dppfd = std::nan("");
    }
  }
  return p;
}

// ---------------------------------------------------------------------------
// Sweeps.
// ---------------------------------------------------------------------------
std::vector<double> log_grid(double lo, double hi, int n) {
  std::vector<double> out;
  out.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    out[std::size_t(i)] = lo * std::pow(hi / lo, double(i) / (n - 1));
  }
  return out;
}

double field(const Point& p, int which) {
  if (which == 0) return p.root_psi;
  if (which == 1) return p.psi_stem;
  if (which == 2) return p.ci;
  if (which == 3) return p.profit;
  if (which == 4) return p.dprofit_dkmax;
  return p.dprofit_dppfd;
}
const char* field_name(int which) {
  static const char* n[6] = {"opt_root_psi", "opt_psi_stem", "ci", "profit",
                             "dprofit*/dkmax", "dprofit*/dppfd"};
  return n[which];
}

// How many times the KIND changes walking the sweep, and where.
int kind_switches(const std::vector<Point>& s, bool report) {
  int n = 0;
  for (std::size_t i = 1; i < s.size(); ++i) {
    if (std::string(s[i].kind) != std::string(s[i - 1].kind)) {
      ++n;
      if (report && n <= 8) {
        std::printf("      switch %s -> %s between h = %.5g and %.5g\n",
                    s[i - 1].kind, s[i].kind, s[i - 1].height, s[i].height);
      }
    }
  }
  return n;
}

void sweep_table(const std::vector<Point>& s, int every) {
  std::printf("  %10s %11s %9s %11s %11s %10s %10s %11s %11s  %s\n", "height",
              "kmax", "ppfd", "root_psi", "psi_stem", "ci", "profit",
              "dPi/dkmax", "dPi/dppfd", "kind");
  for (std::size_t i = 0; i < s.size(); i += std::size_t(every)) {
    const Point& p = s[i];
    std::printf("  %10.4f %11.4e %9.2f %11.6f %11.6f %10.4f %10.5f %11.4e %11.4e  %s\n",
                p.height, p.kmax, p.ppfd, p.root_psi, p.psi_stem, p.ci, p.profit,
                p.dprofit_dkmax, p.dprofit_dppfd, p.kind);
  }
}

// Smoothness, measured three ways and none of them a picture:
//   * the largest LOCAL slope d y / d ln h, which says how fast the answer moves
//     per doubling of height;
//   * the largest second difference relative to the field's own range, which is
//     what an interpolation error is built out of;
//   * whether the field is monotone in height at all.
void smoothness(const std::vector<Point>& s, int which) {
  double lo = 1e300, hi = -1e300, max_slope = 0.0, max_d2 = 0.0;
  bool up = true, down = true;
  for (const Point& p : s) {
    const double y = field(p, which);
    if (!std::isfinite(y)) continue;
    lo = std::min(lo, y);
    hi = std::max(hi, y);
  }
  for (std::size_t i = 1; i + 1 < s.size(); ++i) {
    const double y0 = field(s[i - 1], which), y1 = field(s[i], which),
                 y2 = field(s[i + 1], which);
    if (!(std::isfinite(y0) && std::isfinite(y1) && std::isfinite(y2))) continue;
    const double l0 = std::log(s[i - 1].height), l1 = std::log(s[i].height),
                 l2 = std::log(s[i + 1].height);
    max_slope = std::max(max_slope, std::abs((y2 - y0) / (l2 - l0)));
    // Second difference on the (nearly uniform in log h) grid, as a fraction of
    // the range: this is the term a linear interpolant's error is proportional to.
    const double w = (l1 - l0) / (l2 - l0);
    max_d2 = std::max(max_d2, std::abs(y1 - (y0 + w * (y2 - y0))));
    if (y1 < y0) up = false;
    if (y1 > y0) down = false;
  }
  const double range = hi - lo;
  std::printf("  %-16s range %11.5g   max |dy/dlnh| %11.5g   "
              "max local curvature %10.4g (%.2g of range)   %s\n",
              field_name(which), range, max_slope, max_d2,
              range > 0 ? max_d2 / range : 0.0,
              (up || down) ? "monotone" : "NOT monotone");
}

// ---------------------------------------------------------------------------
// The number the design turns on: a table of K heights, and the error of reading
// the operating point off it instead of solving.
// ---------------------------------------------------------------------------
struct InterpErr { double lin_abs = 0, lin_rel = 0, cub_abs = 0, cub_rel = 0; };

double lerp_at(const std::vector<double>& x, const std::vector<double>& y,
               double xq) {
  std::size_t j = 0;
  while (j + 2 < x.size() && x[j + 1] < xq) ++j;
  const double t = (xq - x[j]) / (x[j + 1] - x[j]);
  return y[j] + t * (y[j + 1] - y[j]);
}

// Catmull-Rom on a non-uniform grid, via the three-point slopes at the bracketing
// nodes -- a cubic Hermite, which is what a table of solves would actually carry
// if the boundary stored the row alongside the value.
double cubic_at(const std::vector<double>& x, const std::vector<double>& y,
                double xq) {
  const std::size_t n = x.size();
  std::size_t j = 0;
  while (j + 2 < n && x[j + 1] < xq) ++j;
  auto slope = [&](std::size_t i) {
    if (i == 0) return (y[1] - y[0]) / (x[1] - x[0]);
    if (i + 1 >= n) return (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
    return (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
  };
  const double h = x[j + 1] - x[j];
  const double t = (xq - x[j]) / h, t2 = t * t, t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * y[j] + (t3 - 2 * t2 + t) * h * slope(j) +
         (-2 * t3 + 3 * t2) * y[j + 1] + (t3 - t2) * h * slope(j + 1);
}

// The reference the tables are scored against, and the tables themselves.
struct Table {
  std::vector<double> x;
  std::vector<std::vector<double> > y;  // [field][node]
};

Table build_table(pl::Leaf& l, const std::vector<double>& xs,
                  const std::vector<double>& psi_soil, Light mode,
                  double held_openness, int nf) {
  Table t;
  t.x = xs;
  t.y.assign(std::size_t(nf), std::vector<double>());
  for (double h : xs) {
    const Point p = place(l, h, psi_soil, mode, held_openness, nf > 4);
    for (int f = 0; f < nf; ++f) t.y[std::size_t(f)].push_back(field(p, f));
  }
  return t;
}

void interp_test(pl::Leaf& l, const std::vector<double>& psi_soil, Light mode,
                 double held_openness, double h_lo, double h_hi,
                 const std::vector<int>& node_counts, bool with_rows) {
  const int nf = with_rows ? 6 : 4;
  const std::vector<double> hq = log_grid(h_lo, h_hi, 601);
  const Table ref = build_table(l, hq, psi_soil, mode, held_openness, nf);

  std::vector<double> range(std::size_t(nf), 1.0);
  for (int f = 0; f < nf; ++f) {
    double lo = 1e300, hi = -1e300;
    for (double y : ref.y[std::size_t(f)]) {
      if (std::isfinite(y)) { lo = std::min(lo, y); hi = std::max(hi, y); }
    }
    range[std::size_t(f)] = (hi > lo) ? hi - lo : 1.0;
  }

  std::printf("  %6s %10s | %-16s %11s %11s | %11s %11s %6s\n", "nodes",
              "spacing", "field", "linear max", "linear rel", "cubic max",
              "cubic rel", "order");
  std::vector<double> prev_cub(std::size_t(nf), 0.0);
  bool have_prev = false;
  for (std::size_t ki = 0; ki < node_counts.size(); ++ki) {
    const int K = node_counts[ki];
    const std::vector<double> xs = log_grid(h_lo, h_hi, K);
    const Table tab = build_table(l, xs, psi_soil, mode, held_openness, nf);
    std::vector<double> gaps;
    for (std::size_t i = 1; i < xs.size(); ++i) gaps.push_back(xs[i] - xs[i - 1]);
    std::sort(gaps.begin(), gaps.end());
    const double med_gap = gaps[gaps.size() / 2];

    for (int f = 0; f < nf; ++f) {
      double lin = 0.0, cub = 0.0;
      for (std::size_t i = 0; i < hq.size(); ++i) {
        const double y = ref.y[std::size_t(f)][i];
        if (!std::isfinite(y)) continue;
        lin = std::max(lin, std::abs(lerp_at(xs, tab.y[std::size_t(f)], hq[i]) - y));
        cub = std::max(cub, std::abs(cubic_at(xs, tab.y[std::size_t(f)], hq[i]) - y));
      }
      const double ord = (have_prev && cub > 0 && prev_cub[std::size_t(f)] > 0)
                             ? std::log(prev_cub[std::size_t(f)] / cub) / std::log(2.0)
                             : 0.0;
      std::printf("  %6d %10.4g | %-16s %11.4e %11.2e | %11.4e %11.2e %6.2f\n",
                  K, med_gap, field_name(f), lin, lin / range[std::size_t(f)],
                  cub, cub / range[std::size_t(f)], ord);
      prev_cub[std::size_t(f)] = cub;
    }
    have_prev = true;
  }
}

// ---------------------------------------------------------------------------
// AGAINST plant's OWN RECORDING, at the stand's own 170 cohorts.
// ---------------------------------------------------------------------------
void against_the_stand(pl::Leaf& l, const std::vector<double>& psi_soil,
                       Light mode, const std::vector<int>& node_counts) {
  std::printf("\n=== 4. AGAINST plant's OWN RECORDING, AT THE STAND'S 170 COHORTS ===\n");
  std::vector<double> hs, plant_rp, plant_ps, plant_pf;
  for (int i = 0; i < kNCohorts; ++i) {
    hs.push_back(kCohorts[i].height);
    plant_rp.push_back(kCohorts[i].root_psi);
    plant_ps.push_back(kCohorts[i].psi_stem);
    plant_pf.push_back(kCohorts[i].profit);
  }
  // The gaps the stand actually has.
  std::vector<double> gaps;
  for (int i = 1; i < kNCohorts; ++i) gaps.push_back(hs[std::size_t(i)] - hs[std::size_t(i - 1)]);
  std::sort(gaps.begin(), gaps.end());
  std::printf("  height range %.4f to %.4f m; adjacent spacing min %.3e "
              "median %.3e max %.3e m\n",
              hs[0], hs[std::size_t(kNCohorts - 1)], gaps[0],
              gaps[gaps.size() / 2], gaps.back());

  // This probe's own solve at those heights, against what plant wrote down.
  std::vector<Point> mine;
  int switches = 0;
  for (int i = 0; i < kNCohorts; ++i) {
    mine.push_back(place(l, hs[std::size_t(i)], psi_soil, mode, 1.0, false));
    if (i && std::string(mine[std::size_t(i)].kind) !=
                 std::string(mine[std::size_t(i - 1)].kind)) ++switches;
  }
  double e_rp = 0, e_ps = 0, e_pf = 0;
  for (int i = 0; i < kNCohorts; ++i) {
    e_rp = std::max(e_rp, std::abs(mine[std::size_t(i)].root_psi - plant_rp[std::size_t(i)]));
    e_ps = std::max(e_ps, std::abs(mine[std::size_t(i)].psi_stem - plant_ps[std::size_t(i)]));
    e_pf = std::max(e_pf, std::abs(mine[std::size_t(i)].profit - plant_pf[std::size_t(i)]));
  }
  std::printf("  this probe vs plant's aux, worst over the 170: opt_root_psi %.3e, "
              "opt_psi_stem %.3e, profit %.3e\n", e_rp, e_ps, e_pf);
  std::printf("  operating-point KIND switches between ADJACENT cohorts: %d of %d pairs\n",
              switches, kNCohorts - 1);

  // What plant's own recording says the collar does across one cohort gap.
  double dmax = 0, dsum = 0;
  std::vector<double> ds;
  for (int i = 1; i < kNCohorts; ++i) {
    const double d = std::abs(plant_rp[std::size_t(i)] - plant_rp[std::size_t(i - 1)]);
    ds.push_back(d);
    dmax = std::max(dmax, d);
    dsum += d;
  }
  std::sort(ds.begin(), ds.end());
  std::printf("  |d opt_root_psi| across ONE cohort gap (plant's own numbers): "
              "median %.3e, 95th %.3e, max %.3e MPa\n",
              ds[ds.size() / 2], ds[std::size_t(0.95 * ds.size())], dmax);
  std::printf("  total variation of opt_root_psi across all 170: %.4f MPa "
              "(its range is %.4f)\n", dsum,
              *std::max_element(plant_rp.begin(), plant_rp.end()) -
                  *std::min_element(plant_rp.begin(), plant_rp.end()));

  // And the deciding number: a K-node table read at those 170 heights.
  std::printf("\n  A K-node table (log-spaced over the stand's own height range),\n"
              "  read at the 170 cohorts, against this probe's exact solve there:\n");
  std::printf("  %6s %10s | %11s %11s | %11s %11s\n", "nodes", "spacing",
              "lin root_psi", "lin profit", "cub root_psi", "cub profit");
  for (int K : node_counts) {
    const std::vector<double> xs =
        log_grid(hs[0], hs[std::size_t(kNCohorts - 1)], K);
    const Table tab = build_table(l, xs, psi_soil, mode, 1.0, 4);
    double lr = 0, lp = 0, cr = 0, cp = 0;
    for (int i = 0; i < kNCohorts; ++i) {
      const double h = hs[std::size_t(i)];
      lr = std::max(lr, std::abs(lerp_at(xs, tab.y[0], h) - mine[std::size_t(i)].root_psi));
      cr = std::max(cr, std::abs(cubic_at(xs, tab.y[0], h) - mine[std::size_t(i)].root_psi));
      lp = std::max(lp, std::abs(lerp_at(xs, tab.y[3], h) - mine[std::size_t(i)].profit));
      cp = std::max(cp, std::abs(cubic_at(xs, tab.y[3], h) - mine[std::size_t(i)].profit));
    }
    std::vector<double> g2;
    for (std::size_t i = 1; i < xs.size(); ++i) g2.push_back(xs[i] - xs[i - 1]);
    std::sort(g2.begin(), g2.end());
    std::printf("  %6d %10.4g | %11.4e %11.4e | %11.4e %11.4e\n", K,
                g2[g2.size() / 2], lr, lp, cr, cp);
  }
}

// ---------------------------------------------------------------------------
// How far the stand sits from a kind boundary. A boundary indexed by
// species x stage is only usable while every cohort answers with the same arm,
// so the question is not "is the century fixture interior" but "how much drier
// would it have to be before it is not".
// ---------------------------------------------------------------------------
void kind_map(pl::Leaf& l) {
  std::printf("\n=== 3. WHERE THE KIND BOUNDARIES ARE ===\n");
  std::printf("  Every cell is a full height sweep (81 heights, 0.3 to 7 m,\n"
              "  crown-centre light off the stand's own canopy).\n\n");
  std::printf("  %10s %11s | %s\n", "moisture", "psi_soil", "kinds over the height sweep");
  std::vector<double> moists = {0.30, 0.214, 0.1770, 0.160, 0.1495, 0.140,
                                0.130, 0.120, 0.110, 0.100, 0.090, 0.080,
                                0.070, 0.060, 0.050};
  // And a fine walk across the two transitions the coarse row straddles, because
  // the question a boundary turns on is not "is this stand interior" but "how
  // thin is the band of soils where ONE species' cohorts do not agree".
  for (int i = 0; i <= 40; ++i) moists.push_back(0.1240 + 0.0001 * i);
  std::sort(moists.begin(), moists.end(), std::greater<double>());
  for (double m : moists) {
    const std::vector<double> psi = uniform_psi(m);
    std::vector<Point> s;
    for (double h : log_grid(0.3, 7.0, 81)) {
      s.push_back(place(l, h, psi, Light::CrownCentre, 1.0, false));
    }
    // Census of kinds in this row.
    std::vector<std::string> seen;
    std::vector<int> n;
    for (const Point& p : s) {
      const std::string k = p.kind;
      std::size_t j = 0;
      for (; j < seen.size(); ++j) if (seen[j] == k) break;
      if (j == seen.size()) { seen.push_back(k); n.push_back(0); }
      ++n[j];
    }
    std::printf("  %10.4f %11.4f | ", m, psi[0]);
    for (std::size_t j = 0; j < seen.size(); ++j) {
      std::printf("%s x%d  ", seen[j].c_str(), n[j]);
    }
    std::printf(" (%d switches)\n", kind_switches(s, false));
  }
  std::printf("\n  The century stand's own last-step soil is moisture 0.1770 in the\n"
              "  top layer and 0.1495 below it; its psi_soil is 0.588 / 1.788 MPa.\n");
}

void one_sweep(pl::Leaf& l, const std::vector<double>& psi_soil, Light mode,
               double held_openness, double h_lo, double h_hi, int n,
               const char* label, bool table) {
  std::vector<Point> s;
  for (double h : log_grid(h_lo, h_hi, n)) {
    s.push_back(place(l, h, psi_soil, mode, held_openness, true));
  }
  std::printf("\n--- %s: light %s", label, light_name(mode));
  if (mode == Light::Held) std::printf(" (openness %.3f)", held_openness);
  std::printf(", %d heights over %.3g to %.3g m ---\n", n, h_lo, h_hi);
  if (table) sweep_table(s, n / 20);
  const int nsw = kind_switches(s, true);
  std::printf("    kind switches over the sweep: %d\n", nsw);
  for (int f = 0; f < 6; ++f) smoothness(s, f);
}


// ---------------------------------------------------------------------------
// 5. AN ADAPTIVE TABLE. A log-spaced grid spends its nodes evenly and the field
// does not need them evenly, so the fair question is not "how many uniform nodes"
// but "how many nodes IF THEY ARE PUT WHERE THE ERROR IS". Bisect any interval
// whose geometric midpoint the linear interpolant misses by more than tol, and
// count what is left. The build cost is counted too: every test is itself an
// exact solve, and a design that hides them is quoting half the price.
// ---------------------------------------------------------------------------
void adaptive_table(pl::Leaf& l, const std::vector<double>& psi_soil, Light mode,
                    double h_lo, double h_hi, int which, double tol,
                    int cap = 4000) {
  std::vector<double> x = log_grid(h_lo, h_hi, 5);
  std::vector<double> y;
  int solves = 0;
  for (double h : x) { y.push_back(field(place(l, h, psi_soil, mode, 1.0, false), which)); ++solves; }
  int rounds = 0;
  while (static_cast<int>(x.size()) < cap) {
    std::vector<double> nx, ny;
    bool split = false;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
      nx.push_back(x[i]);
      ny.push_back(y[i]);
      const double m = std::sqrt(x[i] * x[i + 1]);
      const double exact = field(place(l, m, psi_soil, mode, 1.0, false), which);
      ++solves;
      const double lin = y[i] + (m - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]);
      if (std::abs(exact - lin) > tol) {
        nx.push_back(m);
        ny.push_back(exact);
        split = true;
      }
    }
    nx.push_back(x.back());
    ny.push_back(y.back());
    x.swap(nx);
    y.swap(ny);
    ++rounds;
    if (!split) break;
  }
  std::printf("  %-16s tol %8.1e -> %5d nodes (%d exact solves to build, "
              "%d refinement rounds)\n",
              field_name(which), tol, static_cast<int>(x.size()), solves, rounds);
}

// ---------------------------------------------------------------------------
// 6. SUBSAMPLING THE COHORTS THEMSELVES. The stand's own heights are already
// where the answers are wanted, and they are not evenly spread -- 170 of them
// over 5.5 m with a median gap of 1.4 mm and a maximum of 35 cm. So: solve at
// every m-th cohort and read the rest off, which is the boundary a species x
// stage index would actually build.
// ---------------------------------------------------------------------------
void subsample_cohorts(pl::Leaf& l, const std::vector<double>& psi_soil,
                       Light mode) {
  std::printf("\n=== 6. SOLVING AT EVERY m-TH COHORT AND READING THE REST OFF ===\n");
  std::vector<double> hs;
  for (int i = 0; i < kNCohorts; ++i) hs.push_back(kCohorts[i].height);
  std::vector<Point> exact;
  for (double h : hs) exact.push_back(place(l, h, psi_soil, mode, 1.0, false));

  std::printf("  %4s %7s %11s | %11s %11s | %11s %11s\n", "m", "solves",
              "med spacing", "lin root_psi", "lin profit", "cub root_psi",
              "cub profit");
  for (int m : {2, 3, 4, 6, 8, 12, 16, 24}) {
    std::vector<double> xs;
    std::vector<double> yr, yp;
    for (int i = 0; i < kNCohorts; i += m) {
      xs.push_back(hs[std::size_t(i)]);
      yr.push_back(exact[std::size_t(i)].root_psi);
      yp.push_back(exact[std::size_t(i)].profit);
    }
    if (xs.back() != hs.back()) {
      xs.push_back(hs.back());
      yr.push_back(exact[std::size_t(kNCohorts - 1)].root_psi);
      yp.push_back(exact[std::size_t(kNCohorts - 1)].profit);
    }
    // Strictly increasing knots only: cohorts 2e-8 m apart make a zero-width
    // interval, and a division by it is a NaN rather than an answer.
    std::vector<double> ux, ur, up;
    for (std::size_t i = 0; i < xs.size(); ++i) {
      if (ux.empty() || xs[i] > ux.back() * (1.0 + 1e-12)) {
        ux.push_back(xs[i]); ur.push_back(yr[i]); up.push_back(yp[i]);
      }
    }
    if (ux.size() < 4) continue;
    double lr = 0, lp = 0, cr = 0, cp = 0;
    for (int i = 0; i < kNCohorts; ++i) {
      const double h = std::min(std::max(hs[std::size_t(i)], ux.front()), ux.back());
      lr = std::max(lr, std::abs(lerp_at(ux, ur, h) - exact[std::size_t(i)].root_psi));
      cr = std::max(cr, std::abs(cubic_at(ux, ur, h) - exact[std::size_t(i)].root_psi));
      lp = std::max(lp, std::abs(lerp_at(ux, up, h) - exact[std::size_t(i)].profit));
      cp = std::max(cp, std::abs(cubic_at(ux, up, h) - exact[std::size_t(i)].profit));
    }
    std::vector<double> g;
    for (std::size_t i = 1; i < ux.size(); ++i) g.push_back(ux[i] - ux[i - 1]);
    std::sort(g.begin(), g.end());
    std::printf("  %4d %7d %11.4g | %11.4e %11.4e | %11.4e %11.4e\n", m,
                static_cast<int>(ux.size()), g[g.size() / 2], lr, lp, cr, cp);
  }
}


// ---------------------------------------------------------------------------
// 7. THE SECOND STRUCTURAL AXIS, which the design's list of three channels does
// not name. rooting_depth = min(height, rooting_depth_max) does not merely SCALE
// the root profile: Q(z, rooting_depth) is exactly zero below the rooting depth,
// so a short cohort carries root carbon in FEWER LAYERS, and root_network_from_carbon
// truncates at the deepest non-zero one. The number of supply inputs -- and with
// it the taped shape -- is therefore piecewise constant in height with breaks at
// the soil layer boundaries.
//
// And the other half of the same fact is the prize: above rooting_depth_max the
// whole root profile stops depending on height, so the 2L resistances are not
// merely smooth in height, they are BYTE-IDENTICAL.
// ---------------------------------------------------------------------------
void layer_structure(pl::Leaf& l, const std::vector<double>& psi_soil) {
  std::printf("\n=== 7. HOW MANY SUPPLY INPUTS, AS A FUNCTION OF HEIGHT ===\n");
  std::printf("  %8s %8s %8s   %s\n", "height", "layers", "inputs", "note");
  int last = -1;
  std::vector<double> ref;
  for (double h : log_grid(0.15, 8.0, 900)) {
    const std::vector<double> carbon = root_carbon_at(h);
    const pl::RootNetwork net = fixture::root_network(carbon, kSoilDepth);
    const int n = static_cast<int>(net.r_R_H_min.size());
    if (n != last) {
      std::printf("  %8.4f %8d %8d   first height with this many\n", h, n,
                  12 + 1 + static_cast<int>(kSoilDepth.size()) + 2 * n + 2);
      last = n;
    }
  }
  // Byte-identity above the cap.
  const pl::RootNetwork a = fixture::root_network(root_carbon_at(1.5001), kSoilDepth);
  const pl::RootNetwork b = fixture::root_network(root_carbon_at(30.0), kSoilDepth);
  bool same = a.r_R_H_min.size() == b.r_R_H_min.size();
  for (std::size_t i = 0; same && i < a.r_R_H_min.size(); ++i) {
    same = (a.r_R_H_min[i] == b.r_R_H_min[i]) && (a.r_R_V_sum[i] == b.r_R_V_sum[i]);
  }
  std::printf("\n  resistances at h = 1.5001 m and h = 30 m are %s\n",
              same ? "BYTE-IDENTICAL" : "different");
  std::printf("  so above rooting_depth_max = %.1f m the ONLY cohort-varying leaf\n"
              "  inputs are kmax and ppfd: 2 scalars, not 2L + 2.\n", kRootingDepthMax);

  // Is the ANSWER kinked where the input count changes? Dense local sweeps.
  std::printf("\n  local sweeps across the four layer boundaries and the root cap,\n"
              "  81 heights over +/- 2%% of each, mean-light:\n");
  std::printf("  %8s | %14s %14s %14s\n", "about", "max curv psi",
              "max curv profit", "kinds");
  for (double h0 : {0.3, 0.6, 0.9, 1.2, 1.5}) {
    std::vector<Point> s;
    for (double h : log_grid(h0 * 0.98, h0 * 1.02, 81)) {
      s.push_back(place(l, h, psi_soil, Light::MeanLight, 1.0, true));
    }
    double c_psi = 0, c_pf = 0;
    for (std::size_t i = 1; i + 1 < s.size(); ++i) {
      c_psi = std::max(c_psi, std::abs(s[i].root_psi -
                                       0.5 * (s[i - 1].root_psi + s[i + 1].root_psi)));
      c_pf = std::max(c_pf, std::abs(s[i].profit -
                                     0.5 * (s[i - 1].profit + s[i + 1].profit)));
    }
    std::printf("  %8.2f | %14.4e %14.4e %14d\n", h0, c_psi, c_pf,
                kind_switches(s, false));
  }
}

}  // namespace

int main() {
  std::printf("does the leaf's answer move smoothly and slowly in HEIGHT alone?\n");
  std::printf("century fixture: one TF24 species, eta_c = %.10f, "
              "kmax(h) = %.6e / h\n", eta_c(), kK_s * kThetaHuber / eta_c());
  const std::vector<double> stand = psi_soil_from_moist(kStandMoist);
  std::printf("stand soil at t = 105.32: psi_soil =");
  for (double p : stand) std::printf(" %.4f", p);
  std::printf(" MPa\n");
  std::printf("canopy: openness %.4f at the ground, 1.0 at %.3f m\n",
              kCanopy[0].y, kCanopy[kNCanopy - 1].z);

  pl::Leaf l = make_leaf();

  std::printf("\n=== 1. THE OPERATING POINT ALONG HEIGHT ===\n");
  // The stand's own heights ran 0.392 to 5.901 m. Swept a little wider so the
  // rooting-depth cap at 1.5 m and the canopy top at 5.9 m are both inside.
  one_sweep(l, stand, Light::CrownCentre, 1.0, 0.3, 7.0, 401,
            "the stand's soil and canopy", true);
  one_sweep(l, stand, Light::MeanLight, 1.0, 0.3, 7.0, 401,
            "the same, plant's default mean-light", false);
  one_sweep(l, stand, Light::Held, 0.1883, 0.3, 7.0, 401,
            "light HELD: the kmax and rooting-depth channels alone", false);
  one_sweep(l, uniform_psi(0.11), Light::CrownCentre, 1.0, 0.3, 7.0, 401,
            "a much drier soil (psi_soil 12.6 MPa is past root_psi_crit)", false);

  std::printf("\n=== 2. WHAT A TABLE OF K HEIGHTS COSTS ===\n");
  std::printf("  601 exact solves over 0.3 to 7 m are the reference; the table is\n"
              "  K log-spaced solves, read back linearly and as a cubic Hermite.\n"
              "  'rel' is the max error as a fraction of the field's own range\n"
              "  over the sweep.\n\n");
  std::printf("  --- crown-centre light, 0.3 to 7 m ---\n");
  interp_test(l, stand, Light::CrownCentre, 1.0, 0.3, 7.0,
              {9, 17, 33, 65, 129, 257}, true);
  std::printf("\n  --- plant's own mean-light, over the stand's OWN height range"
              " 0.392 to 5.901 m ---\n");
  interp_test(l, stand, Light::MeanLight, 1.0, 0.39204578, 5.90070301,
              {9, 17, 33, 65, 129, 257}, true);

  kind_map(l);

  against_the_stand(l, stand, Light::MeanLight,
                    {9, 17, 33, 65, 129, 257, 513});

  std::printf("\n=== 5. AN ADAPTIVE TABLE OVER THE STAND'S HEIGHT RANGE ===\n");
  std::printf("  Bisect any interval the linear interpolant misses by more than\n"
              "  tol. Compare the node count against the stand's 170 cohorts.\n");
  for (double tol : {1e-2, 1e-3, 1e-4, 1e-5}) {
    adaptive_table(l, stand, Light::MeanLight, 0.39204578, 5.90070301, 0, tol);
  }
  for (double tol : {1e-2, 1e-3, 1e-4}) {
    adaptive_table(l, stand, Light::MeanLight, 0.39204578, 5.90070301, 3, tol);
  }

  subsample_cohorts(l, stand, Light::MeanLight);

  layer_structure(l, stand);

  std::printf("\ndone\n");
  return 0;
}
