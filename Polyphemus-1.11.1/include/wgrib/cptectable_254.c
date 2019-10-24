#include "cnames.h"

const struct ParmTable parm_table_cptec_254[256] = {
   /* 0 */ {"var0", "undefined"},
   /* 1 */ {"PRES", "Pressure [hPa]"},
   /* 2 */ {"psnm", "Pressure reduced to MSL [hPa]"},
   /* 3 */ {"tsps", "Pressure tendency [Pa/s]"},
   /* 4 */ {"var4", "undefined"},
   /* 5 */ {"var5", "undefined"},
   /* 6 */ {"geop", "Geopotential [dam]"},
   /* 7 */ {"zgeo", "Geopotential height [gpm]"},
   /* 8 */ {"gzge", "Geometric height [m]"},
   /* 9 */ {"var9", "undefined"},
   /* 10 */ {"var10", "undefined"},
   /* 11 */ {"temp", "ABSOLUTE TEMPERATURE [K]"},
   /* 12 */ {"vtmp", "VIRTUAL TEMPERATURE [K]"},
   /* 13 */ {"ptmp", "POTENTIAL TEMPERATURE [K]"},
   /* 14 */ {"psat", "PSEUDO-ADIABATIC POTENTIAL TEMPERATURE [K]"},
   /* 15 */ {"mxtp", "MAXIMUM TEMPERATURE [K]"},
   /* 16 */ {"mntp", "MINIMUM TEMPERATURE [K]"},
   /* 17 */ {"tpor", "DEW POINT TEMPERATURE [K]"},
   /* 18 */ {"dptd", "DEW POINT DEPRESSION [K]"},
   /* 19 */ {"lpsr", "LAPSE RATE [K/m]"},
   /* 20 */ {"var20", "undefined"},
   /* 21 */ {"rds1", "RADAR SPECTRA(1) [non-dim]"},
   /* 22 */ {"rds2", "RADAR SPECTRA(2) [non-dim]"},
   /* 23 */ {"rds3", "RADAR SPECTRA(3) [non-dim]"},
   /* 24 */ {"var24", "undefined"},
   /* 25 */ {"tpan", "TEMPERATURE ANOMALY [K]"},
   /* 26 */ {"psan", "PRESSURE ANOMALY [Pa hPa]"},
   /* 27 */ {"zgan", "GEOPOT HEIGHT ANOMALY [m]"},
   /* 28 */ {"wvs1", "WAVE SPECTRA(1) [non-dim]"},
   /* 29 */ {"wvs2", "WAVE SPECTRA(2) [non-dim]"},
   /* 30 */ {"wvs3", "WAVE SPECTRA(3) [non-dim]"},
   /* 31 */ {"wind", "WIND DIRECTION  [deg]"},
   /* 32 */ {"wins", "WIND SPEED [m/s]"},
   /* 33 */ {"uvel", "ZONAL WIND (U) [m/s]"},
   /* 34 */ {"vvel", "MERIDIONAL WIND (V) [m/s]"},
   /* 35 */ {"fcor", "STREAM FUNCTION [m2/s]"},
   /* 36 */ {"potv", "VELOCITY POTENTIAL [m2/s]"},
   /* 37 */ {"var37", "undefined"},
   /* 38 */ {"sgvv", "SIGMA COORD VERT VEL [sec/sec]"},
   /* 39 */ {"omeg", "OMEGA [Pa/s]"},
   /* 40 */ {"omg2", "VERTICAL VELOCITY [m/s]"},
   /* 41 */ {"abvo", "ABSOLUTE VORTICITY        [10**5/sec]"},
   /* 42 */ {"abdv", "ABSOLUTE DIVERGENCE [10**5/sec]"},
   /* 43 */ {"vort", "VORTICITY  [1/s]"},
   /* 44 */ {"divg", "DIVERGENCE [1/s]"},
   /* 45 */ {"vucs", "VERTICAL U-COMP SHEAR [1/sec]"},
   /* 46 */ {"vvcs", "VERT V-COMP SHEAR [1/sec]"},
   /* 47 */ {"dirc", "DIRECTION OF CURRENT [deg]"},
   /* 48 */ {"spdc", "SPEED OF CURRENT [m/s]"},
   /* 49 */ {"ucpc", "U-COMPONENT OF CURRENT [m/s]"},
   /* 50 */ {"vcpc", "V-COMPONENT OF CURRENT [m/s]"},
   /* 51 */ {"umes", "SPECIFIC HUMIDITY [kg/kg]"},
   /* 52 */ {"umrl", "RELATIVE HUMIDITY [no Dim]"},
   /* 53 */ {"hmxr", "HUMIDITY MIXING RATIO [kg/kg]"},
   /* 54 */ {"agpl", "INST. PRECIPITABLE WATER [Kg/m2]"},
   /* 55 */ {"vapp", "VAPOUR PRESSURE [Pa hpa]"},
   /* 56 */ {"sadf", "SATURATION DEFICIT        [Pa hPa]"},
   /* 57 */ {"evap", "EVAPORATION [Kg/m2/day]"},
   /* 58 */ {"var58", "undefined"},
   /* 59 */ {"prcr", "PRECIPITATION RATE        [kg/m2/day]"},
   /* 60 */ {"thpb", "THUNDER PROBABILITY [%]"},
   /* 61 */ {"prec", "TOTAL PRECIPITATION [Kg/m2/day]"},
   /* 62 */ {"prge", "LARGE SCALE PRECIPITATION [Kg/m2/day]"},
   /* 63 */ {"prcv", "CONVECTIVE PRECIPITATION [Kg/m2/day]"},
   /* 64 */ {"neve", "SNOWFALL [Kg/m2/day]"},
   /* 65 */ {"wenv", "WAT EQUIV ACC SNOW DEPTH [kg/m2]"},
   /* 66 */ {"nvde", "SNOW DEPTH        [cm]"},
   /* 67 */ {"mxld", "MIXED LAYER DEPTH [m cm]"},
   /* 68 */ {"tthd", "TRANS THERMOCLINE DEPTH [m cm]"},
   /* 69 */ {"mthd", "MAIN THERMOCLINE DEPTH [m cm]"},
   /* 70 */ {"mtha", "MAIN THERMOCLINE ANOM [m cm]"},
   /* 71 */ {"cbnv", "CLOUD COVER [0-1]"},
   /* 72 */ {"cvnv", "CONVECTIVE CLOUD COVER [0-1]"},
   /* 73 */ {"lwnv", "LOW CLOUD COVER [0-1]"},
   /* 74 */ {"mdnv", "MEDIUM CLOUD COVER        [0-1]"},
   /* 75 */ {"hinv", "HIGH CLOUD COVER [0-1]"},
   /* 76 */ {"wtnv", "CLOUD WATER [kg/m2]"},
   /* 77 */ {"bli", "BEST LIFTED INDEX (TO 500 HPA) [K]"},
   /* 78 */ {"var78", "undefined"},
   /* 79 */ {"var79", "undefined"},
   /* 80 */ {"var80", "undefined"},
   /* 81 */ {"lsmk", "LAND SEA MASK [0,1]"},
   /* 82 */ {"dslm", "DEV SEA_LEV FROM MEAN [m]"},
   /* 83 */ {"zorl", "ROUGHNESS LENGTH [m]"},
   /* 84 */ {"albe", "ALBEDO [%]"},
   /* 85 */ {"dstp", "DEEP SOIL TEMPERATURE [K]"},
   /* 86 */ {"soic", "SOIL MOISTURE CONTENT [Kg/m2]"},
   /* 87 */ {"vege", "VEGETATION        [%]"},
   /* 88 */ {"var88", "undefined"},
   /* 89 */ {"dens", "DENSITY [kg/m3]"},
   /* 90 */ {"var90", "Undefined"},
   /* 91 */ {"icec", "ICE CONCENTRATION [fraction]"},
   /* 92 */ {"icet", "ICE THICKNESS [m]"},
   /* 93 */ {"iced", "DIRECTION OF ICE DRIFT [deg]"},
   /* 94 */ {"ices", "SPEED OF ICE DRIFT        [m/s]"},
   /* 95 */ {"iceu", "U-COMP OF ICE DRIFT [m/s]"},
   /* 96 */ {"icev", "V-COMP OF ICE DRIFT [m/s]"},
   /* 97 */ {"iceg", "ICE GROWTH        [m]"},
   /* 98 */ {"icdv", "ICE DIVERGENCE [sec/sec]"},
   /* 99 */ {"var99", "undefined"},
   /* 100 */ {"shcw", "SIG HGT COM WAVE/SWELL [m]"},
   /* 101 */ {"wwdi", "DIRECTION OF WIND WAVE [deg]"},
   /* 102 */ {"wwsh", "SIG HGHT OF WIND WAVES [m]"},
   /* 103 */ {"wwmp", "MEAN PERIOD WIND WAVES [sec]"},
   /* 104 */ {"swdi", "DIRECTION OF SWELL WAVE [deg]"},
   /* 105 */ {"swsh", "SIG HEIGHT SWELL WAVES [m]"},
   /* 106 */ {"swmp", "MEAN PERIOD SWELL WAVES [sec]"},
   /* 107 */ {"prwd", "PRIMARY WAVE DIRECTION [deg]"},
   /* 108 */ {"prmp", "PRIM WAVE MEAN PERIOD [s]"},
   /* 109 */ {"swdi", "SECOND WAVE DIRECTION [deg]"},
   /* 110 */ {"swmp", "SECOND WAVE MEAN PERIOD [s]"},
   /* 111 */ {"ocas", "SHORT WAVE ABSORBED AT GROUND [W/m2]"},
   /* 112 */ {"slds", "NET LONG WAVE AT BOTTOM [W/m2]"},
   /* 113 */ {"nswr", "NET SHORT-WAV RAD(TOP) [W/m2]"},
   /* 114 */ {"role", "OUTGOING LONG WAVE AT TOP [W/m2]"},
   /* 115 */ {"lwrd", "LONG-WAV RAD [W/m2]"},
   /* 116 */ {"swea", "SHORT WAVE ABSORBED BY EARTH/ATMOSPHERE  [W/m2]"},
   /* 117 */ {"glbr", "GLOBAL RADIATION [W/m2 ]"},
   /* 118 */ {"var118", "undefined"},
   /* 119 */ {"var119", "undefined"},
   /* 120 */ {"var120", "undefined"},
   /* 121 */ {"clsf", "LATENT HEAT FLUX FROM SURFACE [W/m2]"},
   /* 122 */ {"cssf", "SENSIBLE HEAT FLUX FROM SURFACE [W/m2]"},
   /* 123 */ {"blds", "BOUND LAYER DISSIPATION [W/m2]"},
   /* 124 */ {"var124", "undefined"},
   /* 125 */ {"var125", "undefined"},
   /* 126 */ {"var126", "undefined"},
   /* 127 */ {"imag", "IMAGE [image^data]"},
   /* 128 */ {"tp2m", "2 METRE TEMPERATURE [K]"},
   /* 129 */ {"dp2m", "2 METRE DEWPOINT TEMPERATURE [K]"},
   /* 130 */ {"u10m", "10 METRE U-WIND COMPONENT [m/s]"},
   /* 131 */ {"v10m", "10 METRE V-WIND COMPONENT [m/s]"},
   /* 132 */ {"topo", "TOPOGRAPHY [m]"},
   /* 133 */ {"gsfp", "GEOMETRIC MEAN SURFACE PRESSURE [hPa]"},
   /* 134 */ {"lnsp", "LN SURFACE PRESSURE [hPa]"},
   /* 135 */ {"pslc", "SURFACE PRESSURE [hPa]"},
   /* 136 */ {"pslm", "M S L PRESSURE (MESINGER METHOD) [hPa]"},
   /* 137 */ {"mask", "MASK  [-/+]"},
   /* 138 */ {"mxwu", "MAXIMUM U-WIND [m/s]"},
   /* 139 */ {"mxwv", "MAXIMUM V-WIND [m/s]"},
   /* 140 */ {"cape", "CONVECTIVE AVAIL. POT.ENERGY [m2/s2]"},
   /* 141 */ {"cine", "CONVECTIVE INHIB. ENERGY [m2/s2]"},
   /* 142 */ {"lhcv", "CONVECTIVE LATENT HEATING [K/s]"},
   /* 143 */ {"mscv", "CONVECTIVE MOISTURE SOURCE [1/s]"},
   /* 144 */ {"scvm", "SHALLOW CONV. MOISTURE SOURCE [1/s]"},
   /* 145 */ {"scvh", "SHALLOW CONVECTIVE HEATING [K/s]"},
   /* 146 */ {"mxwp", "MAXIMUM WIND PRESS. LVL  [hPa]"},
   /* 147 */ {"ustr", "STORM MOTION U-COMPONENT [m/s]"},
   /* 148 */ {"vstr", "STORM MOTION V-COMPONENT [m/s]"},
   /* 149 */ {"cbnt", "MEAN CLOUD COVER [0-1]"},
   /* 150 */ {"pcbs", "PRESSURE AT CLOUD BASE [hPa]"},
   /* 151 */ {"pctp", "PRESSURE AT CLOUD TOP [hPa]"},
   /* 152 */ {"fzht", "FREEZING LEVEL HEIGHT [m]"},
   /* 153 */ {"fzrh", "FREEZING LEVEL RELATIVE HUMIDITY [%]"},
   /* 154 */ {"fdlt", "FLIGHT LEVELS TEMPERATURE [K]"},
   /* 155 */ {"fdlu", "FLIGHT LEVELS U-WIND [m/s]"},
   /* 156 */ {"fdlv", "FLIGHT LEVELS V-WIND [m/s]"},
   /* 157 */ {"tppp", "TROPOPAUSE PRESSURE   [hPa]"},
   /* 158 */ {"tppt", "TROPOPAUSE TEMPERATURE [K]"},
   /* 159 */ {"tppu", "TROPOPAUSE U-WIND COMPONENT [m/s]"},
   /* 160 */ {"tppv", "TROPOPAUSE v-WIND COMPONENT [m/s]"},
   /* 161 */ {"var161", "undefined"},
   /* 162 */ {"gvdu", "GRAVITY WAVE DRAG DU/DT [m/s2]"},
   /* 163 */ {"gvdv", "GRAVITY WAVE DRAG DV/DT [m/s2]"},
   /* 164 */ {"gvus", "GRAVITY WAVE DRAG SFC ZONAL STRESS  [Pa]"},
   /* 165 */ {"gvvs", "GRAVITY WAVE DRAG SFC MERIDIONAL STRESS [Pa]"},
   /* 166 */ {"var166", "undefined"},
   /* 167 */ {"dvsh", "DIVERGENCE OF SPECIFIC HUMIDITY [1/s]"},
   /* 168 */ {"hmfc", "HORIZ. MOISTURE FLUX CONV.  [1/s]"},
   /* 169 */ {"vmfl", "VERT. INTEGRATED MOISTURE FLUX CONV. [kg/(m2*s)]"},
   /* 170 */ {"vadv", "VERTICAL MOISTURE ADVECTION  [kg/(kg*s)]"},
   /* 171 */ {"nhcm", "NEG. HUM. CORR. MOISTURE SOURCE [kg/(kg*s)]"},
   /* 172 */ {"lglh", "LARGE SCALE LATENT HEATING   [K/s]"},
   /* 173 */ {"lgms", "LARGE SCALE MOISTURE SOURCE  [1/s]"},
   /* 174 */ {"smav", "SOIL MOISTURE AVAILABILITY  [0-1]"},
   /* 175 */ {"tgrz", "SOIL TEMPERATURE OF ROOT ZONE [K]"},
   /* 176 */ {"bslh", "BARE SOIL LATENT HEAT [Ws/m2]"},
   /* 177 */ {"evpp", "POTENTIAL SFC EVAPORATION [m]"},
   /* 178 */ {"rnof", "RUNOFF [kg/m2/s)]"},
   /* 179 */ {"pitp", "INTERCEPTION LOSS [W/m2]"},
   /* 180 */ {"vpca", "VAPOR PRESSURE OF CANOPY AIR SPACE [mb]"},
   /* 181 */ {"qsfc", "SURFACE SPEC HUMIDITY   [kg/kg]"},
   /* 182 */ {"ussl", "SOIL WETNESS OF SURFACE [0-1]"},
   /* 183 */ {"uzrs", "SOIL WETNESS OF ROOT ZONE [0-1]"},
   /* 184 */ {"uzds", "SOIL WETNESS OF DRAINAGE ZONE [0-1]"},
   /* 185 */ {"amdl", "STORAGE ON CANOPY [m]"},
   /* 186 */ {"amsl", "STORAGE ON GROUND [m]"},
   /* 187 */ {"tsfc", "SURFACE TEMPERATURE [K]"},
   /* 188 */ {"tems", "SURFACE ABSOLUTE TEMPERATURE [K]"},
   /* 189 */ {"tcas", "TEMPERATURE OF CANOPY AIR SPACE [K]"},
   /* 190 */ {"ctmp", "TEMPERATURE AT CANOPY [K]"},
   /* 191 */ {"tgsc", "GROUND/SURFACE COVER TEMPERATURE [K]"},
   /* 192 */ {"uves", "SURFACE ZONAL WIND (U) [m/s]"},
   /* 193 */ {"usst", "SURFACE ZONAL WIND STRESS [Pa]"},
   /* 194 */ {"vves", "SURFACE MERIDIONAL WIND (V) [m/s]"},
   /* 195 */ {"vsst", "SURFACE MERIDIONAL WIND STRESS [Pa]"},
   /* 196 */ {"suvf", "SURFACE MOMENTUM FLUX [W/m2]"},
   /* 197 */ {"iswf", "INCIDENT SHORT WAVE FLUX [W/m2]"},
   /* 198 */ {"ghfl", "TIME AVE GROUND HT FLX   [W/m2]"},
   /* 199 */ {"var199", "undefined"},
   /* 200 */ {"lwbc", "NET LONG WAVE AT BOTTOM (CLEAR) [W/m2]"},
   /* 201 */ {"lwtc", "OUTGOING LONG WAVE AT TOP (CLEAR) [W/m2]"},
   /* 202 */ {"swec", "SHORT WV ABSRBD BY EARTH/ATMOS (CLEAR) [W/m2]"},
   /* 203 */ {"ocac", "SHORT WAVE ABSORBED AT GROUND (CLEAR) [W/m2]"},
   /* 204 */ {"var204", "undefined"},
   /* 205 */ {"lwrh", "LONG WAVE RADIATIVE HEATING  [K/s]"},
   /* 206 */ {"swrh", "SHORT WAVE RADIATIVE HEATING [K/s]"},
   /* 207 */ {"olis", "DOWNWARD LONG WAVE AT BOTTOM [W/m2]"},
   /* 208 */ {"olic", "DOWNWARD LONG WAVE AT BOTTOM (CLEAR) [W/m2]"},
   /* 209 */ {"ocis", "DOWNWARD SHORT WAVE AT GROUND [W/m2]"},
   /* 210 */ {"ocic", "DOWNWARD SHORT WAVE AT GROUND (CLEAR) [W/m2]"},
   /* 211 */ {"oles", "UPWARD LONG WAVE AT BOTTOM [W/m2]"},
   /* 212 */ {"oces", "UPWARD SHORT WAVE AT GROUND [W/m2]"},
   /* 213 */ {"swgc", "UPWARD SHORT WAVE AT GROUND (CLEAR) [W/m2]"},
   /* 214 */ {"roce", "UPWARD SHORT WAVE AT TOP [W/m2]"},
   /* 215 */ {"swtc", "UPWARD SHORT WAVE AT TOP (CLEAR) [W/m2]"},
   /* 216 */ {"var216", "undefined"},
   /* 217 */ {"var217", "undefined"},
   /* 218 */ {"hhdf", "HORIZONTAL HEATING DIFFUSION [K/s]"},
   /* 219 */ {"hmdf", "HORIZONTAL MOISTURE DIFFUSION [1/s]"},
   /* 220 */ {"hddf", "HORIZONTAL DIVERGENCE DIFFUSION [1/s2]"},
   /* 221 */ {"hvdf", "HORIZONTAL VORTICITY DIFFUSION [1/s2]"},
   /* 222 */ {"vdms", "VERTICAL DIFF. MOISTURE SOURCE [1/s]"},
   /* 223 */ {"vdfu", "VERTICAL DIFFUSION DU/DT [m/s2]"},
   /* 224 */ {"vdfv", "VERTICAL DIFFUSION DV/DT [m/s2]"},
   /* 225 */ {"vdfh", "VERTICAL DIFFUSION HEATING [K/s]"},
   /* 226 */ {"umrs", "SURFACE RELATIVE HUMIDITY [no Dim]"},
   /* 227 */ {"vdcc", "VERTICAL DIST TOTAL CLOUD COVER [no Dim]"},
   /* 228 */ {"var228", "undefined"},
   /* 229 */ {"var229", "undefined"},
   /* 230 */ {"usmt", "TIME MEAN SURFACE ZONAL WIND (U) [m/s]"},
   /* 231 */ {"vsmt", "TIME MEAN SURFACE MERIDIONAL WIND (V) [m/s]"},
   /* 232 */ {"tsmt", "TIME MEAN SURFACE ABSOLUTE TEMPERATURE [K]"},
   /* 233 */ {"rsmt", "TIME MEAN SURFACE RELATIVE HUMIDITY [no Dim]"},
   /* 234 */ {"atmt", "TIME MEAN ABSOLUTE TEMPERATURE [K]"},
   /* 235 */ {"stmt", "TIME MEAN DEEP SOIL TEMPERATURE [K]"},
   /* 236 */ {"ommt", "TIME MEAN DERIVED OMEGA [Pa/s]"},
   /* 237 */ {"dvmt", "TIME MEAN DIVERGENCE [1/s]"},
   /* 238 */ {"zhmt", "TIME MEAN GEOPOTENTIAL HEIGHT [m]"},
   /* 239 */ {"lnmt", "TIME MEAN LOG SURFACE PRESSURE [ln(cbar)]"},
   /* 240 */ {"mkmt", "TIME MEAN MASK [-/+]"},
   /* 241 */ {"vvmt", "TIME MEAN MERIDIONAL WIND (V) [m/s]"},
   /* 242 */ {"omtm", "TIME MEAN OMEGA  [cbar/s]"},
   /* 243 */ {"ptmt", "TIME MEAN POTENTIAL TEMPERATURE [K]"},
   /* 244 */ {"pcmt", "TIME MEAN PRECIP. WATER  [kg/m2]"},
   /* 245 */ {"rhmt", "TIME MEAN RELATIVE HUMIDITY [%]"},
   /* 246 */ {"mpmt", "TIME MEAN SEA LEVEL PRESSURE [hPa]"},
   /* 247 */ {"simt", "TIME MEAN SIGMADOT [1/s]"},
   /* 248 */ {"uemt", "TIME MEAN SPECIFIC HUMIDITY [kg/kg]"},
   /* 249 */ {"fcmt", "TIME MEAN STREAM FUNCTION| m2/s]"},
   /* 250 */ {"psmt", "TIME MEAN SURFACE PRESSURE [hPa]"},
   /* 251 */ {"tmmt", "TIME MEAN SURFACE TEMPERATURE [K]"},
   /* 252 */ {"pvmt", "TIME MEAN VELOCITY POTENTIAL [m2/s]"},
   /* 253 */ {"tvmt", "TIME MEAN VIRTUAL TEMPERATURE [K]"},
   /* 254 */ {"vtmt", "TIME MEAN VORTICITY [1/s]"},
   /* 255 */ {"uvmt", "TIME MEAN ZONAL WIND (U) [m/s]"},
};

