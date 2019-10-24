#include "cnames.h"

const struct ParmTable parm_table_ecmwf_128[256] = {
      /* 0 */ {"var0", "undefined"},
      /* 1 */ {"STRF", "Stream function [m**2 s**-1]"},
      /* 2 */ {"VPOT", "Velocity potential [m**2 s**-1]"},
      /* 3 */ {"PT", "Potential temperature [K]"},
      /* 4 */ {"EQPT", "Equivalent potential temperature [K]"},
      /* 5 */ {"SEPT", "Saturated equivalent potential temperature [K]"},
      /* 6 */ {"SSFR", "Soil sand fraction [(0 - 1)]"},
      /* 7 */ {"SCFR", "Soil clay fraction [(0 - 1)]"},
      /* 8 */ {"SRO", "Surface runoff [m]"},
      /* 9 */ {"SSRO", "Sub-surface runoff [m]"},
      /* 10 */ {"WIND", "Wind speed [m s**-1]"},
      /* 11 */ {"UDVW", "U component of divergent wind [m s**-1]"},
      /* 12 */ {"VDVW", "V component of divergent wind [m s**-1]"},
      /* 13 */ {"URTW", "U component of rotational wind [m s**-1]"},
      /* 14 */ {"VRTW", "V component of rotational wind [m s**-1]"},
      /* 15 */ {"ALUVP", "UV visible albedo for direct radiation [(0 - 1)]"},
      /* 16 */ {"ALUVD", "UV visible albedo for diffuse radiation [(0 - 1)]"},
      /* 17 */ {"ALNIP", "Near IR albedo for direct radiation [(0 - 1)]"},
      /* 18 */ {"ALNID", "Near IR albedo for diffuse radiation [(0 - 1)]"},
      /* 19 */ {"UVCS", "Clear sky surface UV [W m**-2 s]"},
      /* 20 */ {"PARCS", "Clear sky surface PAR [W m**-2 s]"},
      /* 21 */ {"UCTP", "Unbalanced component of temperature [K]"},
      /* 22 */ {"UCLN", "Unbalanced component of logarithm of surface pressure []"},
      /* 23 */ {"UCDV", "Unbalanced component of divergence [s**-1]"},
      /* 24 */ {"var24", "Reserved for future unbalanced components []"},
      /* 25 */ {"var25", "Reserved for future unbalanced components []"},
      /* 26 */ {"CL", "Lake cover [(0 - 1)]"},
      /* 27 */ {"CVL", "Low vegetation cover [(0 - 1)]"},
      /* 28 */ {"CVH", "High vegetation cover [(0 - 1)]"},
      /* 29 */ {"TVL", "Type of low vegetation []"},
      /* 30 */ {"TVH", "Type of high vegetation []"},
      /* 31 */ {"CI", "Sea-ice cover [(0 - 1)]"},
      /* 32 */ {"ASN", "Snow albedo [(0 - 1)]"},
      /* 33 */ {"RSN", "Snow density [kg m**-3]"},
      /* 34 */ {"SSTK", "Sea surface temperature [K]"},
      /* 35 */ {"ISTL1", "Ice surface temperature layer 1 [K]"},
      /* 36 */ {"ISTL2", "Ice surface temperature layer 2 [K]"},
      /* 37 */ {"ISTL3", "Ice surface temperature layer 3 [K]"},
      /* 38 */ {"ISTL4", "Ice surface temperature layer 4 [K]"},
      /* 39 */ {"SWVL1", "Volumetric soil water layer 1 [m**3 m**-3]"},
      /* 40 */ {"SWVL2", "Volumetric soil water layer 2 [m**3 m**-3]"},
      /* 41 */ {"SWVL3", "Volumetric soil water layer 3 [m**3 m**-3]"},
      /* 42 */ {"SWVL4", "Volumetric soil water layer 4 [m**3 m**-3]"},
      /* 43 */ {"SLT", "Soil type []"},
      /* 44 */ {"ES", "Snow evaporation [m of water]"},
      /* 45 */ {"SMLT", "Snowmelt [m of water]"},
      /* 46 */ {"SDUR", "Solar duration [s]"},
      /* 47 */ {"DSRP", "Direct solar radiation [w m**-2]"},
      /* 48 */ {"MAGSS", "Magnitude of surface stress [N m**-2 s]"},
      /* 49 */ {"10FG", "10 metre wind gust [m s**-1]"},
      /* 50 */ {"LSPF", "Large-scale precipitation fraction [s]"},
      /* 51 */ {"MX2T24", "Maximum temperature at 2 metres since last 24 hours [K]"},
      /* 52 */ {"MN2T24", "Minimum temperature at 2 metres since last 24 hours [K]"},
      /* 53 */ {"MONT", "Montgomery potential [m**2 s**-2]"},
      /* 54 */ {"PRES", "Pressure [Pa]"},
      /* 55 */ {"MEAN2T24", "Mean temperature at 2 metres since last 24 hours [K]"},
      /* 56 */ {"MN2D24", "Mean 2 metre dewpoint temperature in past 24 hours [K]"},
      /* 57 */ {"UVB", "Downward UV radiation at the surface [w m**-2 s]"},
      /* 58 */ {"PAR", "Photosynthetically active radiation at the surface [w m**-2 s]"},
      /* 59 */ {"CAPE", "Convective available potential energy [J kg**-1]"},
      /* 60 */ {"PV", "Potential vorticity [K m**2 kg**-1 s**-1]"},
      /* 61 */ {"var61", "undefined"},
      /* 62 */ {"OBCT", "Observation count []"},
      /* 63 */ {"var63", "Start time for skin temperature difference [s]"},
      /* 64 */ {"var64", "Finish time for skin temperature difference [s]"},
      /* 65 */ {"var65", "Skin temperature difference [K]"},
      /* 66 */ {"var66", "Leaf area index, low vegetation [m**2 / m**2]"},
      /* 67 */ {"var67", "Leaf area index, high vegetation [m**2 / m**2]"},
      /* 68 */ {"var68", "Minimum stomatal resistance, low vegetation [s m**-1]"},
      /* 69 */ {"var69", "Minimum stomatal resistance, high vegetation [s m**-1]"},
      /* 70 */ {"var70", "Biome cover, low vegetation [(0 - 1)]"},
      /* 71 */ {"var71", "Biome cover, high vegetation [(0 - 1)]"},
      /* 72 */ {"ISSRD", "Instantaneous surface solar radiation downwards [w m**-2]"},
      /* 73 */ {"ISTRD", "Instantaneous surface thermal radiation downwards [w m**-2]"},
      /* 74 */ {"SDFOR", "Standard deviation of filtered subgrid orography [m]"},
      /* 75 */ {"CRWC", "Cloud rain water content [kg kg**-1]"},
      /* 76 */ {"CSWC", "Cloud snow water content [kg kg**-1]"},
      /* 77 */ {"ETADOT", "Eta-coordinate vertical velocity [s**-1]"},
      /* 78 */ {"TCLW", "Total column liquid water [kg m**-2]"},
      /* 79 */ {"TCIW", "Total column ice water [kg m**-2]"},
      /* 80 */ {"var80", "Experimental product []"},
      /* 81 */ {"var81", "Experimental product []"},
      /* 82 */ {"var82", "Experimental product []"},
      /* 83 */ {"var83", "Experimental product []"},
      /* 84 */ {"var84", "Experimental product []"},
      /* 85 */ {"var85", "Experimental product []"},
      /* 86 */ {"var86", "Experimental product []"},
      /* 87 */ {"var87", "Experimental product []"},
      /* 88 */ {"var88", "Experimental product []"},
      /* 89 */ {"var89", "Experimental product []"},
      /* 90 */ {"var90", "Experimental product []"},
      /* 91 */ {"var91", "Experimental product []"},
      /* 92 */ {"var92", "Experimental product []"},
      /* 93 */ {"var93", "Experimental product []"},
      /* 94 */ {"var94", "Experimental product []"},
      /* 95 */ {"var95", "Experimental product []"},
      /* 96 */ {"var96", "Experimental product []"},
      /* 97 */ {"var97", "Experimental product []"},
      /* 98 */ {"var98", "Experimental product []"},
      /* 99 */ {"var99", "Experimental product []"},
      /* 100 */ {"var100", "Experimental product []"},
      /* 101 */ {"var101", "Experimental product []"},
      /* 102 */ {"var102", "Experimental product []"},
      /* 103 */ {"var103", "Experimental product []"},
      /* 104 */ {"var104", "Experimental product []"},
      /* 105 */ {"var105", "Experimental product []"},
      /* 106 */ {"var106", "Experimental product []"},
      /* 107 */ {"var107", "Experimental product []"},
      /* 108 */ {"var108", "Experimental product []"},
      /* 109 */ {"var109", "Experimental product []"},
      /* 110 */ {"var110", "Experimental product []"},
      /* 111 */ {"var111", "Experimental product []"},
      /* 112 */ {"var112", "Experimental product []"},
      /* 113 */ {"var113", "Experimental product []"},
      /* 114 */ {"var114", "Experimental product []"},
      /* 115 */ {"var115", "Experimental product []"},
      /* 116 */ {"var116", "Experimental product []"},
      /* 117 */ {"var117", "Experimental product []"},
      /* 118 */ {"var118", "Experimental product []"},
      /* 119 */ {"var119", "Experimental product []"},
      /* 120 */ {"var120", "Experimental product []"},
      /* 121 */ {"MX2T6", "Maximum temperature at 2 metres since last 6 hours [K]"},
      /* 122 */ {"MN2T6", "Minimum temperature at 2 metres since last 6 hours [K]"},
      /* 123 */ {"10FG6", "10 metre wind gust in the past 6 hours [m s**-1]"},
      /* 124 */ {"EMIS", "Surface emissivity [dimensionless]"},
      /* 125 */ {"var125", "Vertically integrated total energy [J m**-2]"},
      /* 126 */ {"var126", "Generic parameter for sensitive area prediction [Various]"},
      /* 127 */ {"AT", "Atmospheric tide []"},
      /* 128 */ {"BV", "Budget values []"},
      /* 129 */ {"Z", "Geopotential [m**2 s**-2]"},
      /* 130 */ {"T", "Temperature [K]"},
      /* 131 */ {"U", "U velocity [m s**-1]"},
      /* 132 */ {"V", "V velocity [m s**-1]"},
      /* 133 */ {"Q", "Specific humidity [kg kg**-1]"},
      /* 134 */ {"SP", "Surface pressure [Pa]"},
      /* 135 */ {"W", "Vertical velocity [Pa s**-1]"},
      /* 136 */ {"TCW", "Total column water [kg m**-2]"},
      /* 137 */ {"TCWV", "Total column water vapour [kg m**-2]"},
      /* 138 */ {"VO", "Vorticity (relative) [s**-1]"},
      /* 139 */ {"STL1", "Soil temperature level 1 [K]"},
      /* 140 */ {"SWL1", "Soil wetness level 1 [m of water]"},
      /* 141 */ {"SD", "Snow depth [m of water equivalent]"},
      /* 142 */ {"LSP", "Stratiform precipitation (Large-scale precipitation) [m]"},
      /* 143 */ {"CP", "Convective precipitation [m]"},
      /* 144 */ {"SF", "Snowfall [m of water equivalent]"},
      /* 145 */ {"BLD", "Boundary layer dissipation [W m**-2 s]"},
      /* 146 */ {"SSHF", "Surface sensible heat flux [W m**-2 s]"},
      /* 147 */ {"SLHF", "Surface latent heat flux [W m**-2 s]"},
      /* 148 */ {"CHNK", "Charnock []"},
      /* 149 */ {"SNR", "Surface net radiation [W m**-2 s]"},
      /* 150 */ {"TNR", "Top net radiation []"},
      /* 151 */ {"MSL", "Mean sea level pressure [Pa]"},
      /* 152 */ {"LNSP", "Logarithm of surface pressure []"},
      /* 153 */ {"SWHR", "Short-wave heating rate [K]"},
      /* 154 */ {"LWHR", "Long-wave heating rate [K]"},
      /* 155 */ {"D", "Divergence [s**-1]"},
      /* 156 */ {"GH", "Height [gpm]"},
      /* 157 */ {"R", "Relative humidity [%]"},
      /* 158 */ {"TSP", "Tendency of surface pressure [Pa s**-1]"},
      /* 159 */ {"BLH", "Boundary layer height [m]"},
      /* 160 */ {"SDOR", "Standard deviation of orography []"},
      /* 161 */ {"ISOR", "Anisotropy of sub-gridscale orography []"},
      /* 162 */ {"ANOR", "Angle of sub-gridscale orography [rad]"},
      /* 163 */ {"SLOR", "Slope of sub-gridscale orography []"},
      /* 164 */ {"TCC", "Total cloud cover [(0 - 1)]"},
      /* 165 */ {"10U", "10 metre U wind component [m s**-1]"},
      /* 166 */ {"10V", "10 metre V wind component [m s**-1]"},
      /* 167 */ {"2T", "2 metre temperature [K]"},
      /* 168 */ {"2D", "2 metre dewpoint temperature [K]"},
      /* 169 */ {"SSRD", "Surface solar radiation downwards [W m**-2 s]"},
      /* 170 */ {"STL2", "Soil temperature level 2 [K]"},
      /* 171 */ {"SWL2", "Soil wetness level 2 [m of water]"},
      /* 172 */ {"LSM", "Land-sea mask [(0 - 1)]"},
      /* 173 */ {"SR", "Surface roughness [m]"},
      /* 174 */ {"AL", "Albedo [(0 - 1)]"},
      /* 175 */ {"STRD", "Surface thermal radiation downwards [W m**-2 s]"},
      /* 176 */ {"SSR", "Surface solar radiation [W m**-2 s]"},
      /* 177 */ {"STR", "Surface thermal radiation [W m**-2 s]"},
      /* 178 */ {"TSR", "Top solar radiation [W m**-2 s]"},
      /* 179 */ {"TTR", "Top thermal radiation [W m**-2 s]"},
      /* 180 */ {"EWSS", "East-West surface stress [N m**-2 s]"},
      /* 181 */ {"NSSS", "North-South surface stress [N m**-2 s]"},
      /* 182 */ {"E", "Evaporation [m of water]"},
      /* 183 */ {"STL3", "Soil temperature level 3 [K]"},
      /* 184 */ {"SWL3", "Soil wetness level 3 [m of water]"},
      /* 185 */ {"CCC", "Convective cloud cover [(0 - 1)]"},
      /* 186 */ {"LCC", "Low cloud cover [(0 - 1)]"},
      /* 187 */ {"MCC", "Medium cloud cover [(0 - 1)]"},
      /* 188 */ {"HCC", "High cloud cover [(0 - 1)]"},
      /* 189 */ {"SUND", "Sunshine duration [s]"},
      /* 190 */ {"EWOV", "East-West component of sub-gridscale orographic variance [m**2]"},
      /* 191 */ {"NSOV", "North-South component of sub-gridscale orographic variance [m**2]"},
      /* 192 */ {"NWOV", "North-West/South-East component of sub-gridscale orographic variance [m**2]"},
      /* 193 */ {"NEOV", "North-East/South-West component of sub-gridscale orographic variance [m**2]"},
      /* 194 */ {"BTMP", "Brightness temperature [K]"},
      /* 195 */ {"LGWS", "Latitudinal component of gravity wave stress [N m**-2 s]"},
      /* 196 */ {"MGWS", "Meridional component of gravity wave stress [N m**-2 s]"},
      /* 197 */ {"GWD", "Gravity wave dissipation [W m**-2 s]"},
      /* 198 */ {"SRC", "Skin reservoir content [m of water]"},
      /* 199 */ {"VEG", "Vegetation fraction [(0 - 1)]"},
      /* 200 */ {"VSO", "Variance of sub-gridscale orography [m**2]"},
      /* 201 */ {"MX2T", "Maximum temperature at 2 metres since previous post-processing [K]"},
      /* 202 */ {"MN2T", "Minimum temperature at 2 metres since previous post-processing [K]"},
      /* 203 */ {"O3", "Ozone mass mixing ratio [kg kg**-1]"},
      /* 204 */ {"PAW", "Precipitation analysis weights []"},
      /* 205 */ {"RO", "Runoff [m]"},
      /* 206 */ {"TCO3", "Total column ozone [kg m**-2]"},
      /* 207 */ {"10SI", "10 metre wind speed [m s**-1]"},
      /* 208 */ {"TSRC", "Top net solar radiation, clear sky [W m**-2 s]"},
      /* 209 */ {"TTRC", "Top net thermal radiation, clear sky [W m**-2 s]"},
      /* 210 */ {"SSRC", "Surface net solar radiation, clear sky [W m**-2 s]"},
      /* 211 */ {"STRC", "Surface net thermal radiation, clear sky [W m**-2 s]"},
      /* 212 */ {"TISR", "TOA incident solar radiation [W m**-2 s]"},
      /* 213 */ {"VIMD", "Vertically integrated moisture divergence [kg m**-2]"},
      /* 214 */ {"DHR", "Diabatic heating by radiation [K]"},
      /* 215 */ {"DHVD", "Diabatic heating by vertical diffusion [K]"},
      /* 216 */ {"DHCC", "Diabatic heating by cumulus convection [K]"},
      /* 217 */ {"DHLC", "Diabatic heating large-scale condensation [K]"},
      /* 218 */ {"VDZW", "Vertical diffusion of zonal wind [m s**-1]"},
      /* 219 */ {"VDMW", "Vertical diffusion of meridional wind [m s**-1]"},
      /* 220 */ {"EWGD", "East-West gravity wave drag tendency [m s**-1]"},
      /* 221 */ {"NSGD", "North-South gravity wave drag tendency [m s**-1]"},
      /* 222 */ {"CTZW", "Convective tendency of zonal wind [m s**-1]"},
      /* 223 */ {"CTMW", "Convective tendency of meridional wind [m s**-1]"},
      /* 224 */ {"VDH", "Vertical diffusion of humidity [kg kg**-1]"},
      /* 225 */ {"HTCC", "Humidity tendency by cumulus convection [kg kg**-1]"},
      /* 226 */ {"HTLC", "Humidity tendency by large-scale condensation [kg kg**-1]"},
      /* 227 */ {"CRNH", "Change from removal of negative humidity [kg kg**-1]"},
      /* 228 */ {"TP", "Total precipitation [m]"},
      /* 229 */ {"IEWS", "Instantaneous X surface stress [N m**-2]"},
      /* 230 */ {"INSS", "Instantaneous Y surface stress [N m**-2]"},
      /* 231 */ {"ISHF", "Instantaneous surface heat flux [W m**-2]"},
      /* 232 */ {"IE", "Instantaneous moisture flux [kg m**-2 s**-1]"},
      /* 233 */ {"ASQ", "Apparent surface humidity [kg kg**-1]"},
      /* 234 */ {"LSRH", "Logarithm of surface roughness length for heat []"},
      /* 235 */ {"SKT", "Skin temperature [K]"},
      /* 236 */ {"STL4", "Soil temperature level 4 [K]"},
      /* 237 */ {"SWL4", "Soil wetness level 4 [m]"},
      /* 238 */ {"TSN", "Temperature of snow layer [K]"},
      /* 239 */ {"CSF", "Convective snowfall [m of water equivalent]"},
      /* 240 */ {"LSF", "Large-scale snowfall [m of water equivalent]"},
      /* 241 */ {"ACF", "Accumulated cloud fraction tendency [(-1 to 1)]"},
      /* 242 */ {"ALW", "Accumulated liquid water tendency [(-1 to 1)]"},
      /* 243 */ {"FAL", "Forecast albedo [(0 - 1)]"},
      /* 244 */ {"FSR", "Forecast surface roughness [m]"},
      /* 245 */ {"FLSR", "Forecast logarithm of surface roughness for heat []"},
      /* 246 */ {"CLWC", "Cloud liquid water content [kg kg**-1]"},
      /* 247 */ {"CIWC", "Cloud ice water content [kg kg**-1]"},
      /* 248 */ {"CC", "Cloud cover [(0 - 1)]"},
      /* 249 */ {"AIW", "Accumulated ice water tendency [(-1 to 1)]"},
      /* 250 */ {"ICE", "Ice age [(0 - 1)]"},
      /* 251 */ {"ATTE", "Adiabatic tendency of temperature [K]"},
      /* 252 */ {"ATHE", "Adiabatic tendency of humidity [kg kg**-1]"},
      /* 253 */ {"ATZE", "Adiabatic tendency of zonal wind [m s**-1]"},
      /* 254 */ {"ATMW", "Adiabatic tendency of meridional wind [m s**-1]"},
      /* 255 */ {"var255", "Indicates a missing value []"},
};
