#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -25.46*x up 19.18*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {< -3.03,  -5.34,  -3.05>, < -2.71,  -4.55,  -6.25>, Rcell pigment {Black}}
cylinder {< -2.61,   5.25,  -0.37>, < -2.28,   6.05,  -3.57>, Rcell pigment {Black}}
cylinder {<  1.75,   4.98,   0.01>, <  2.07,   5.78,  -3.19>, Rcell pigment {Black}}
cylinder {<  1.32,  -5.61,  -2.68>, <  1.65,  -4.81,  -5.88>, Rcell pigment {Black}}
cylinder {< -3.03,  -5.34,  -3.05>, < -2.61,   5.25,  -0.37>, Rcell pigment {Black}}
cylinder {< -2.71,  -4.55,  -6.25>, < -2.28,   6.05,  -3.57>, Rcell pigment {Black}}
cylinder {<  1.65,  -4.81,  -5.88>, <  2.07,   5.78,  -3.19>, Rcell pigment {Black}}
cylinder {<  1.32,  -5.61,  -2.68>, <  1.75,   4.98,   0.01>, Rcell pigment {Black}}
cylinder {< -3.03,  -5.34,  -3.05>, <  1.32,  -5.61,  -2.68>, Rcell pigment {Black}}
cylinder {< -2.71,  -4.55,  -6.25>, <  1.65,  -4.81,  -5.88>, Rcell pigment {Black}}
cylinder {< -2.28,   6.05,  -3.57>, <  2.07,   5.78,  -3.19>, Rcell pigment {Black}}
cylinder {< -2.61,   5.25,  -0.37>, <  1.75,   4.98,   0.01>, Rcell pigment {Black}}
atom(< -2.48,  -3.93,  -4.36>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #0 
atom(<  1.09,  -6.22,  -4.57>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #1 
atom(< -1.16,  -4.42,  -2.63>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #2 
atom(< -0.54,  -6.53,  -3.10>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #3 
atom(< -2.43,   0.96,  -1.42>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #4 
atom(<  1.14,  -1.32,  -1.63>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #5 
atom(< -0.79,   1.27,  -2.89>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #6 
atom(< -0.17,  -0.83,  -3.35>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #7 
atom(<  1.88,  -4.20,  -3.99>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #8 
atom(<  5.44,  -6.49,  -4.20>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #9 
atom(<  3.19,  -4.69,  -2.26>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #10 
atom(<  3.81,  -6.80,  -2.72>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #11 
atom(<  1.93,   0.70,  -1.04>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #12 
atom(<  5.49,  -1.59,  -1.25>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #13 
atom(<  3.56,   1.00,  -2.52>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #14 
atom(<  4.18,  -1.10,  -2.98>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #15 
atom(< -2.15,  -3.14,  -7.56>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #16 
atom(<  1.42,  -5.43,  -7.77>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #17 
atom(< -0.84,  -3.63,  -5.83>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #18 
atom(< -0.22,  -5.73,  -6.30>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #19 
atom(< -2.10,   1.76,  -4.62>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #20 
atom(<  1.47,  -0.53,  -4.83>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #21 
atom(< -0.47,   2.07,  -6.09>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #22 
atom(<  0.15,  -0.04,  -6.55>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #23 
atom(<  2.20,  -3.41,  -7.19>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #24 
atom(<  5.77,  -5.69,  -7.40>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #25 
atom(<  3.51,  -3.90,  -5.46>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #26 
atom(<  4.13,  -6.00,  -5.92>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #27 
atom(<  2.25,   1.49,  -4.24>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #28 
atom(<  5.82,  -0.80,  -4.45>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #29 
atom(<  3.89,   1.80,  -5.72>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #30 
atom(<  4.51,  -0.31,  -6.18>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #31 
atom(< -1.83,  -2.34, -10.76>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #32 
atom(<  1.74,  -4.63, -10.97>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #33 
atom(< -0.52,  -2.83,  -9.03>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #34 
atom(<  0.10,  -4.94,  -9.50>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #35 
atom(< -1.78,   2.56,  -7.82>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #36 
atom(<  1.79,   0.27,  -8.03>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #37 
atom(< -0.14,   2.87,  -9.29>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #38 
atom(<  0.48,   0.76,  -9.75>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #39 
atom(<  2.52,  -2.61, -10.39>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #40 
atom(<  6.09,  -4.90, -10.60>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #41 
atom(<  3.84,  -3.10,  -8.66>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #42 
atom(<  4.46,  -5.20,  -9.12>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #43 
atom(<  2.57,   2.29,  -7.44>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #44 
atom(<  6.14,   0.00,  -7.65>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #45 
atom(<  4.21,   2.60,  -8.92>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #46 
atom(<  4.83,   0.49,  -9.38>, 0.95, rgb <1.00, 0.50, 0.00>, 0.0, ase2) // #47 
