

// ------ GEN SPECTRA BINS ------
const double simpbins_pt_gen[] = {//simple 15 GeV pt bins
  25.,//junk
  40., 55., 70., 85., 
  100., 115., 130., 145., 160., 175., 190., 
  205., 220., 235., 250., 265., 280., 295., 
  310., 325., 340., 355., 370., 385., 
  400., 415., 430., 445., 460., 475., 490., 
  505., 520., 535., 550., 565., 580., 595., 
  610., 625., 640., 655., 670., 685., 
  700., 715., 730., 745., 760., 775., 790., 
  805., 820., 835., 850., 865., 880., 895., 
  910., 925., 940., 955., 970., 985., 
  1000., 1015., 1030., 1045., 1060., 1075., 1090., 
  1105., 1120., 1135., 1150., 1165., 1180., 1195., 
  1210., 1225., 1240., 1255., 1270., 1285., 
  1300., 1315., 1330., 1345., 1360., 1375., 1390., 
  1405., 1420., 1435., 1450., 1465., 1480., 1495., 
  1510., 1525., 1540., 1555., 1570., 1585., 
  1600., 1615., 1630., 1645., 1660., 1675., 1690., 
  1705., 1720., 1735., 1750., 1765., 1780., 1795., 
  1810., 1825., 1840., 1855., 1870., 1885., 
  1900., 1915., 1930., 1945., 1960., 1975., 1990., 
  2005.
//last bin? with content
};
const int n_simpbins_pt_gen = sizeof(simpbins_pt_gen)/sizeof(double)-1;



// ------ RECO SPECTRA BINS ------
const double simpbins_pt_reco[] = {//simple 10 GeV pt bins
  40., //junk
  55., 70., 85., 
  100., 115., 130., 145., 160., 175., 190., 
  205., 220., 235., 250., 265., 280., 295., 
  310., 325., 340., 355., 370., 385., 
  400., 415., 430., 445., 460., 475., 490., 
  505., 520., 535., 550., 565., 580., 595., 
  610., 625., 640., 655., 670., 685., 
  700., 715., 730., 745., 760., 775., 790., 
  805., 820., 835., 850., 865., 880., 895., 
  910., 925., 940., 955., 970., 985., 
  1000., 1015., 1030., 1045., 1060., 1075., 1090., 
  1105., 1120., 1135., 1150., 1165., 1180., 1195., 
  1210., 1225., 1240., 1255., 1270., 1285., 
  1300., 1315., 1330., 1345., 1360., 1375., 1390., 
  1405., 1420., 1435., 1450., 1465., 1480., 1495., 
  1510., 1525., 1540., 1555., 1570., 1585., 
  1600., 1615., 1630., 1645., 1660., 1675., 1690., 
  1705., 1720., 1735., 1750., 1765., 1780., 1795., 
  1810., 1825., 1840., 1855., 1870., 1885., 
  1900., 1915., 1930., 1945., 1960., 1975., 1990., 
  2005.
};
const int n_simpbins_pt_reco = sizeof(simpbins_pt_reco)/sizeof(double)-1;

/* NOTE
For easy bin changing, use regex:
find ([0-9])([0-9])0
replace $1$2(2)

will change **0 to **2

*/



















