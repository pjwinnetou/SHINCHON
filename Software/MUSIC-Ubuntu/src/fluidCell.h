#ifndef fluidCell_H_
#define fluidCell_H_

struct fluidCell {
   double ed, sd, temperature, pressure;
   double vx, vy, vz;
   double pi[4][4];
   double bulkPi;
};

struct fluidCell_ideal {
    float sd, ed, pressure, temperature;
    float ux, uy, ueta;
};

#endif  // fluidCell_H_
