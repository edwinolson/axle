#ifndef MATHUTIL_H
#define MATHUTIL_H

#ifndef M_TWOPI
#define M_TWOPI    6.2831853071795862319959 /* 2*pi */
#endif

#ifndef M_PI
#define M_PI    3.141592653589793238462643383279502884196
#endif

/** Map vin to [0, 2*PI) **/
static inline double mod2pi_positive(double vin)
{
    return vin - M_TWOPI * floor(vin / M_TWOPI);
}

/** Map vin to [-PI, PI) **/
static inline double mod2pi(double vin)
{
    return mod2pi_positive(vin + M_PI) - M_PI;
}

#endif
