#include <math.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
extern long dadt_counter;
extern double InfusionRate[99];
extern double *par_ptr;
extern double podo;
extern double tlast;

// prj-specific differential eqns
void RxODE_mod_ivsc_2cmtc_full_dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double
	D,
	Ac,
	Vc,
	Ad,
	ka,
	F,
	k12,
	k21,
	Ap,
	keD,
	kon,
	T,
	koff,
	DT,
	ksyn,
	keT,
	keDT;

	Vc = par_ptr[0];
	ka = par_ptr[1];
	F = par_ptr[2];
	k12 = par_ptr[3];
	k21 = par_ptr[4];
	keD = par_ptr[5];
	kon = par_ptr[6];
	koff = par_ptr[7];
	ksyn = par_ptr[8];
	keT = par_ptr[9];
	keDT = par_ptr[10];

	Ad = __zzStateVar__[0];
	Ac = __zzStateVar__[1];
	Ap = __zzStateVar__[2];
	T = __zzStateVar__[3];
	DT = __zzStateVar__[4];

	D = Ac / Vc;
	__DDtStateVar__[0] = InfusionRate[0] + - ka * Ad;
	__DDtStateVar__[1] = InfusionRate[1] + F * ka * Ad - k12 * Ac + k21 * Ap - keD * Ac +( - kon * D * T + koff * DT) * Vc;
	__DDtStateVar__[2] = InfusionRate[2] + k12 * Ac - k21 * Ap;
	__DDtStateVar__[3] = InfusionRate[3] + ksyn - keT * T - kon * D * T + koff * DT;
	__DDtStateVar__[4] = InfusionRate[4] + - keDT * DT + kon * D * T - koff * DT;
    dadt_counter++;
}

// prj-specific derived vars
void RxODE_mod_ivsc_2cmtc_full_calc_lhs(double t, double *__zzStateVar__, double *lhs) {
double
	D,
	Ac,
	Vc,
	Ad,
	ka,
	F,
	k12,
	k21,
	Ap,
	keD,
	kon,
	T,
	koff,
	DT,
	ksyn,
	keT,
	keDT;

	Vc = par_ptr[0];
	ka = par_ptr[1];
	F = par_ptr[2];
	k12 = par_ptr[3];
	k21 = par_ptr[4];
	keD = par_ptr[5];
	kon = par_ptr[6];
	koff = par_ptr[7];
	ksyn = par_ptr[8];
	keT = par_ptr[9];
	keDT = par_ptr[10];

	Ad = __zzStateVar__[0];
	Ac = __zzStateVar__[1];
	Ap = __zzStateVar__[2];
	T = __zzStateVar__[3];
	DT = __zzStateVar__[4];

	D = Ac / Vc;

	lhs[0]=D;
}
