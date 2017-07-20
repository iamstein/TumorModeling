#include <math.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
extern long dadt_counter;
extern double InfusionRate[99];
extern double *par_ptr;
extern double podo;
extern double tlast;

// prj-specific differential eqns
void RxODE_mod_ivsc_3cmtct_full_dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double
	AmtD0,
	ka,
	D1,
	F,
	VD1,
	k13D,
	k31D,
	VD3,
	D3,
	keD1,
	kon1,
	T1,
	koff1,
	DT1,
	k12D,
	k21D,
	VD2,
	D2,
	keD3,
	kon3,
	T3,
	koff3,
	DT3,
	ksyn1,
	k13T,
	k31T,
	VT3,
	VT1,
	keT1,
	ksyn3,
	keT3,
	k13DT,
	k31DT,
	VDT3,
	VDT1,
	keDT1,
	keDT3;

	ka = par_ptr[0];
	F = par_ptr[1];
	VD1 = par_ptr[2];
	k13D = par_ptr[3];
	k31D = par_ptr[4];
	VD3 = par_ptr[5];
	keD1 = par_ptr[6];
	kon1 = par_ptr[7];
	koff1 = par_ptr[8];
	k12D = par_ptr[9];
	k21D = par_ptr[10];
	VD2 = par_ptr[11];
	keD3 = par_ptr[12];
	kon3 = par_ptr[13];
	koff3 = par_ptr[14];
	ksyn1 = par_ptr[15];
	k13T = par_ptr[16];
	k31T = par_ptr[17];
	VT3 = par_ptr[18];
	VT1 = par_ptr[19];
	keT1 = par_ptr[20];
	ksyn3 = par_ptr[21];
	keT3 = par_ptr[22];
	k13DT = par_ptr[23];
	k31DT = par_ptr[24];
	VDT3 = par_ptr[25];
	VDT1 = par_ptr[26];
	keDT1 = par_ptr[27];
	keDT3 = par_ptr[28];

	AmtD0 = __zzStateVar__[0];
	D1 = __zzStateVar__[1];
	D2 = __zzStateVar__[2];
	D3 = __zzStateVar__[3];
	T1 = __zzStateVar__[4];
	T3 = __zzStateVar__[5];
	DT1 = __zzStateVar__[6];
	DT3 = __zzStateVar__[7];

	__DDtStateVar__[0] = InfusionRate[0] + - ka * AmtD0;
	__DDtStateVar__[1] = InfusionRate[1] + F * ka * AmtD0 / VD1 - k13D * D1 + k31D * VD3 / VD1 * D3 - keD1 * D1 - kon1 * D1 * T1 + koff1 * DT1 - k12D * D1 + k21D * VD2 / VD1 * D2;
	__DDtStateVar__[2] = InfusionRate[2] + k12D * VD1 / VD2 * D1 - k21D * D2;
	__DDtStateVar__[3] = InfusionRate[3] + k13D * VD1 / VD3 * D1 - k31D * D3 - keD3 * D3 - kon3 * D3 * T3 + koff3 * DT3;
	__DDtStateVar__[4] = InfusionRate[4] + ksyn1 - k13T * T1 + k31T * VT3 / VT1 * T3 - keT1 * T1 - kon1 * D1 * T1 + koff1 * DT1;
	__DDtStateVar__[5] = InfusionRate[5] + ksyn3 + k13T * VT1 / VT3 * T1 - k31T * T3 - keT3 * T3 - kon3 * D3 * T3 + koff3 * DT3;
	__DDtStateVar__[6] = InfusionRate[6] + - k13DT * DT1 + k31DT * VDT3 / VDT1 * DT3 - keDT1 * DT1 + kon1 * D1 * T1 - koff1 * DT1;
	__DDtStateVar__[7] = InfusionRate[7] + k13DT * VDT1 / VDT3 * DT1 - k31DT * DT3 - keDT3 * DT3 + kon3 * D3 * T3 - koff3 * DT3;
    dadt_counter++;
}

// prj-specific derived vars
void RxODE_mod_ivsc_3cmtct_full_calc_lhs(double t, double *__zzStateVar__, double *lhs) {
double
	AmtD0,
	ka,
	D1,
	F,
	VD1,
	k13D,
	k31D,
	VD3,
	D3,
	keD1,
	kon1,
	T1,
	koff1,
	DT1,
	k12D,
	k21D,
	VD2,
	D2,
	keD3,
	kon3,
	T3,
	koff3,
	DT3,
	ksyn1,
	k13T,
	k31T,
	VT3,
	VT1,
	keT1,
	ksyn3,
	keT3,
	k13DT,
	k31DT,
	VDT3,
	VDT1,
	keDT1,
	keDT3;

	ka = par_ptr[0];
	F = par_ptr[1];
	VD1 = par_ptr[2];
	k13D = par_ptr[3];
	k31D = par_ptr[4];
	VD3 = par_ptr[5];
	keD1 = par_ptr[6];
	kon1 = par_ptr[7];
	koff1 = par_ptr[8];
	k12D = par_ptr[9];
	k21D = par_ptr[10];
	VD2 = par_ptr[11];
	keD3 = par_ptr[12];
	kon3 = par_ptr[13];
	koff3 = par_ptr[14];
	ksyn1 = par_ptr[15];
	k13T = par_ptr[16];
	k31T = par_ptr[17];
	VT3 = par_ptr[18];
	VT1 = par_ptr[19];
	keT1 = par_ptr[20];
	ksyn3 = par_ptr[21];
	keT3 = par_ptr[22];
	k13DT = par_ptr[23];
	k31DT = par_ptr[24];
	VDT3 = par_ptr[25];
	VDT1 = par_ptr[26];
	keDT1 = par_ptr[27];
	keDT3 = par_ptr[28];

	AmtD0 = __zzStateVar__[0];
	D1 = __zzStateVar__[1];
	D2 = __zzStateVar__[2];
	D3 = __zzStateVar__[3];
	T1 = __zzStateVar__[4];
	T3 = __zzStateVar__[5];
	DT1 = __zzStateVar__[6];
	DT3 = __zzStateVar__[7];


}
