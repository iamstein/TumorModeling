#include <math.h>
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
extern long dadt_counter;
extern double InfusionRate[99];
extern double *par_ptr;
extern double podo;
extern double tlast;

// prj-specific differential eqns
void RxODE_mod_ivsc_3cmtct_shed3_dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
double
	D1,
	AmtD1,
	VD1,
	AmtD0,
	ka,
	F,
	k13D,
	k31D,
	VD3,
	D3,
	keD1,
	kon1,
	S1,
	koff1,
	DS1,
	k12D,
	k21D,
	VD2,
	D2,
	keD3,
	kon3,
	S3,
	M3,
	koff3,
	DS3,
	DM3,
	ksynS1,
	k13S,
	k31S,
	VS3,
	VS1,
	keS1,
	ksynS3,
	kshedM3,
	keS3,
	ksynM3,
	keM3,
	k13DS,
	k31DS,
	VDS3,
	VDS1,
	keDS1,
	kshedDM3,
	keDS3,
	keDM3;

	VD1 = par_ptr[0];
	ka = par_ptr[1];
	F = par_ptr[2];
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
	ksynS1 = par_ptr[15];
	k13S = par_ptr[16];
	k31S = par_ptr[17];
	VS3 = par_ptr[18];
	VS1 = par_ptr[19];
	keS1 = par_ptr[20];
	ksynS3 = par_ptr[21];
	kshedM3 = par_ptr[22];
	keS3 = par_ptr[23];
	ksynM3 = par_ptr[24];
	keM3 = par_ptr[25];
	k13DS = par_ptr[26];
	k31DS = par_ptr[27];
	VDS3 = par_ptr[28];
	VDS1 = par_ptr[29];
	keDS1 = par_ptr[30];
	kshedDM3 = par_ptr[31];
	keDS3 = par_ptr[32];
	keDM3 = par_ptr[33];

	AmtD0 = __zzStateVar__[0];
	AmtD1 = __zzStateVar__[1];
	D2 = __zzStateVar__[2];
	D3 = __zzStateVar__[3];
	S1 = __zzStateVar__[4];
	S3 = __zzStateVar__[5];
	M3 = __zzStateVar__[6];
	DS1 = __zzStateVar__[7];
	DS3 = __zzStateVar__[8];
	DM3 = __zzStateVar__[9];

	D1 = AmtD1 / VD1;
	__DDtStateVar__[0] = InfusionRate[0] + - ka * AmtD0;
	__DDtStateVar__[1] = InfusionRate[1] +( F * ka * AmtD0 / VD1 - k13D * D1 + k31D * VD3 / VD1 * D3 - keD1 * D1 - kon1 * D1 * S1 + koff1 * DS1 - k12D * D1 + k21D * VD2 / VD1 * D2) * VD1;
	__DDtStateVar__[2] = InfusionRate[2] + k12D * VD1 / VD2 * D1 - k21D * D2;
	__DDtStateVar__[3] = InfusionRate[3] + k13D * VD1 / VD3 * D1 - k31D * D3 - keD3 * D3 - kon3 * D3 *( S3 + M3) + koff3 *( DS3 + DM3);
	__DDtStateVar__[4] = InfusionRate[4] + ksynS1 - k13S * S1 + k31S * VS3 / VS1 * S3 - keS1 * S1 - kon1 * D1 * S1 + koff1 * DS1;
	__DDtStateVar__[5] = InfusionRate[5] + ksynS3 + kshedM3 * M3 + k13S * VS1 / VS3 * S1 - k31S * S3 - keS3 * S3 - kon3 * D3 * S3 + koff3 * DS3;
	__DDtStateVar__[6] = InfusionRate[6] + ksynM3 - kshedM3 * M3 - keM3 * M3 - kon3 * D3 * M3 + koff3 * DM3;
	__DDtStateVar__[7] = InfusionRate[7] + - k13DS * DS1 + k31DS * VDS3 / VDS1 * DS3 - keDS1 * DS1 + kon1 * D1 * S1 - koff1 * DS1;
	__DDtStateVar__[8] = InfusionRate[8] + kshedDM3 * DM3 + k13DS * VDS1 / VDS3 * DS1 - k31DS * DS3 - keDS3 * DS3 + kon3 * D3 * S3 - koff3 * DS3;
	__DDtStateVar__[9] = InfusionRate[9] + - kshedDM3 * DM3 - keDM3 * DM3 + kon3 * D3 * M3 - koff3 * DM3;
    dadt_counter++;
}

// prj-specific derived vars
void RxODE_mod_ivsc_3cmtct_shed3_calc_lhs(double t, double *__zzStateVar__, double *lhs) {
double
	D1,
	AmtD1,
	VD1,
	AmtD0,
	ka,
	F,
	k13D,
	k31D,
	VD3,
	D3,
	keD1,
	kon1,
	S1,
	koff1,
	DS1,
	k12D,
	k21D,
	VD2,
	D2,
	keD3,
	kon3,
	S3,
	M3,
	koff3,
	DS3,
	DM3,
	ksynS1,
	k13S,
	k31S,
	VS3,
	VS1,
	keS1,
	ksynS3,
	kshedM3,
	keS3,
	ksynM3,
	keM3,
	k13DS,
	k31DS,
	VDS3,
	VDS1,
	keDS1,
	kshedDM3,
	keDS3,
	keDM3;

	VD1 = par_ptr[0];
	ka = par_ptr[1];
	F = par_ptr[2];
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
	ksynS1 = par_ptr[15];
	k13S = par_ptr[16];
	k31S = par_ptr[17];
	VS3 = par_ptr[18];
	VS1 = par_ptr[19];
	keS1 = par_ptr[20];
	ksynS3 = par_ptr[21];
	kshedM3 = par_ptr[22];
	keS3 = par_ptr[23];
	ksynM3 = par_ptr[24];
	keM3 = par_ptr[25];
	k13DS = par_ptr[26];
	k31DS = par_ptr[27];
	VDS3 = par_ptr[28];
	VDS1 = par_ptr[29];
	keDS1 = par_ptr[30];
	kshedDM3 = par_ptr[31];
	keDS3 = par_ptr[32];
	keDM3 = par_ptr[33];

	AmtD0 = __zzStateVar__[0];
	AmtD1 = __zzStateVar__[1];
	D2 = __zzStateVar__[2];
	D3 = __zzStateVar__[3];
	S1 = __zzStateVar__[4];
	S3 = __zzStateVar__[5];
	M3 = __zzStateVar__[6];
	DS1 = __zzStateVar__[7];
	DS3 = __zzStateVar__[8];
	DM3 = __zzStateVar__[9];

	D1 = AmtD1 / VD1;

	lhs[0]=D1;
}
