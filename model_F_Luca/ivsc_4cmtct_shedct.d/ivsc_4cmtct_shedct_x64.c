#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if (_jac_counter_val() == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if (_dadt_counter_val() == 0) Rprintf
#define LHS_Rprintf Rprintf
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t);
typedef double (*RxODE_transit3)(double t, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_transit4)(double t, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val);
typedef long (*RxODE_cnt) ();
typedef void (*RxODE_inc) ();
typedef double (*RxODE_val) ();
typedef SEXP (*RxODE_ode_solver) (SEXP sexp_theta, SEXP sexp_inits, SEXP sexp_lhs, SEXP sexp_time, SEXP sexp_evid,SEXP sexp_dose, SEXP sexp_pcov, SEXP sexp_cov, SEXP sexp_locf, SEXP sexp_atol, SEXP sexp_rtol, SEXP sexp_hmin, SEXP sexp_hmax, SEXP sexp_h0, SEXP sexp_mxordn, SEXP sexp_mxords, SEXP sexp_mx,SEXP sexp_stiff, SEXP sexp_transit_abs, SEXP sexp_object, SEXP sexp_extra_args);
typedef void (*RxODE_assign_fn_pointers)(void (*fun_dydt)(unsigned int, double, double *, double *),void (*fun_calc_lhs)(double, double *, double *),void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),void (*fun_update_inis)(SEXP _ini_sexp),int fun_jt,int fun_mf, int fun_debug);

typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef void (*RxODE_ode_solver_0_6_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc,double hmin, double hmax,double h0,int mxordn,int mxords,int mxstep);
// Give par pointers
RxODE_vec _par_ptr, _InfusionRate;
RxODE_update_par_ptr _update_par_ptr;
RxODE_cnt _dadt_counter_val, _jac_counter_val;
RxODE_inc _dadt_counter_inc, _jac_counter_inc;
RxODE_val podo, tlast;
RxODE_transit4 _transit4;
RxODE_transit3 _transit3;
RxODE_fn _safe_log, _safe_zero, factorial, _as_zero;
RxODE_assign_fn_pointers _assign_fn_pointers;
RxODE_ode_solver_old_c _old_c;
RxODE_ode_solver_0_6_c _c_0_6;

extern void __ODE_SOLVER_PTR__();



// prj-specific differential eqns
void ivsc_4cmtct_shedct_x64_dydt(unsigned int _neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
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
	M1,
	DM1,
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
	kshedM1,
	k13S,
	k31S,
	VS3,
	VS1,
	keS1,
	ksynS3,
	kshedM3,
	keS3,
	ksynM3,
	k31M,
	k13M,
	keM3,
	k13DS,
	k31DS,
	VDS3,
	VDS1,
	keDS1,
	kshedDM3,
	keDS3,
	keDM3,
	k31DM,
	k13DM,
	keDM1,
	kshedDM1,
	ksynM1,
	keM1;

	_update_par_ptr(t);
	VD1 = _par_ptr(0);
	ka = _par_ptr(1);
	F = _par_ptr(2);
	k13D = _par_ptr(3);
	k31D = _par_ptr(4);
	VD3 = _par_ptr(5);
	keD1 = _par_ptr(6);
	kon1 = _par_ptr(7);
	koff1 = _par_ptr(8);
	k12D = _par_ptr(9);
	k21D = _par_ptr(10);
	VD2 = _par_ptr(11);
	keD3 = _par_ptr(12);
	kon3 = _par_ptr(13);
	koff3 = _par_ptr(14);
	ksynS1 = _par_ptr(15);
	kshedM1 = _par_ptr(16);
	k13S = _par_ptr(17);
	k31S = _par_ptr(18);
	VS3 = _par_ptr(19);
	VS1 = _par_ptr(20);
	keS1 = _par_ptr(21);
	ksynS3 = _par_ptr(22);
	kshedM3 = _par_ptr(23);
	keS3 = _par_ptr(24);
	ksynM3 = _par_ptr(25);
	k31M = _par_ptr(26);
	k13M = _par_ptr(27);
	keM3 = _par_ptr(28);
	k13DS = _par_ptr(29);
	k31DS = _par_ptr(30);
	VDS3 = _par_ptr(31);
	VDS1 = _par_ptr(32);
	keDS1 = _par_ptr(33);
	kshedDM3 = _par_ptr(34);
	keDS3 = _par_ptr(35);
	keDM3 = _par_ptr(36);
	k31DM = _par_ptr(37);
	k13DM = _par_ptr(38);
	keDM1 = _par_ptr(39);
	kshedDM1 = _par_ptr(40);
	ksynM1 = _par_ptr(41);
	keM1 = _par_ptr(42);

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
	DM1 = __zzStateVar__[10];
	M1 = __zzStateVar__[11];

	D1=AmtD1/_safe_zero(VD1);
	__DDtStateVar__[0] = _InfusionRate(0) + -ka*AmtD0;
	__DDtStateVar__[1] = _InfusionRate(1) + (F*ka*AmtD0/_safe_zero(VD1)-k13D*D1+k31D*VD3/_safe_zero(VD1)*D3-keD1*D1-kon1*D1*S1+koff1*DS1-kon1*D1*M1+koff1*DM1-k12D*D1+k21D*VD2/_safe_zero(VD1)*D2)*VD1;
	__DDtStateVar__[2] = _InfusionRate(2) + k12D*VD1/_safe_zero(VD2)*D1-k21D*D2;
	__DDtStateVar__[3] = _InfusionRate(3) + k13D*VD1/_safe_zero(VD3)*D1-k31D*D3-keD3*D3-kon3*D3*(S3+M3)+koff3*(DS3+DM3);
	__DDtStateVar__[4] = _InfusionRate(4) + ksynS1+kshedM1*M1-k13S*S1+k31S*VS3/_safe_zero(VS1)*S3-keS1*S1-kon1*D1*S1+koff1*DS1;
	__DDtStateVar__[5] = _InfusionRate(5) + ksynS3+kshedM3*M3+k13S*VS1/_safe_zero(VS3)*S1-k31S*S3-keS3*S3-kon3*D3*S3+koff3*DS3;
	__DDtStateVar__[6] = _InfusionRate(6) + ksynM3-kshedM3*M3-k31M*M3+k13M*VD1/_safe_zero(VD3)*M1-keM3*M3-kon3*D3*M3+koff3*DM3;
	__DDtStateVar__[7] = _InfusionRate(7) + +kshedM1*DM1-k13DS*DS1+k31DS*VDS3/_safe_zero(VDS1)*DS3-keDS1*DS1+kon1*D1*S1-koff1*DS1;
	__DDtStateVar__[8] = _InfusionRate(8) + kshedDM3*DM3+k13DS*VDS1/_safe_zero(VDS3)*DS1-k31DS*DS3-keDS3*DS3+kon3*D3*S3-koff3*DS3;
	__DDtStateVar__[9] = _InfusionRate(9) + -kshedDM3*DM3-keDM3*DM3+kon3*D3*M3-koff3*DM3-k31DM*DM3+k13DM*(VD1/_safe_zero(VD3))*DM1;
	__DDtStateVar__[10] = _InfusionRate(10) + -keDM1*DM1-kshedDM1*DM1+kon1*D1*M1-koff1*DM1-k13DM*DM1+k31DM*(VD3/_safe_zero(VD1))*DM3;
	__DDtStateVar__[11] = _InfusionRate(11) + ksynM1-kshedM1*M1-keM1*M1+k31M*VD3/_safe_zero(VD1)*M3-k13M*M1-kon1*D1*M1+koff1*DM1;
    _dadt_counter_inc();
}

// Jacobian derived vars
void ivsc_4cmtct_shedct_x64_calc_jac(unsigned int _neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  _jac_counter_inc();
}
// Functional based initial conditions.
void ivsc_4cmtct_shedct_x64_inis(SEXP _ini_sexp){
	double *__zzStateVar__ = REAL(_ini_sexp);
	double t=0;
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
	M1,
	DM1,
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
	kshedM1,
	k13S,
	k31S,
	VS3,
	VS1,
	keS1,
	ksynS3,
	kshedM3,
	keS3,
	ksynM3,
	k31M,
	k13M,
	keM3,
	k13DS,
	k31DS,
	VDS3,
	VDS1,
	keDS1,
	kshedDM3,
	keDS3,
	keDM3,
	k31DM,
	k13DM,
	keDM1,
	kshedDM1,
	ksynM1,
	keM1;

	_update_par_ptr(0.0);
	VD1 = _par_ptr(0);
	ka = _par_ptr(1);
	F = _par_ptr(2);
	k13D = _par_ptr(3);
	k31D = _par_ptr(4);
	VD3 = _par_ptr(5);
	keD1 = _par_ptr(6);
	kon1 = _par_ptr(7);
	koff1 = _par_ptr(8);
	k12D = _par_ptr(9);
	k21D = _par_ptr(10);
	VD2 = _par_ptr(11);
	keD3 = _par_ptr(12);
	kon3 = _par_ptr(13);
	koff3 = _par_ptr(14);
	ksynS1 = _par_ptr(15);
	kshedM1 = _par_ptr(16);
	k13S = _par_ptr(17);
	k31S = _par_ptr(18);
	VS3 = _par_ptr(19);
	VS1 = _par_ptr(20);
	keS1 = _par_ptr(21);
	ksynS3 = _par_ptr(22);
	kshedM3 = _par_ptr(23);
	keS3 = _par_ptr(24);
	ksynM3 = _par_ptr(25);
	k31M = _par_ptr(26);
	k13M = _par_ptr(27);
	keM3 = _par_ptr(28);
	k13DS = _par_ptr(29);
	k31DS = _par_ptr(30);
	VDS3 = _par_ptr(31);
	VDS1 = _par_ptr(32);
	keDS1 = _par_ptr(33);
	kshedDM3 = _par_ptr(34);
	keDS3 = _par_ptr(35);
	keDM3 = _par_ptr(36);
	k31DM = _par_ptr(37);
	k13DM = _par_ptr(38);
	keDM1 = _par_ptr(39);
	kshedDM1 = _par_ptr(40);
	ksynM1 = _par_ptr(41);
	keM1 = _par_ptr(42);

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
	DM1 = __zzStateVar__[10];
	M1 = __zzStateVar__[11];

	D1=AmtD1/_safe_zero(VD1);
  __zzStateVar__[0]=AmtD0;
  __zzStateVar__[1]=AmtD1;
  __zzStateVar__[2]=D2;
  __zzStateVar__[3]=D3;
  __zzStateVar__[4]=S1;
  __zzStateVar__[5]=S3;
  __zzStateVar__[6]=M3;
  __zzStateVar__[7]=DS1;
  __zzStateVar__[8]=DS3;
  __zzStateVar__[9]=DM3;
  __zzStateVar__[10]=DM1;
  __zzStateVar__[11]=M1;
}
// prj-specific derived vars
void ivsc_4cmtct_shedct_x64_calc_lhs(double t, double *__zzStateVar__, double *_lhs) {
double 
		__DDtStateVar_0__,
	__DDtStateVar_1__,
	__DDtStateVar_2__,
	__DDtStateVar_3__,
	__DDtStateVar_4__,
	__DDtStateVar_5__,
	__DDtStateVar_6__,
	__DDtStateVar_7__,
	__DDtStateVar_8__,
	__DDtStateVar_9__,
	__DDtStateVar_10__,
	__DDtStateVar_11__,
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
	M1,
	DM1,
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
	kshedM1,
	k13S,
	k31S,
	VS3,
	VS1,
	keS1,
	ksynS3,
	kshedM3,
	keS3,
	ksynM3,
	k31M,
	k13M,
	keM3,
	k13DS,
	k31DS,
	VDS3,
	VDS1,
	keDS1,
	kshedDM3,
	keDS3,
	keDM3,
	k31DM,
	k13DM,
	keDM1,
	kshedDM1,
	ksynM1,
	keM1;

	_update_par_ptr(t);
	VD1 = _par_ptr(0);
	ka = _par_ptr(1);
	F = _par_ptr(2);
	k13D = _par_ptr(3);
	k31D = _par_ptr(4);
	VD3 = _par_ptr(5);
	keD1 = _par_ptr(6);
	kon1 = _par_ptr(7);
	koff1 = _par_ptr(8);
	k12D = _par_ptr(9);
	k21D = _par_ptr(10);
	VD2 = _par_ptr(11);
	keD3 = _par_ptr(12);
	kon3 = _par_ptr(13);
	koff3 = _par_ptr(14);
	ksynS1 = _par_ptr(15);
	kshedM1 = _par_ptr(16);
	k13S = _par_ptr(17);
	k31S = _par_ptr(18);
	VS3 = _par_ptr(19);
	VS1 = _par_ptr(20);
	keS1 = _par_ptr(21);
	ksynS3 = _par_ptr(22);
	kshedM3 = _par_ptr(23);
	keS3 = _par_ptr(24);
	ksynM3 = _par_ptr(25);
	k31M = _par_ptr(26);
	k13M = _par_ptr(27);
	keM3 = _par_ptr(28);
	k13DS = _par_ptr(29);
	k31DS = _par_ptr(30);
	VDS3 = _par_ptr(31);
	VDS1 = _par_ptr(32);
	keDS1 = _par_ptr(33);
	kshedDM3 = _par_ptr(34);
	keDS3 = _par_ptr(35);
	keDM3 = _par_ptr(36);
	k31DM = _par_ptr(37);
	k13DM = _par_ptr(38);
	keDM1 = _par_ptr(39);
	kshedDM1 = _par_ptr(40);
	ksynM1 = _par_ptr(41);
	keM1 = _par_ptr(42);

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
	DM1 = __zzStateVar__[10];
	M1 = __zzStateVar__[11];

	D1=AmtD1/_safe_zero(VD1);
	__DDtStateVar_0__ = _InfusionRate(0) + -ka*AmtD0;
	__DDtStateVar_1__ = _InfusionRate(1) + (F*ka*AmtD0/_safe_zero(VD1)-k13D*D1+k31D*VD3/_safe_zero(VD1)*D3-keD1*D1-kon1*D1*S1+koff1*DS1-kon1*D1*M1+koff1*DM1-k12D*D1+k21D*VD2/_safe_zero(VD1)*D2)*VD1;
	__DDtStateVar_2__ = _InfusionRate(2) + k12D*VD1/_safe_zero(VD2)*D1-k21D*D2;
	__DDtStateVar_3__ = _InfusionRate(3) + k13D*VD1/_safe_zero(VD3)*D1-k31D*D3-keD3*D3-kon3*D3*(S3+M3)+koff3*(DS3+DM3);
	__DDtStateVar_4__ = _InfusionRate(4) + ksynS1+kshedM1*M1-k13S*S1+k31S*VS3/_safe_zero(VS1)*S3-keS1*S1-kon1*D1*S1+koff1*DS1;
	__DDtStateVar_5__ = _InfusionRate(5) + ksynS3+kshedM3*M3+k13S*VS1/_safe_zero(VS3)*S1-k31S*S3-keS3*S3-kon3*D3*S3+koff3*DS3;
	__DDtStateVar_6__ = _InfusionRate(6) + ksynM3-kshedM3*M3-k31M*M3+k13M*VD1/_safe_zero(VD3)*M1-keM3*M3-kon3*D3*M3+koff3*DM3;
	__DDtStateVar_7__ = _InfusionRate(7) + +kshedM1*DM1-k13DS*DS1+k31DS*VDS3/_safe_zero(VDS1)*DS3-keDS1*DS1+kon1*D1*S1-koff1*DS1;
	__DDtStateVar_8__ = _InfusionRate(8) + kshedDM3*DM3+k13DS*VDS1/_safe_zero(VDS3)*DS1-k31DS*DS3-keDS3*DS3+kon3*D3*S3-koff3*DS3;
	__DDtStateVar_9__ = _InfusionRate(9) + -kshedDM3*DM3-keDM3*DM3+kon3*D3*M3-koff3*DM3-k31DM*DM3+k13DM*(VD1/_safe_zero(VD3))*DM1;
	__DDtStateVar_10__ = _InfusionRate(10) + -keDM1*DM1-kshedDM1*DM1+kon1*D1*M1-koff1*DM1-k13DM*DM1+k31DM*(VD3/_safe_zero(VD1))*DM3;
	__DDtStateVar_11__ = _InfusionRate(11) + ksynM1-kshedM1*M1-keM1*M1+k31M*VD3/_safe_zero(VD1)*M3-k13M*M1-kon1*D1*M1+koff1*DM1;

	_lhs[0]=D1;
}
extern SEXP ivsc_4cmtct_shedct_x64_model_vars(){
	SEXP lst    = PROTECT(allocVector(VECSXP, 11));
	SEXP names  = PROTECT(allocVector(STRSXP, 11));
	SEXP params = PROTECT(allocVector(STRSXP, 43));
	SEXP lhs    = PROTECT(allocVector(STRSXP, 1));
	SEXP state  = PROTECT(allocVector(STRSXP, 12));
	SEXP sens   = PROTECT(allocVector(STRSXP, 0));
	SEXP fn_ini = PROTECT(allocVector(STRSXP, 0));
	SEXP dfdy   = PROTECT(allocVector(STRSXP, 0));
	SEXP tran   = PROTECT(allocVector(STRSXP, 12));
	SEXP trann  = PROTECT(allocVector(STRSXP, 12));
	SEXP mmd5   = PROTECT(allocVector(STRSXP, 2));
	SEXP mmd5n  = PROTECT(allocVector(STRSXP, 2));
	SEXP model  = PROTECT(allocVector(STRSXP, 4));
	SEXP modeln = PROTECT(allocVector(STRSXP, 4));
	SET_STRING_ELT(lhs,0,mkChar("D1"));
	SET_STRING_ELT(params,0,mkChar("VD1"));
	SET_STRING_ELT(params,1,mkChar("ka"));
	SET_STRING_ELT(params,2,mkChar("F"));
	SET_STRING_ELT(params,3,mkChar("k13D"));
	SET_STRING_ELT(params,4,mkChar("k31D"));
	SET_STRING_ELT(params,5,mkChar("VD3"));
	SET_STRING_ELT(params,6,mkChar("keD1"));
	SET_STRING_ELT(params,7,mkChar("kon1"));
	SET_STRING_ELT(params,8,mkChar("koff1"));
	SET_STRING_ELT(params,9,mkChar("k12D"));
	SET_STRING_ELT(params,10,mkChar("k21D"));
	SET_STRING_ELT(params,11,mkChar("VD2"));
	SET_STRING_ELT(params,12,mkChar("keD3"));
	SET_STRING_ELT(params,13,mkChar("kon3"));
	SET_STRING_ELT(params,14,mkChar("koff3"));
	SET_STRING_ELT(params,15,mkChar("ksynS1"));
	SET_STRING_ELT(params,16,mkChar("kshedM1"));
	SET_STRING_ELT(params,17,mkChar("k13S"));
	SET_STRING_ELT(params,18,mkChar("k31S"));
	SET_STRING_ELT(params,19,mkChar("VS3"));
	SET_STRING_ELT(params,20,mkChar("VS1"));
	SET_STRING_ELT(params,21,mkChar("keS1"));
	SET_STRING_ELT(params,22,mkChar("ksynS3"));
	SET_STRING_ELT(params,23,mkChar("kshedM3"));
	SET_STRING_ELT(params,24,mkChar("keS3"));
	SET_STRING_ELT(params,25,mkChar("ksynM3"));
	SET_STRING_ELT(params,26,mkChar("k31M"));
	SET_STRING_ELT(params,27,mkChar("k13M"));
	SET_STRING_ELT(params,28,mkChar("keM3"));
	SET_STRING_ELT(params,29,mkChar("k13DS"));
	SET_STRING_ELT(params,30,mkChar("k31DS"));
	SET_STRING_ELT(params,31,mkChar("VDS3"));
	SET_STRING_ELT(params,32,mkChar("VDS1"));
	SET_STRING_ELT(params,33,mkChar("keDS1"));
	SET_STRING_ELT(params,34,mkChar("kshedDM3"));
	SET_STRING_ELT(params,35,mkChar("keDS3"));
	SET_STRING_ELT(params,36,mkChar("keDM3"));
	SET_STRING_ELT(params,37,mkChar("k31DM"));
	SET_STRING_ELT(params,38,mkChar("k13DM"));
	SET_STRING_ELT(params,39,mkChar("keDM1"));
	SET_STRING_ELT(params,40,mkChar("kshedDM1"));
	SET_STRING_ELT(params,41,mkChar("ksynM1"));
	SET_STRING_ELT(params,42,mkChar("keM1"));
	SET_STRING_ELT(state,0,mkChar("AmtD0"));
	SET_STRING_ELT(state,1,mkChar("AmtD1"));
	SET_STRING_ELT(state,2,mkChar("D2"));
	SET_STRING_ELT(state,3,mkChar("D3"));
	SET_STRING_ELT(state,4,mkChar("S1"));
	SET_STRING_ELT(state,5,mkChar("S3"));
	SET_STRING_ELT(state,6,mkChar("M3"));
	SET_STRING_ELT(state,7,mkChar("DS1"));
	SET_STRING_ELT(state,8,mkChar("DS3"));
	SET_STRING_ELT(state,9,mkChar("DM3"));
	SET_STRING_ELT(state,10,mkChar("DM1"));
	SET_STRING_ELT(state,11,mkChar("M1"));
	SET_STRING_ELT(modeln,0,mkChar("model"));
	SET_STRING_ELT(model,0,mkChar("\n     D1           = AmtD1/VD1;\n     d/dt(AmtD0)  =  -ka *AmtD0;\n     d/dt(AmtD1)  =(F*ka *AmtD0/VD1      - k13D *D1 +    k31D *VD3/VD1*D3  - keD1 *D1  - kon1*D1*S1 + koff1*DS1  - kon1*D1*M1 + koff1*DM1 - k12D*D1 + k21D*VD2/VD1*D2)*VD1;\n     d/dt(D2)     =                                                                                               k12D*VD1/VD2*D1 - k21D*D2;\n     d/dt(D3)     =                        k13D *VD1/VD3*D1     - k31D*D3  - keD3 *D3  - kon3*D3*(S3+M3) + koff3*(DS3+DM3);\n     d/dt(S1)     = ksynS1+kshedM1*M1               - k13S *S1 +     k31S*VS3/VS1*S3  - keS1 *S1  - kon1*D1*S1      + koff1*DS1;\n     d/dt(S3)     = ksynS3 +kshedM3*M3   + k13S *VS1/VS3*S1     - k31S*S3  - keS3 *S3  - kon3*D3*S3      + koff3*DS3;\n     d/dt(M3)     = ksynM3 -kshedM3*M3  -k31M*M3+k13M*VD1/VD3*M1             - keM3 *M3  - kon3*D3*M3      + koff3*DM3;\n     d/dt(DS1)    =       +kshedM1*DM1               - k13DS*DS1 + k31DS*VDS3/VDS1*DS3 - keDS1*DS1 + kon1*D1*S1      - koff1*DS1;\n     d/dt(DS3)    =         kshedDM3*DM3 + k13DS*VDS1/VDS3*DS1 - k31DS*DS3 - keDS3*DS3 + kon3*D3*S3      - koff3*DS3;\n     d/dt(DM3)    = -kshedDM3*DM3  - keDM3*DM3 + kon3*D3*M3 - koff3*DM3-k31DM*DM3+k13DM*(VD1/VD3)*DM1;\n     d/dt(DM1)    = -keDM1*DM1 -kshedDM1*DM1 +kon1*D1*M1 -koff1*DM1-k13DM*DM1+k31DM*(VD3/VD1)*DM3;\n     d/dt(M1)     = ksynM1 -kshedM1*M1 -keM1*M1 +k31M*VD3/VD1*M3 -k13M*M1 -kon1*D1*M1 +koff1*DM1;\n  \n"));
	SET_STRING_ELT(modeln,1,mkChar("normModel"));
	SET_STRING_ELT(model,1,mkChar("D1=AmtD1/VD1;\nd/dt(AmtD0)=-ka*AmtD0;\nd/dt(AmtD1)=(F*ka*AmtD0/VD1-k13D*D1+k31D*VD3/VD1*D3-keD1*D1-kon1*D1*S1+koff1*DS1-kon1*D1*M1+koff1*DM1-k12D*D1+k21D*VD2/VD1*D2)*VD1;\nd/dt(D2)=k12D*VD1/VD2*D1-k21D*D2;\nd/dt(D3)=k13D*VD1/VD3*D1-k31D*D3-keD3*D3-kon3*D3*(S3+M3)+koff3*(DS3+DM3);\nd/dt(S1)=ksynS1+kshedM1*M1-k13S*S1+k31S*VS3/VS1*S3-keS1*S1-kon1*D1*S1+koff1*DS1;\nd/dt(S3)=ksynS3+kshedM3*M3+k13S*VS1/VS3*S1-k31S*S3-keS3*S3-kon3*D3*S3+koff3*DS3;\nd/dt(M3)=ksynM3-kshedM3*M3-k31M*M3+k13M*VD1/VD3*M1-keM3*M3-kon3*D3*M3+koff3*DM3;\nd/dt(DS1)=+kshedM1*DM1-k13DS*DS1+k31DS*VDS3/VDS1*DS3-keDS1*DS1+kon1*D1*S1-koff1*DS1;\nd/dt(DS3)=kshedDM3*DM3+k13DS*VDS1/VDS3*DS1-k31DS*DS3-keDS3*DS3+kon3*D3*S3-koff3*DS3;\nd/dt(DM3)=-kshedDM3*DM3-keDM3*DM3+kon3*D3*M3-koff3*DM3-k31DM*DM3+k13DM*(VD1/VD3)*DM1;\nd/dt(DM1)=-keDM1*DM1-kshedDM1*DM1+kon1*D1*M1-koff1*DM1-k13DM*DM1+k31DM*(VD3/VD1)*DM3;\nd/dt(M1)=ksynM1-kshedM1*M1-keM1*M1+k31M*VD3/VD1*M3-k13M*M1-kon1*D1*M1+koff1*DM1;\n"));
	SET_STRING_ELT(modeln,2,mkChar("parseModel"));
	SET_STRING_ELT(model,2,mkChar("D1=AmtD1/_safe_zero(VD1);\n__DDtStateVar__[0] = _InfusionRate(0) + -ka*AmtD0;\n__DDtStateVar__[1] = _InfusionRate(1) + (F*ka*AmtD0/_safe_zero(VD1)-k13D*D1+k31D*VD3/_safe_zero(VD1)*D3-keD1*D1-kon1*D1*S1+koff1*DS1-kon1*D1*M1+koff1*DM1-k12D*D1+k21D*VD2/_safe_zero(VD1)*D2)*VD1;\n__DDtStateVar__[2] = _InfusionRate(2) + k12D*VD1/_safe_zero(VD2)*D1-k21D*D2;\n__DDtStateVar__[3] = _InfusionRate(3) + k13D*VD1/_safe_zero(VD3)*D1-k31D*D3-keD3*D3-kon3*D3*(S3+M3)+koff3*(DS3+DM3);\n__DDtStateVar__[4] = _InfusionRate(4) + ksynS1+kshedM1*M1-k13S*S1+k31S*VS3/_safe_zero(VS1)*S3-keS1*S1-kon1*D1*S1+koff1*DS1;\n__DDtStateVar__[5] = _InfusionRate(5) + ksynS3+kshedM3*M3+k13S*VS1/_safe_zero(VS3)*S1-k31S*S3-keS3*S3-kon3*D3*S3+koff3*DS3;\n__DDtStateVar__[6] = _InfusionRate(6) + ksynM3-kshedM3*M3-k31M*M3+k13M*VD1/_safe_zero(VD3)*M1-keM3*M3-kon3*D3*M3+koff3*DM3;\n__DDtStateVar__[7] = _InfusionRate(7) + +kshedM1*DM1-k13DS*DS1+k31DS*VDS3/_safe_zero(VDS1)*DS3-keDS1*DS1+kon1*D1*S1-koff1*DS1;\n__DDtStateVar__[8] = _InfusionRate(8) + kshedDM3*DM3+k13DS*VDS1/_safe_zero(VDS3)*DS1-k31DS*DS3-keDS3*DS3+kon3*D3*S3-koff3*DS3;\n__DDtStateVar__[9] = _InfusionRate(9) + -kshedDM3*DM3-keDM3*DM3+kon3*D3*M3-koff3*DM3-k31DM*DM3+k13DM*(VD1/_safe_zero(VD3))*DM1;\n__DDtStateVar__[10] = _InfusionRate(10) + -keDM1*DM1-kshedDM1*DM1+kon1*D1*M1-koff1*DM1-k13DM*DM1+k31DM*(VD3/_safe_zero(VD1))*DM3;\n__DDtStateVar__[11] = _InfusionRate(11) + ksynM1-kshedM1*M1-keM1*M1+k31M*VD3/_safe_zero(VD1)*M3-k13M*M1-kon1*D1*M1+koff1*DM1;\n"));
	SET_STRING_ELT(modeln,3,mkChar("expandModel"));
	SET_STRING_ELT(model,3,mkChar("\n     D1           = AmtD1/VD1;\n     d/dt(AmtD0)  =  -ka *AmtD0;\n     d/dt(AmtD1)  =(F*ka *AmtD0/VD1      - k13D *D1 +    k31D *VD3/VD1*D3  - keD1 *D1  - kon1*D1*S1 + koff1*DS1  - kon1*D1*M1 + koff1*DM1 - k12D*D1 + k21D*VD2/VD1*D2)*VD1;\n     d/dt(D2)     =                                                                                               k12D*VD1/VD2*D1 - k21D*D2;\n     d/dt(D3)     =                        k13D *VD1/VD3*D1     - k31D*D3  - keD3 *D3  - kon3*D3*(S3+M3) + koff3*(DS3+DM3);\n     d/dt(S1)     = ksynS1+kshedM1*M1               - k13S *S1 +     k31S*VS3/VS1*S3  - keS1 *S1  - kon1*D1*S1      + koff1*DS1;\n     d/dt(S3)     = ksynS3 +kshedM3*M3   + k13S *VS1/VS3*S1     - k31S*S3  - keS3 *S3  - kon3*D3*S3      + koff3*DS3;\n     d/dt(M3)     = ksynM3 -kshedM3*M3  -k31M*M3+k13M*VD1/VD3*M1             - keM3 *M3  - kon3*D3*M3      + koff3*DM3;\n     d/dt(DS1)    =       +kshedM1*DM1               - k13DS*DS1 + k31DS*VDS3/VDS1*DS3 - keDS1*DS1 + kon1*D1*S1      - koff1*DS1;\n     d/dt(DS3)    =         kshedDM3*DM3 + k13DS*VDS1/VDS3*DS1 - k31DS*DS3 - keDS3*DS3 + kon3*D3*S3      - koff3*DS3;\n     d/dt(DM3)    = -kshedDM3*DM3  - keDM3*DM3 + kon3*D3*M3 - koff3*DM3-k31DM*DM3+k13DM*(VD1/VD3)*DM1;\n     d/dt(DM1)    = -keDM1*DM1 -kshedDM1*DM1 +kon1*D1*M1 -koff1*DM1-k13DM*DM1+k31DM*(VD3/VD1)*DM3;\n     d/dt(M1)     = ksynM1 -kshedM1*M1 -keM1*M1 +k31M*VD3/VD1*M3 -k13M*M1 -kon1*D1*M1 +koff1*DM1;\n  \n"));
	SEXP ini    = PROTECT(allocVector(REALSXP,0));
	SEXP inin   = PROTECT(allocVector(STRSXP, 0));
	SET_STRING_ELT(names,0,mkChar("params"));
	SET_VECTOR_ELT(lst,  0,params);
	SET_STRING_ELT(names,1,mkChar("lhs"));
	SET_VECTOR_ELT(lst,  1,lhs);
	SET_STRING_ELT(names,2,mkChar("state"));
	SET_VECTOR_ELT(lst,  2,state);
	SET_STRING_ELT(names,3,mkChar("trans"));
	SET_VECTOR_ELT(lst,  3,tran);
	SET_STRING_ELT(names,5,mkChar("model"));
	SET_VECTOR_ELT(lst,  5,model);
	SET_STRING_ELT(names,4,mkChar("ini"));
	SET_VECTOR_ELT(lst,  4,ini);
	SET_STRING_ELT(names,6,mkChar("md5"));
	SET_VECTOR_ELT(lst,  6,mmd5);
	SET_STRING_ELT(names,7,mkChar("podo"));
	SET_VECTOR_ELT(lst,  7,ScalarLogical(0));
	SET_STRING_ELT(names,8,mkChar("dfdy"));
	SET_VECTOR_ELT(lst,  8,dfdy);
	SET_STRING_ELT(names,9,mkChar("sens"));
	SET_VECTOR_ELT(lst,  9,sens);
	SET_STRING_ELT(names,10,mkChar("fn.ini"));
	SET_VECTOR_ELT(lst,  10,fn_ini);
	SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
	SET_STRING_ELT(mmd5,0,mkChar("6d4412ab51b586d5e8122daafa60ab4c"));
	SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
	SET_STRING_ELT(mmd5,1,mkChar(__PARSED_MD5_STR__));
	SET_STRING_ELT(trann,0,mkChar("jac"));
	SET_STRING_ELT(tran,0,mkChar("fullint"));
	SET_STRING_ELT(trann,1,mkChar("prefix"));
	SET_STRING_ELT(tran, 1,mkChar("ivsc_4cmtct_shedct_x64_"));
	SET_STRING_ELT(trann,2,mkChar("dydt"));
	SET_STRING_ELT(tran, 2,mkChar("ivsc_4cmtct_shedct_x64_dydt"));
	SET_STRING_ELT(trann,3,mkChar("calc_jac"));
	SET_STRING_ELT(tran, 3,mkChar("ivsc_4cmtct_shedct_x64_calc_jac"));
	SET_STRING_ELT(trann,4,mkChar("calc_lhs"));
	SET_STRING_ELT(tran, 4,mkChar("ivsc_4cmtct_shedct_x64_calc_lhs"));
	SET_STRING_ELT(trann,5,mkChar("model_vars"));
	SET_STRING_ELT(tran, 5,mkChar("ivsc_4cmtct_shedct_x64_model_vars"));
	SET_STRING_ELT(trann,6,mkChar("ode_solver"));
	SET_STRING_ELT(tran, 6,mkChar("ivsc_4cmtct_shedct_x64_ode_solver"));
	SET_STRING_ELT(trann,7,mkChar("ode_solver_sexp"));
	SET_STRING_ELT(tran, 7,mkChar("ivsc_4cmtct_shedct_x64_ode_solver_sexp"));
	SET_STRING_ELT(trann,8,mkChar("ode_solver_0_6"));
	SET_STRING_ELT(tran, 8,mkChar("ivsc_4cmtct_shedct_x64_ode_solver_0_6"));
	SET_STRING_ELT(trann,9,mkChar("ode_solver_focei_eta"));
	SET_STRING_ELT(tran, 9,mkChar("ivsc_4cmtct_shedct_x64_ode_solver_focei_eta"));
	SET_STRING_ELT(trann,10,mkChar("ode_solver_ptr"));
	SET_STRING_ELT(tran, 10,mkChar("ivsc_4cmtct_shedct_x64_ode_solver_ptr"));
	SET_STRING_ELT(trann,11,mkChar("inis"));
	SET_STRING_ELT(tran, 11,mkChar("ivsc_4cmtct_shedct_x64_inis"));
	setAttrib(tran, R_NamesSymbol, trann);
	setAttrib(mmd5, R_NamesSymbol, mmd5n);
	setAttrib(model, R_NamesSymbol, modeln);
	setAttrib(ini, R_NamesSymbol, inin);
	setAttrib(lst, R_NamesSymbol, names);
	UNPROTECT(16);
	return lst;
}
void __ODE_SOLVER__(
                    int *neq,
                    double *theta,      //order:
                    double *time,
                    int *evid,
                    int *ntime,
                    double *inits,
                    double *dose,
                    double *ret,
                    double *atol,
                    double *rtol,
                    int *stiff,
                    int *transit_abs,
                    int *nlhs,
                    double *lhs,
                    int *rc
                    ){
  // Backward compatible ode solver for 0.5* C interface
  __ODE_SOLVER_PTR__();
  _old_c(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

void __ODE_SOLVER_0_6__(int *neq,
                        double *theta,  //order:
                        double *time,
                        int *evid,
                        int *ntime,
                        double *inits,
                        double *dose,
                        double *ret,
                        double *atol,
                        double *rtol,
                        int *stiff,
                        int *transit_abs,
                        int *nlhs,
                        double *lhs,
                        int *rc,
                        double hmin,
                        double hmax,
                        double h0,
                        int mxordn,
                        int mxords,
                        int mxstep) {
  // Backward compatible ode solver for 0.5* C interface
  __ODE_SOLVER_PTR__();
  _c_0_6(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc,
	hmin, hmax, h0, mxordn, mxords, mxstep);
}

extern void __ODE_SOLVER_PTR__  (){
  _assign_fn_pointers(__DYDT__ , __CALC_LHS__ , __CALC_JAC__, __INIS__, __JT__ , __MF__,
#ifdef __DEBUG__
                      1
#else
                      0
#endif
                      );
}

extern SEXP __ODE_SOLVER_SEXP__ (// Parameters
                                 SEXP sexp_theta,
                                 SEXP sexp_inits,
                                 SEXP sexp_lhs,
				 // Events
				 SEXP sexp_time,
				 SEXP sexp_evid,
				 SEXP sexp_dose,
				 // Covariates
				 SEXP sexp_pcov,
				 SEXP sexp_cov,
				 SEXP sexp_locf,
				 // Solver Options
				 SEXP sexp_atol,
				 SEXP sexp_rtol,
				 SEXP sexp_hmin,
				 SEXP sexp_hmax,
				 SEXP sexp_h0,
				 SEXP sexp_mxordn,
				 SEXP sexp_mxords,
				 SEXP sexp_mx,
				 SEXP sexp_stiff,
				 SEXP sexp_transit_abs,
				 // Object Creation
				 SEXP sexp_object,
				 SEXP sexp_extra_args){
  RxODE_ode_solver ode_solver=(RxODE_ode_solver) R_GetCCallable("RxODE","RxODE_ode_solver");
  __ODE_SOLVER_PTR__();
  ode_solver(sexp_theta,sexp_inits,sexp_lhs,sexp_time,sexp_evid,sexp_dose,sexp_pcov,sexp_cov,sexp_locf,sexp_atol,
	     sexp_rtol,sexp_hmin,sexp_hmax,sexp_h0,sexp_mxordn,sexp_mxords,sexp_mx,sexp_stiff,sexp_transit_abs,
	     sexp_object,sexp_extra_args);
}

//Initilize the dll to match RxODE's calls
void __R_INIT__ (DllInfo *info){
  // Get the RxODE calling interfaces
  _InfusionRate   = (RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRate");
  _update_par_ptr = (RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptr");
  _par_ptr = (RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptr");
  _dadt_counter_val = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_val");
  _jac_counter_val  = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_val");
  _dadt_counter_inc = (RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_inc");
  _jac_counter_inc  = (RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_inc");
  podo  = (RxODE_val) R_GetCCallable("RxODE","RxODE_podo");
  tlast = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlast");
  factorial=(RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _transit3 = (RxODE_transit3) R_GetCCallable("RxODE","RxODE_transit3");
  _transit4 = (RxODE_transit4) R_GetCCallable("RxODE","RxODE_transit4");
  _safe_log =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  _safe_zero =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  _as_zero =(RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _assign_fn_pointers=(RxODE_assign_fn_pointers) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","RxODE_ode_solver_old_c");
  _c_0_6 = (RxODE_ode_solver_0_6_c)R_GetCCallable("RxODE","RxODE_ode_solver_0_6_c");
  // Register the outside functions
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_STR__,       (DL_FUNC) __ODE_SOLVER__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_SEXP_STR__,  (DL_FUNC) __ODE_SOLVER_SEXP__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_0_6_STR__,   (DL_FUNC) __ODE_SOLVER_0_6__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_PTR_STR__,   (DL_FUNC) __ODE_SOLVER_PTR__);

  /* R_CallMethodDef callMethods[]  = { */
  /*   {__ODE_SOLVER_PTR_STR__, (DL_FUNC) &__ODE_SOLVER_PTR__, 0}, */
  /*   {__ODE_SOLVER_SEXP_STR__, (DL_FUNC) &__ODE_SOLVER_SEXP__, 21}, */
  /*   {NULL, NULL, 0} */
  /* }; */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info,TRUE);
  // Register the function pointers so if someone directly calls the
  // ode solvers directly, they use the last loaded RxODE model.
  __ODE_SOLVER_PTR__();
}
