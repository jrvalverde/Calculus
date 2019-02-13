/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	COMPLEX.H
*
*	Definitions to use with the module of complex calculus
*
*	Designed by: J. R. Valverde Carrillo
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

/*-------------------------------------------
	Definition of complex as a
	data structure
--------------------------------------------*/

typedef struct
    {
	double r, i;
    } complex;

/*-------------------------------------------
	Data structure for polar form
	manipulation
--------------------------------------------*/

typedef struct
    {
	double mod,	/* mod = sqrt(r * r + i * i)	*/
	       ang;	/* ang = atan(i / r)		*/
    } polar;

extern void c_set(complex *, real, real);
extern void p_set(polar *, real, real);
extern void c_topolar(polar *, complex *);

extern bool c_comp(complex *, complex *);

extern void c_copy(complex*, complex *);
extern void c_add(complex *, complex *, complex *);
extern void c_subs(complex *, complex *, complex *);
extern void c_mult(complex *, complex *, complex *);
extern void c_conj(complex *, complex *);
extern void c_inv(complex *, complex *);
extern void c_div(complex *, complex *, complex *);
extern double c_mod(real *, complex *);
extern void c_exp(complex *, complex *);
extern void c_log(complex *, complex *);
extern void c_logx(real, complex *, complex *);
extern void c_pow(complex *, complex *, complex *);
extern void c_radix(complex *, complex *, complex *);
extern void c_rpow(complex *, complex *, real);
extern void c_rradix(complex *, complex *, real);
extern void c_sin(complex *, complex *);
extern void c_cos(complex *, complex *);
extern void c_tan(real *, complex *);
extern void p_tocomplex(complex *, polar *);
extern void p_copy(polar *, polar *);
extern void p_mult(polar *, polar *, polar *);
extern void p_div(polar *, polar *, polar *);
extern void p_mod(real *, polar *);
extern void p_rpow(polar *, polar *, real);
extern void p_rradix(polar *, polar *, real);
extern int ComplError();
extern int getComplError();
extern void resetComplError();
extern void setComplError(int);
