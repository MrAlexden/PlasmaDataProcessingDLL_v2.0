#pragma once

/*  This file exists in case u need to adjust some settings
* while including this library as .lib in ur project */

#ifndef INTERFACE_H
#define INTERFACE_H

#define DYNAMIC		// outside my project (needed to enter DLL_API definition directive)
/*         ^---------- manually define STATIC if u include this library as .lib or
							define DYNAMIC if u use it as .dll */
//#define DLLEXPORT	// uncomment this definition if u want _declspec(dllexport)

/* directives below works only if u include this interface file in ur project
(as soon as in mine project INSIDE_MY_PDPv2_PROJECT is defined in proper compilation purposes) */
#ifndef INSIDE_MY_PDPv2_PROJECT
	#ifdef STATIC // you need to define this var by yourself (if u include this lib static, not .dll)
		#ifdef _DEBUG
			#pragma comment(lib, __FILE__"\\..\\x64\\Debug\\PlasmaDataProcessingDLL_v2.0.lib")
		#else
			#pragma comment(lib, __FILE__"\\..\\x64\\Release\\PlasmaDataProcessingDLL_v2.0.lib")
		#endif
	#endif
#endif

namespace myplasmadll
{

#ifndef DLL_API
	#ifdef STATIC
		#define DLL_API
	#elif defined (DYNAMIC)
		#ifdef DLLEXPORT
			#define DLL_API _declspec(dllexport)
		#else
			#define DLL_API _declspec(dllimport)
		#endif
	#endif
#endif

#ifdef __cplusplus
	EXTERN_C // defined in <windows.h>
	{
#endif

	/* sets the voltage corresponding to the data */
	DLL_API int set_voltage(_In_ const double* in, _In_ const int size);

	/* sets the data which needs to be processed into workspace */
	DLL_API int set_signal(_In_ const double* in, _In_ const int size);

	/*  1st param type of diagnostic (calculation process depends of this param)
		2nd param sets the additional parameters for calculation (THE SIZE OF INPUT IS ALWAYS 12),
		3rd param (optional) sets the initial values for string expression variables,
		4th param (optional) number of initial values,
		5th param (optional) tells if initial value is fixed and not need to be changed (size is same as 2nd param) */
	DLL_API int set_params(	_In_ const unsigned int calculation_type,
							_In_ const double* in, 
							_In_opt_ const double* init_var_vals = nullptr, 
							_In_opt_ const int init_vals_size = NULL, 
							_In_opt_ const int* is_init_val_fixed = nullptr);

	/* sets up passed string as expression (1st string),
		sets up passed string as parameters of this expression (2nd string),
		parameters should be separated from each other by delimeter (4 exmpl: _,;|-),
		delimeter can be set as 3rd (optional) string (by default delimeters are set of : "\t \n,./*|_+-:;'()!%^&=<>[]"),
		4th (optional) parameter of this function is independent variable of expression (default name is : 'x') */
	DLL_API int set_approximation_equation(	_In_ const char* s_expression,					_In_ const unsigned int expression_s_size,
											_In_ const char* s_variables,					_In_ const unsigned int variables_s_size,
											_In_opt_ const char* s_delimeter = nullptr,		_In_opt_ const unsigned int delimeter_s_size = NULL,
											_In_opt_ const char independent_variable = 'x');

	/* extract signal out from the whole dataset and perform calculation with particular type of diagnostic */
	DLL_API int calculate();

	/* perform calculation at one segment. all the data must be cleaned out of noise and prepared already.
		this is the function for extern usage (wrapper of _make_one_segment) */
	DLL_API int make_one_segment(	_In_ const double* vx, _In_ const unsigned int x_size,
									_In_ const double* vy, _In_ const unsigned int y_size);

#ifdef __cplusplus
	}
#endif

	typedef enum
	{
		LANGMUIR_PROBE,
		DOUBLE_LANGMUIR_PROBE,
		RETARDING_FIELD_ANALYZER,
		RADIAL_CYLINDRICAL_ANALYZER,
		MAGNETIC_PLASMA_ANALYSER
	} diagnostics;
}

#endif /* INTERFACE_H */