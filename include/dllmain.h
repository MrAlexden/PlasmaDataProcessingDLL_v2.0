#pragma once

#ifndef MY_PLASMA_DLL_H
#define MY_PLASMA_DLL_H

#define WIN32_LEAN_AND_MEAN					// исключает редко используемые компоненты из заголовков Windows
#define ERR(f) if (err = (f), err < 0) goto EndBlock;

/* Файлы заголовков Windows */
#include <windows.h>				// for windows interface

/* Директивы компиляции */
#define INSIDE_MY_PDPv2_PROJECT				// директива, чтобы не включать данную библиотеку в код (см. interface.h)
#define GPU_COMPUTATION				1		// opencl gpu computation wil be used to speed up processing
#define STRING_EXPRESSION			1		// definition which means that the computation of approximation will be performed via string parsing and evaluating
#define MESSAGEBOX					1		// definition which determines whether the 'MessageBoxA' will be used in case of error or not
#define SHOW_PROGRESS				1		// show progress bar in a pop up window or not

/* Макросы */
#if MESSAGEBOX > 0
	#define SHOUT_ERR(f) MessageBoxA(NULL, ERR_GetErrorDescription(f).c_str(), "Error!", MB_ICONINFORMATION | MB_OK), err_num(f)
#else
	#define SHOUT_ERR(f) err_num(f)
#endif

/* Файл взаимодействия с внешним миром */
#include "interface.h"

/* Мои файлы заголовков */
#include "error_codes.h"
#include "../myspace/mynamespace.h"				// usefull functions which can be used by any file
#include "workspace.h"							// global workspace
#include "results_container.h"					// global results
#include "usefull_signal_extractor.h"
#include "../muparser/include/muParser.h"
#include "../opencl/gpu_compute.h"
#include "../levenberg_marquardt/levenberg_marquardt.h"
#include "support_functions.h"
#include "../progress_window/resource.h"
#include "../progress_window/progress_window.h"
//#include "../diagnostics/langmuir_probe.h"

/* Глобальные переменные */
extern mu::Parser expression_parser;			// string expression parser and evaluator class

#endif /* MY_PLASMA_DLL_H */
