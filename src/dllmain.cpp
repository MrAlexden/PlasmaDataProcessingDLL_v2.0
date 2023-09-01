#include "../headers/dllmain.h"

Data<double> workspace;                 // structure all the input data and global variables are stored in
Results<double> results;                // the class all the results are stored in
mu::Parser expression_parser;           // string expression parser and evaluator class
static int execution_order = NULL;      // bitwise mask which used to determine whether the module was executed or not

namespace myplasmadll
{

    /* sets the voltage corresponding to the data */
    int set_voltage(_In_ const double* in, _In_ const int size)
    {
        if (in == nullptr || size == 0) return SHOUT_ERR(ERR_BadInputVecs);

        workspace.set_voltage(in, size);

        /* now we know that the voltge is set */
        execution_order |= 1 << 3;

        return 0;
    }

    /* sets the data which needs to be processed into workspace */
    int set_signal(_In_ const double* in, _In_ const int size)
    {
        if (in == nullptr || size == 0) return SHOUT_ERR(ERR_BadInputVecs);

        workspace.set_signal(in, size);

        /* now we know that the signal is set */
        execution_order |= 1 << 2;

        return 0;
    }

    /* sets the additional parameters (1st param) for calculation (THE SIZE OF INPUT IS ALWAYS 12),
    2nd param (optional) sets the initial values for string expression variables,
    3rd param (optional) number of initial values */
    int set_params(_In_ const double* in, _In_ const double* init_var_vals, _In_ const int init_vals_size)
    {
        if (in == nullptr) return SHOUT_ERR(ERR_BadInputVecs);
        if ((in[1] < 0 || in[2] < 0 || in[1] > in[2])
            || in[6] < 0
            || in[7] <= 0
            || in[8] == 0 // может быть < 0 для переворачивания сигнала
            || in[9] == 0 // может быть < 0 для переворачивания сигнала
            || in[10] <= 0
            || in[11] <= 0) return SHOUT_ERR(ERR_ZeroInputVals);
        if (in[3] + in[4] > 0.9
            || in[3] < 0
            || in[4] < 0) return SHOUT_ERR(ERR_BadCutOff);
        if (in[5] < 0 || in[5] > 0.9) return SHOUT_ERR(ERR_BadLinFit);
#if STRING_EXPRESSION > 0
        if (init_var_vals != nullptr && init_vals_size != NULL && (execution_order & 1) == 0) return SHOUT_ERR(ERR_ExpressionNotSet); // string expression not set
        else if (init_var_vals != nullptr && init_vals_size != NULL) memcpy(&workspace.variables_values[1], init_var_vals, sizeof(double) * init_vals_size);
#else
        if (init_var_vals != nullptr && init_vals_size != NULL)
        {
            if (workspace.variables_values.size() < 2) workspace.variables_values.resize(init_vals_size + 1);
            memcpy(&workspace.variables_values[1], init_var_vals, sizeof(double) * init_vals_size);
        }
#endif

        workspace.probe_area = in[0];		            // площадь поверхности зонда
        workspace.begin_time = in[1];		            // время начала обработки (опционально)
        workspace.end_time = in[2];		                // время конца обработки (опционально)
        workspace.left_cut_off = in[3];		            // часть точек отсечки слева
        workspace.right_cut_off = in[4];	            // часть точек отсечки справа
        workspace.linerize_probe_part = in[5];	        // часть точек линейно аппроксимации
        workspace.filtration_window_width = in[6];	    // ширина окна фильтрации сигнала
        workspace.ramp_freq = (int)in[7];	            // частота пилы
        workspace.measuring_resistance = in[8];         // сопротивление, через которое измерялось падение напряжения
        workspace.amplification_factor = in[9];         // коэффициент усиления напряжения на АЦП
        workspace.ion_mass = in[10];					// масса иона рабочего вещества
        workspace.niter = (int)in[11];				    // количество итераций аппроксимации(сильно влияет на скорость работы программы)

        /* now we know that all the additional parameters are set */
        execution_order |= 1 << 1;

        return 0;
    }

    /* sets up passed string as expression (1st string),
    sets up passed string as parameters of this expression (2nd string),
    parameters should be separated from each other by delimeter (example: _,;|-),
	delimeter can be set as third (optional) string (by default delimeters are set of : "\t \n,./*|_+-:;'()!%^&=<>[]"),
	4th (optional) parameter of this function is independent variable of expression (default name is : 'x') */
    int set_approximation_equation( _In_ const char* s_expression,      _In_ const unsigned int expression_s_size,
                                    _In_ const char* s_variables,       _In_ const unsigned int variables_s_size,
                                    _In_opt_ const char* s_delimeter,   _In_opt_ const unsigned int delimeter_s_size,
                                    _In_opt_ const char independent_variable)
    {
        /* getting variables names */
        std::vector<std::string> params_names = myspace::separate_string_by_delimeter(std::string(s_variables, s_variables + variables_s_size),
                                                    s_delimeter == nullptr ? "" : std::string(s_delimeter, s_delimeter + delimeter_s_size));
        /* independent variable */
        params_names.insert(params_names.begin(), std::string(&independent_variable, &independent_variable + 1));
        /* duplicate name of independent variable check */
        params_names = myspace::erase_duplicates(params_names);

        /* set up expression and check for errors */
        try
        {
            expression_parser.SetExpr(myspace::utf8_decode(std::string(s_expression, s_expression + expression_s_size)));

            workspace.variables_values = std::vector<double>(params_names.size()); // params values are all = 0, just to execute .Eval() and find mistakes, if any exist

            /* initialising parameters */
            for (int i = 0; i < workspace.variables_values.size(); ++i)
                expression_parser.DefineVar(myspace::utf8_decode(params_names[i]), &workspace.variables_values[i]);

            /* .Eval() will throw an exception if any mistake in passed expression appear */
            double result = expression_parser.Eval();
        }
        catch (const std::exception& ex)
        {
            return 
#if MESSAGEBOX > 0
                MessageBoxA(NULL, ex.what(), "Error at passed expression!", MB_ICONINFORMATION | MB_OK), 
#endif
                ERR_BadExpression;
        }

        /* if no error occured -> seting up std::string expression */
        workspace.s_expression = std::string(s_expression, s_expression + expression_s_size);
        workspace.names_of_expression_variables = params_names;

        /* now we know that the expression is set */
        execution_order |= 1;

        return 0;
    }

    /* extract signal out from the whole dataset and perform calculation with particular type of diagnostic */
    int calculate(_In_ const int calculation_type)
    {
#if STRING_EXPRESSION > 0
        if (execution_order < 15) return SHOUT_ERR(ERR_NoetAllDataSet); // means that bitwise is less than 1.1.1.1
#else
        if (execution_order < 14) return SHOUT_ERR(ERR_NoetAllDataSet); // means that bitwise is less than 1.1.1.0
#endif

        int err = 0;

        workspace.voltage *= workspace.amplification_factor;
        workspace.signal /= workspace.measuring_resistance;

        ERR(erase_noise());
        ERR(prepare_ramp());

#define ISIGN workspace.segments_beginning_indices
#define LEN results.get_SizeOfSegment()

#if SHOW_PROGRESS > 0
        {
            ThreadArgs args; HANDLE progress_window_thread;
            if (progress_window_thread = CreateThread(NULL, NULL, (LPTHREAD_START_ROUTINE)&DialogBoxParamWrapper, &args, NULL, NULL), progress_window_thread)
                WaitForSingleObject(progress_window_thread, 10), SendMessage(GetDlgItem(workspace.progress_window, IDC_PROGRESS1), PBM_SETRANGE, 0, MAKELPARAM(0, ISIGN.size() - 1));
#endif

            for (int n = 0; n < ISIGN.size(); ++n)
            {
#if SHOW_PROGRESS > 0
                if (n % 4 == 0) // отправляю мэсседж раз в 4 итерации чтобы меньше нагружать
                    SendMessage(GetDlgItem(workspace.progress_window, IDC_PROGRESS1), PBM_SETPOS, n, NULL),
                    SetDlgItemText(workspace.progress_window, IDC_STATIC, std::wstring(L"In Progress... " + std::to_wstring(int(((double)n / ISIGN.size()) * 100)) + L" %").c_str());
#endif

                ERR(make_one_segment(calculation_type, 
                    results.get_ramp().data(), LEN,
                    workspace.signal.data() + ISIGN[n], LEN));
            }

#if SHOW_PROGRESS > 0
            if (progress_window_thread)
            {
                EndDialog(workspace.progress_window, NULL);
                TerminateThread(progress_window_thread, -1);
                CloseHandle(progress_window_thread);
            }
        }
#endif

    EndBlock:

#if MESSAGEBOX > 0
        if (err < 0) MessageBoxA(NULL, ERR_GetErrorDescription(err).c_str(), "Error!", MB_ICONINFORMATION | MB_OK);
#endif

        return err;

#undef ISIGN
#undef LEN
    }

    /* perform calculation at one segment. all the data must be cleaned out of noise and prepared already */
    int make_one_segment(_In_ const int calculation_type,
        _In_ const double* vx, _In_ const unsigned int x_size,
        _In_ const double* vy, _In_ const unsigned int y_size)
    {
#if STRING_EXPRESSION > 0
        if (execution_order < 3) return SHOUT_ERR(ERR_NoetAllDataSet); // means that bitwise is less than 0.0.1.1
#else
        if (execution_order < 2) return SHOUT_ERR(ERR_NoetAllDataSet); // means that bitwise is less than 0.0.1.0
#endif

        int err = 0;

        std::vector <double> vparams = (workspace.variables_values.size() > 1) ? 
            std::vector <double>(&workspace.variables_values[1], &workspace.variables_values[1] + workspace.variables_values.size() - 1) :
            std::vector <double>();

#if STRING_EXPRESSION > 0
        ERR(mystringcompute::LevenbergMarquardt(std::vector<double>(vx, vx + x_size), std::vector<double>(vy, vy + y_size), vparams));
#else
        ERR(myspace::LevenbergMarquardt(std::vector<double>(vx, vx + x_size), std::vector<double>(vy, vy + y_size), vparams, fx));
#endif

    EndBlock:

        return 0;

#undef NUM_OF_SEG
    }

    /* returns the size of stored data as bytes amount */
    //DLL_API int getbytes()
    //{
    //    return workspace.get_bytes();
    //}

    /* returns the size of stored data as int number */
    //DLL_API int getsize()
    //{
    //    return workspace.size();
    //}

    /* returns the pointer to the stored data */
    //DLL_API double* getdata()
    //{
    //    return workspace.data();
    //}
}

BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
    switch (ul_reason_for_call)
    {
        case DLL_PROCESS_ATTACH:
        {
#ifdef DEBUG
            MessageBoxA(NULL, "process attached", "", MB_ICONINFORMATION | MB_OK);
#endif
#if SHOW_PROGRESS > 0
            workspace.hInstThisDll = hModule;
#endif
            break;
        }
        case DLL_PROCESS_DETACH:
        {
#ifdef DEBUG
            MessageBoxA(NULL, "process detached", "", MB_ICONINFORMATION | MB_OK);
#endif
#if GPU_COMPUTATION > 0
            mygpucompute::clear_gpu();
#endif
            workspace.clean();
            results.~Results();

#if STRING_EXPRESSION > 0
            expression_parser.ClearVar();
            expression_parser.ClearFun();
            expression_parser.ClearConst();
            expression_parser.ClearInfixOprt();
            expression_parser.ClearPostfixOprt();
            expression_parser.ClearOprt();
#endif

            execution_order = NULL;

            break;
        }
        case DLL_THREAD_ATTACH:
            break;
        case DLL_THREAD_DETACH:
            break;
        default:
            __assume(0);
    }
    return TRUE;
}