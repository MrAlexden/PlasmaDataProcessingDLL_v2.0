#pragma once

#if GPU_COMPUTATION > 0
    #define STRING_EXPRESSION 1
#endif

#ifndef WORKSPACE_H
#define WORKSPACE_H

#include <windows.h>    // HWND / HINSTANCE
#include <type_traits>  // for std::is_convertible
#include <cstring>      // for std::memcpy

/* структура в которой хранятся все данные, используемые в ходе работы. Очищается при отсоединении процесса с этой dll */
template <typename T>
    requires std::is_convertible_v<T, double>
struct Data
{

public:

    std::vector<T> voltage;                                     // measured voltage
    std::vector<T> signal;                                      // measured signal

    int diagnostic = 0;                                         // type of diagnostic (uses enum "diagnostics")

    T left_cut_off = 0.0;                                       // a part of a segment which will be ignored in noise deleting and calculation (at the beginning of a segment)
    T right_cut_off = 0.0;                                      // a part of a segment which will be ignored in noise deleting and calculation (at the end of a segment)
    T begin_time = 0.0;                                         // the time of plasma-discharge beginning in seconds
    T end_time = 0.0;                                           // the time of plasma-discharge ending in seconds
    T probe_area = 0.0;                                         // area of the probe plasma-interface area in m^2
    T linerize_probe_part = 0.0;                                // a part of the probe segment which linear approximatiob will be applied to (from 0 to 0.9)
    T filtration_window_width = 0.0;                            // width of Savitsky_Golay filter which will be applied to the signal
    T measuring_resistance = 0.0;                               // the resistnce on which the voltage difference profile is measured (signal)
    T amplification_factor = 0.0;                               // the amplification parameter of a voltage amplifier before which the voltage writen to A-to-D
    T ion_mass = 0.0;                                           // mass of fuel ion in kilogramms

    size_t one_segment_length = 0;                              // original length of one segment (number of points, without cutoff) of the input signal
    size_t ramp_freq = 0;                                       // frequancy of voltage in Herz
    size_t niter = 0;                                           // number of Levenberg_Mrqurdt least squares reduction steps

    std::vector <long long> ramp_peaks_indices;                 // std::vector which stroes the indices of local peaks of saw-shaped voltage
    std::vector <long long> segments_beginning_indices;         // std::vector which stores the indices of segments in which the signal being separated by saw-shaped voltage

    HINSTANCE hInstThisDll = NULL;                              // HINSTANCE of this DLL/LIB. Needed to create the progress bar window
    HWND progress_window = NULL;                                // Holder of the progress bar window

    std::string s_expression;                                   // std::string mathematical expression which the input signal will be approximated with
    std::vector <std::string> names_of_expression_variables;    // std::string names of variables of s_expression (first variable is independent)
    std::vector <T> variables_values;                           // the vector which is linking to a variables names an being used as a temporary storage of value of each variable (first variable is independent)
    std::vector <bool> is_variable_fixed;                       // is initial value fixed and not need to be changed or not
    
    ~Data()
    {
        voltage.~vector();
        signal.~vector();

        diagnostic = 0;

        left_cut_off = 0.0;
        right_cut_off = 0.0;
        begin_time = 0.0;
        end_time = 0.0;
        probe_area = 0.0;
        linerize_probe_part = 0.0;
        filtration_window_width = 0.0;
        measuring_resistance = 0.0;
        amplification_factor = 0.0;
        ion_mass = 0.0;

        one_segment_length = 0;
        ramp_freq = 0;
        niter = 0;

        ramp_peaks_indices.~vector();
        segments_beginning_indices.~vector();

        hInstThisDll = NULL;
        progress_window = NULL;

        s_expression.~basic_string();
        names_of_expression_variables.~vector();
        variables_values.~vector();
        is_variable_fixed.~vector();

        return;
    }

public:

    void clean()
    {
        this->~Data();
    }

    bool is_empty() const
    {
        return (voltage.empty() && signal.empty()) ? true : false;
    }
};

extern Data<double> workspace;

#endif /* WORKSPACE_H */
