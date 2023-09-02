#pragma once

#ifndef SUPPORT_FUNCTIONS_H
#define SUPPORT_FUNCTIONS_H

#include <vector>       // for vector
#include <algorithm>    // for min_element|max_element
#include <functional>   // for function

extern mu::Parser expression_parser;

/* makes linerization of a segment of an input voltage (ramp or sine shaped).
the result is set on results.rump */
inline int prepare_ramp()
{

#define IVOLT workspace.ramp_peaks_indices	// это поле в workspace должно быть заполнено заранее
#define LEN workspace.one_segment_length	// это поле в workspace должно быть заполнено заранее
#define L workspace.left_cut_off			// это поле в workspace должно быть заполнено заранее
#define R workspace.right_cut_off			// это поле в workspace должно быть заполнено заранее

	std::vector <double> voltage(workspace.voltage), ramp;

	int todel = (int)abs(IVOLT[2] - IVOLT[1]) % 100; double length_between_voltage_peaks = abs(IVOLT[2] - IVOLT[1]) + ((double)todel / 100 <= 0.5 ? 0 : 100) - todel;

	/* +5 потому что мы нашли пики пилы, а в этой функции пила строится от нижней точки,
	соответственно нужно спуститься вниз пилы (пила выглядит как "/|/|/")
	--------------------------------------\/--------------------------------------------------\/ */
	int ind_beg = (L >= 0.2) ? (IVOLT[1] + 5) + L * length_between_voltage_peaks : (IVOLT[1] + 5) + 0.2 * length_between_voltage_peaks,
		ind_end = (R >= 0.2) ? (IVOLT[1] + 5) + (1 - R) * length_between_voltage_peaks : (IVOLT[1] + 5) + (1 - 0.2) * length_between_voltage_peaks;
	ramp.resize(LEN * (1 - L - R));

	/* определяем приращение в одной точке пилы от предыдущей */
	double amp_beg = voltage[ind_beg],
		amp_end = voltage[ind_end],
		segment_length = amp_end - amp_beg,
		segment_points_amount = abs(ind_beg - ind_end),
		transform_factor = (double)LEN / length_between_voltage_peaks,
		delta = segment_length / (segment_points_amount * transform_factor); // <------ ВОЗМОЖНО ДОЛЖНО БЫТЬ abs(...)

#define iterator std::vector <double>::iterator
	std::function<iterator(iterator, iterator)> minormax;

	if (amp_beg > amp_end) /* пила выглядит как '\' */
	{
		minormax = [](iterator begin, iterator end) { return max_element(begin, end); };
	}
	else /* пила выглядит как '/' */
	{
		minormax = [](iterator begin, iterator end) { return min_element(begin, end); };
	}
#undef iterator

	/* Нужен держатель иначе дельта функция ругается на константный итератор */
	std::vector <double> voltage_holder(voltage.begin() + (IVOLT[1] + 5), voltage.begin() + (IVOLT[1] + 5) + length_between_voltage_peaks);
	int index_of_minormax = minormax(voltage_holder.begin(), voltage_holder.end()) - voltage_holder.begin();

	/* после этого пила всегда будет выглядеть как '/' */
	ramp[0] = *(voltage.begin() + (IVOLT[1] + 5) + index_of_minormax + L * length_between_voltage_peaks);
	for (int i = 1; i < ramp.size(); ++i)
		ramp[i] = ramp[i - 1] + delta;
	results.set_ramp(ramp);

	return 0;

#undef IVOLT
#undef LEN
#undef L
#undef R

};

/* function to pass into lev-marq if string expression being used */
inline double s_fx(const double x, const std::vector<double>& vparams = {})
{
	if (!vparams.empty() && workspace.variables_values.size() > 1) memcpy(&workspace.variables_values[1], vparams.data(), sizeof(double) * vparams.size());

	workspace.variables_values[0] = x;

	return expression_parser.Eval();
};

#endif