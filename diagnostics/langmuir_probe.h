#pragma once

#ifndef LANGMUIR_PROBE
#define LANGMUIR_PROBE

#define ERR(f) if (err = (f), err < 0) goto EndBlock;

#include <vector>
#include "../include/error_codes.h"
#include "../myspace/mynamespace.h"

#define LANGMUIR_PROBE_RESULTS_NUMBER	2
#define FILTRATION						0
#define APPROXIMATION					1
#if STRING_EXPRESSION > 0
	#define EXPRESSION_VARIABELS_NUMBER 0
#else
	#define EXPRESSION_VARIABELS_NUMBER 4
#endif

template <typename T>
	requires std::is_convertible_v<T, double>
class langmuir_probe
{
private:

	std::vector <T> vx, vy;
	std::vector <bool>				is_fixed(EXPRESSION_VARIABELS_NUMBER);
	std::vector <T>					exp_vars(EXPRESSION_VARIABELS_NUMBER);
	std::vector <std::vector <T>>	_results(LANGMUIR_PROBE_RESULTS_NUMBER);

	T probe_area;
	T linerize_probe_part;
	T filtration_window_width;
	T ion_mass;

	constexpr T fx(const T x, const std::vector<T> & vparams) const 
	{
#define A 0
#define B 1
		return vparams[A] + vparams[B] * x;
	}

	int initialize_variables(const std::vector<T>& vparams = {})
	{
#if STRING_EXPRESSION > 0
		exp_vars.assign(vparams.begin(), vparams.end());
#else
		if (vparams.empty())
		{
			exp_vars[0] = vy[0];
			exp_vars[1] = 1;
			exp_vars[3] = 15;
			exp_vars[2] = abs((vy.back() - exp_vars[0] - exp_vars[1] * vx.back()) / (exp(vx.back() / exp_vars[3])));
		}
		else
		{
			exp_vars.assign(vparams.begin(), vparams.end());
		}
#endif

		return 0;
	}

public:

	langmuir_probe()
	{
		initialize_variables();
	}

	~langmuir_probe()
	{}

	int calculate()
	{
		int err = 0;

		if (filtration_window_width > 0)
			_results[FILTRATION] = myspace::sg_smooth(vy, vx.size() * (filtration_window_width / 2), (5/*OriginLab poly order*/ - 1));

#if STRING_EXPRESSION <= 0
		if (linerize_probe_part > 0) // TODO: подумать что делать если стринговое уравнение и внем нет линейной части
		{
			std::vector<T> params_buff(2);
			ERR(myspace::LevenbergMarquardt(std::vector<double>(vx.begin(), vx.begin() + vx.size() * linerize_probe_part),
											std::vector<double>(vy.begin(), vy.begin() + vx.size() * linerize_probe_part),
											params_buff, fx));

			exp_vars[0] = params_buff[0];
			exp_vars[1] = params_buff[1];

			is_fixed[0] = is_fixed[1] = true;
		}
		else
		{
			exp_vars[0] = vy[0];
			exp_vars[1] = 1;
		}
#endif

#if STRING_EXPRESSION > 0
		ERR(mystringcompute::LevenbergMarquardt(std::vector<double>(vx), _results[FILTRATION], vparams));
#else

#endif

	EndBlock:

		return err;
	}
}

static inline init()
{

}

inline int langmuir_probe(std::vector<double> & vx, std::vector<double> & vy)
{

#if STRING_EXPRESSION > 0

#else

#endif

	int i = 0;
	vector <bool> is_fixed = { false, false, false, false };
	vector <myflo> vParams(4);

	vres.resize(vx.size());
	vfilt.resize(vx.size());
	vdiff.resize(vx.size());

	if (filtS != 0)
	{
		myspace::sg_smooth(vy, vx.size() * (filtS / 2), (5/*OriginLab poly order*/ - 1));
		vy = vfilt;
	}

	if (linfitP != 0)
	{
		vector <myflo> vx, vy;

		vx.assign(vx.begin(), vx.begin() + vx.size() * linfitP);
		vy.assign(vy.begin(), vy.begin() + vx.size() * linfitP);

		vector <myflo> vAB = linear_fit(vx, vy);

		vParams[0] = vAB[0];
		vParams[1] = vAB[1];

		is_fixed[0] = is_fixed[1] = true;
	}
	else
	{
		vParams[0] = vy[0];
		vParams[1] = 1;
	}

	vParams[3] = 15;
	vParams[2] = abs((vy.back() - vParams[0] - vParams[1] * vx.back()) / (exp(vx.back() / vParams[3])));
	//vParams[2] = -vParams[0];

	if (linfitP != 0)
	{
		vector <myflo> vx, vy;

		vx.assign(vx.begin() + vx.size() * linfitP, vx.end());
		vy.assign(vy.begin() + vx.size() * linfitP, vy.end());

		levmarq(vx, vy, vParams, is_fixed, Num_iter, fx_STEP);
		//LMfit(vx, vy, vParams, is_fixed, Num_iter, fx_STEP, 1/*boost*/);
	}
	else
		levmarq(vx, vy, vParams, is_fixed, Num_iter, fx_STEP);
	//LMfit(vx, vy, vParams, is_fixed, Num_iter, fx_STEP, 1/*boost*/);

	if (vParams[3] < 0.0
		|| vParams[3] > 100.0
		|| is_invalid(vParams[3]))
	{
		vector <myflo> vx, vy;

		vx.assign(vx.begin(), vx.begin() + vx.size() * 0.5);
		vy.assign(vy.begin(), vy.begin() + vx.size() * 0.5);

		vector <myflo> vAB = linear_fit(vx, vy);

		vParams[0] = vAB[0];
		vParams[1] = vAB[1];
		vParams[3] = 15;
		vParams[2] = abs((vy.back() - vParams[0] - vParams[1] * vx.back()) / (exp(vx.back() / vParams[3])));
		//vParams[2] = -vParams[0];

		is_fixed[0] = is_fixed[1] = true;

		vx.assign(vx.begin() + vx.size() * (1 - 0.5), vx.end());
		vy.assign(vy.begin() + vx.size() * (1 - 0.5), vy.end());

		levmarq(vx, vy, vParams, is_fixed, max(Num_iter, 400), fx_STEP);
		//LMfit(vx, vy, vParams, is_fixed, max(Num_iter, 400), fx_STEP, 1/*boost*/);

		if (vParams[3] < 0.0
			|| vParams[3] > 100.0
			|| is_invalid(vParams[3]))
		{
			vParams[3] = 15;
			vParams[2] = abs((vy.back() - vParams[0] - vParams[1] * vx.back()) / (exp(vx.back() / vParams[3])));
			//vParams[2] = - vAB[0] / 1E5;
		}
	}

	/* записываем в vres */
	for (size_t i = 0; i < vx.size(); ++i)
		vres[i] = fx_STEP(vx[i], vParams);

	/* записываем в vcoeffs */
	vcoeffs.resize(4);
	vcoeffs[0] = find_floating_potential(vx, vParams);							// Floating potential
	vcoeffs[1] = find_ion_current(vParams[0], vParams[1], vcoeffs[0]);				// Saturate ion current
	vcoeffs[2] = vParams[3];														// Temperature
	vcoeffs[3] = find_density(vcoeffs[1], vcoeffs[2]);								// Density ne


	return 0;
}

#endif