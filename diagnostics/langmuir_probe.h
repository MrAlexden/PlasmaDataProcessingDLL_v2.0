#pragma once

#ifdef LANGMUIR_PROBE
#define LANGMUIR_PROBE

#include <vector>
#include "../headers/error_codes.h"
#include "../myspace/mynamespace.h"

inline int langmuir_probe(std::vector<double> & vx, std::vector<double> & vy, )
{

	if (vx.size() != vy.size()
		|| vx.empty()
		|| vy.empty())
	{
		vcoeffs.resize(4);
		vcoeffs[0] = vcoeffs[1] = vcoeffs[2] = vcoeffs[3] = 0;
		return ERR_BadSegInput;
	}

	int i = 0;
	vector <bool> vFixed = { false, false, false, false };
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

		vFixed[0] = vFixed[1] = true;
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

		levmarq(vx, vy, vParams, vFixed, Num_iter, fx_STEP);
		//LMfit(vx, vy, vParams, vFixed, Num_iter, fx_STEP, 1/*boost*/);
	}
	else
		levmarq(vx, vy, vParams, vFixed, Num_iter, fx_STEP);
	//LMfit(vx, vy, vParams, vFixed, Num_iter, fx_STEP, 1/*boost*/);

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

		vFixed[0] = vFixed[1] = true;

		vx.assign(vx.begin() + vx.size() * (1 - 0.5), vx.end());
		vy.assign(vy.begin() + vx.size() * (1 - 0.5), vy.end());

		levmarq(vx, vy, vParams, vFixed, max(Num_iter, 400), fx_STEP);
		//LMfit(vx, vy, vParams, vFixed, max(Num_iter, 400), fx_STEP, 1/*boost*/);

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
	for (i = 0; i < vx.size(); ++i)
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