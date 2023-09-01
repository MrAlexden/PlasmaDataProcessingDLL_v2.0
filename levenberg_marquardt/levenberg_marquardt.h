#pragma once

#ifndef LEVENBERG_MARQURDT_H
#define LEVENBERG_MARQURDT_H

#if GPU_COMPUTATION > 0 || STRING_EXPRESSION > 0

#include "../myspace/mynamespace.h"

namespace mystringcompute
{

#define VECTOR std::vector<T>

	template <typename T>
		requires std::is_convertible_v<T, double>
	inline int LevenbergMarquardt(_In_ const VECTOR& vx,								// independent data
								  _In_ const VECTOR& vy,								// dependent data
								  _Inout_ VECTOR& vp,									// vector of parameters we are searching for
#if STRING_EXPRESSION <= 0
								  _In_ T(*fx)(const T, const VECTOR&),					// the function that the data will be approximated by
#endif
								  _In_opt_ const std::vector <bool>& vF = {},			// vector of param's fixed status 
								  _In_opt_ const unsigned int niter = 400,				// number of max iterations this method does
								  _In_opt_ const unsigned int myboost = 1)				// speed up calculation coefficient
	{
		int x, i, j, it, npar = 0; std::vector <bool> vFixed(vF);
		T lambda = 0.01f, up = 10, down = 1 / up, mult, err = 0, newerr = 0, derr = 0, target_derr = 1E-12;

		/* check for input mistakes */
		if (vx.size() != vy.size()
			|| vx.size() <= myboost
			|| vy.size() <= myboost)
			return -1;
		if (vFixed.empty())
			for (i = 0; i < vp.size(); ++i)
				vFixed.push_back(false);
		if (vp.size() != vFixed.size())
			return -2;
		if (niter == NULL)
			return -3;
#if STRING_EXPRESSION > 0
		if (workspace.variables_values.size() - 1 != vp.size())
			return -4;
#endif

		/* params counting, fixed params are not calculated */
		for (i = 0; i < vFixed.size(); ++i)
		{
			if (!vFixed[i])
				++npar;
			if (myspace::is_invalid(vp[i]))
				vp[i] = 1;
		}

		myspace::matrix <T> Hessian(npar, npar, 0.0), ch(npar, npar, 0.0);
		VECTOR grad(npar), drvtv(npar), delta(npar, 0.0), newvP = vp;

		/* calculate the initial error ("chi-squared") */
#if GPU_COMPUTATION > 0
		mygpucompute::create_gpu_kernel_for_Chi_sqr(vy.size(), workspace.names_of_expression_variables.size() - 1);

		mygpucompute::set_gpu_independent_data(vx.data());
		mygpucompute::set_gpu_dependent_data(vy.data());

		mygpucompute::set_gpu_variables_data(vp.data());

		err = mygpucompute::calculate_gpu_Chi_sqr();
#elif STRING_EXPRESSION > 0
		err = s_Chi_sqr(vx, vy, vp, myboost);
#else
		err = myspace::Chi_sqr(vx, vy, vp, fx, myboost);
#endif

		/* main iteration */
		for (it = 0; it < niter; ++it)
		{
			/* zeroing */
			for (i = 0; i < npar; ++i)
			{
				drvtv[i] = 0.0;
				for (j = 0; j < npar; ++j)
					Hessian[i][j] = 0;
			}

#if STRING_EXPRESSION > 0
			if (workspace.variables_values.size() > 1) memcpy(&workspace.variables_values[1], vp.data(), sizeof(T) * vp.size());
#endif

			/* calculate the approximation to the Hessian and the "derivative" drvtv */
			for (x = 0 + (myboost - 1); x < vy.size() - (myboost - 1); x += myboost)
			{
#if STRING_EXPRESSION > 0
				workspace.variables_values[0] = vx[x];  T f = expression_parser.Eval();
#endif

				/* calculate gradient */
				for (i = 0, j = 0; i < vFixed.size(); ++i, ++j)
					if (!vFixed[i])
#if STRING_EXPRESSION > 0
	#define derivative_step 0.001
					{
						workspace.variables_values[i + 1] += derivative_step;
						grad[j] = (expression_parser.Eval() - f) / derivative_step;
						workspace.variables_values[i + 1] = vp[i];
					}
	#undef derivative_step
#else
						grad[j] = partial_derivative(vx[x], vp, fx, i);
#endif

				for (i = 0; i < npar; ++i)
				{
#if STRING_EXPRESSION > 0
					drvtv[i] += (vy[x] - f) * grad[i] * myboost;
#else
					drvtv[i] += (vy[x] - fx(vx[x], vp)) * grad[i] * myboost;
#endif

					for (j = 0; j < npar; ++j)
						Hessian[i][j] += grad[i] * grad[j] * myboost;
				}
			}

			/*  make a step "delta."  If the step is rejected, increase lambda and try again */
			mult = 1 + lambda;
			bool ill = true; /* ill-conditioned? */
			while (ill && (it < niter))
			{
				for (i = 0; i < npar; ++i)
					Hessian[i][i] = Hessian[i][i] * mult;

				ill = myspace::cholesky_decomp(ch, Hessian);

				if (!ill)
				{
					myspace::solve_axb_cholesky(ch, delta, drvtv);

					for (i = 0, j = 0; i < vFixed.size(); ++i)
						if (!vFixed[i])
							newvP[i] = vp[i] + delta[j++];

#if GPU_COMPUTATION > 0
					mygpucompute::set_gpu_variables_data(newvP.data());
					newerr = mygpucompute::calculate_gpu_Chi_sqr();
#elif STRING_EXPRESSION > 0
					newerr = s_Chi_sqr(vx, vy, newvP, myboost);
#else
					newerr = myspace::Chi_sqr(vx, vy, newvP, fx, myboost);
#endif
					derr = newerr - err;
					ill = (derr > 0);
				}

				if (ill)
				{
					mult = (1 + lambda * up) / (1 + lambda);
					lambda *= up;
					++it;
				}
			}

			for (i = 0; i < vFixed.size(); ++i)
				if (!vFixed[i])
					vp[i] = newvP[i];

			err = newerr;
			lambda *= down;

			if ((!ill) && (-derr < target_derr))
				//if ((!ill) && (abs(err) < target_derr))
				break;
		}

		return 0;
	};

#undef VECTOR

}

#endif /* GPU_COMPUTATION > 0 || STRING_EXPRESSION > 0 */

#endif /* LEVENBERG_MARQURDT_H */