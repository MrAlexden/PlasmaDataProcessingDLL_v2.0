#pragma once

#ifndef GPU_COMPUTE_H
#define GPU_COMPUTE_H

#if GPU_COMPUTATION > 0

#include "opencl.hpp"

namespace mygpucompute
{
	static Device device;		// compile OpenCL C code for the fastest available device

	static Memory<double> vx_global;
	static Memory<double> vy_global;
	static Memory<double> vp_global;
	static Memory<double> kernel_output;

	static std::map<std::string, Kernel> processes;		// kernel that runs on the device

	static std::string current;

	inline int create_gpu_kernel_for_Chi_sqr(const uint size_of_data, const uint num_of_variables)
	{
		if (!device.is_initialized())
			device.construct(select_device_with_most_flops());

		/* allocate memory on both hostand device */
		if (vx_global.length() != size_of_data)									vx_global = Memory<double>(device, size_of_data);
		if (vy_global.length() != size_of_data)									vy_global = Memory<double>(device, size_of_data);
		if (vp_global.length() != num_of_variables)								vp_global = Memory<double>(device, num_of_variables);
		if (kernel_output.length() != std::ceil(size_of_data / WORKGROUP_SIZE)) kernel_output = Memory<double>(device, std::ceil(size_of_data / WORKGROUP_SIZE));

		current = "Chi_sqr";

		/* фактически, второй параметр здесь это количество итераций цикла for */
		processes[current].construct(device, size_of_data, current, vx_global, vy_global, size_of_data, vp_global, kernel_output);

		return 0;
	}

	inline int set_gpu_independent_data(const double* vx)
	{
		std::memcpy(vx_global.data(), vx, vx_global.capacity());
		vx_global.write_to_device();
		return 0;
	}

	inline int set_gpu_dependent_data(const double* vy)
	{
		std::memcpy(vy_global.data(), vy, vy_global.capacity());
		vy_global.write_to_device();
		return 0;
	}

	inline int set_gpu_variables_data(const double* vp)
	{
		std::memcpy(vp_global.data(), vp, vp_global.capacity());
		vp_global.write_to_device();
		return 0;
	}

	inline double calculate_gpu_Chi_sqr()
	{
		processes[current].run();

		kernel_output.read_from_device();

		double result = 0;
		for (int i = 0; i < kernel_output.length(); ++i)
			result += kernel_output.data()[i];

		return result;
	}

	inline int clear_gpu()
	{
		processes.clear();

		vx_global.~Memory();
		vy_global.~Memory();
		vp_global.~Memory();
		kernel_output.~Memory();

		return 0;
	}

	/* sets up computing device and script by passed name of operation (example: "multiplication", "division", "subtraction", "addition").
	operation between two arrays will be performed */
	inline int create_gpu_kernel_for_arrays_ariphmetics(const uint size_of_data, const std::string operation)
	{
		if (!device.is_initialized())
			device.construct(select_device_with_most_flops());

		/* allocate memory on both hostand device */
		if (vx_global.length() != size_of_data) vx_global = Memory<double>(device, size_of_data);
		if (vy_global.length() != size_of_data) vy_global = Memory<double>(device, size_of_data);

		current = operation + "_arrarr";

		/* фактически, второй параметр здесь это количество итераций цикла for */
		processes[current].construct(device, size_of_data, current, vx_global, vy_global);

		return 0;
	}

	/* sets up computing device and script by passed name of operation (example: "multipliction", "division", "subtraction", "addition").
	operation between an array and a constant will be performed */
	inline int create_gpu_kernel_for_array_constant_ariphmetics(const uint size_of_data, const std::string operation)
	{
		if (!device.is_initialized())
			device.construct(select_device_with_most_flops());

		/* allocate memory on both hostand device */
		if (vx_global.length() != size_of_data) vx_global = Memory<double>(device, size_of_data);
		if (vy_global.length() != 1)			vy_global = Memory<double>(device, 1);

		current = operation + "_arrconst";

		/* фактически, второй параметр здесь это количество итераций цикла for */
		processes[current].construct(device, size_of_data, current, vx_global, vy_global);

		return 0;
	}

	inline double* calculate_gpu_ariphmetics()
	{
		processes[current].run();

		vx_global.read_from_device();

		return vx_global.data();
	}

}

#endif /* GPU_COMPUTATION > 0 */

#endif /* GPU_COMPUTE_H */
