#define FUNCTIONS_AND_CONSTANTS_REDEFINITIONS 0

#include "opencl.hpp" // for WORKGROUP_SIZE definition
#include "kernel.hpp" // note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(
#include "../headers/workspace.h"

static std::string get_expression_as_string_accounting_variables(std::string&& str_x);
struct name_replacement_position { std::string name; std::string replacement; size_t position; };
inline void quick_sort_for_name_replacement_position_struct(name_replacement_position* mas, size_t size);

std::string opencl_c_container()
{ 
	return R( // begin of OpenCL C code

			__kernel void Chi_sqr(__global const double* vx, __global const double* vy, const int size, __global const double* vparams, __global double* chi_sqr_total)
			{
				// l_id - position of work_item in current work_group
				// g_id - position of work_item in oveall problem
				int l_id = get_local_id(0), g_id = get_global_id(0);
				int group_size = get_local_size(0);

				double v1 = g_id * 2 < size ? ) + get_expression_as_string_accounting_variables("vx[g_id*2]") + R( - vy[g_id * 2] : 0;
				double v2 = g_id * 2 + 1 < size ? ) + get_expression_as_string_accounting_variables("vx[g_id*2+1]") + R( - vy[g_id * 2 + 1] : 0;

				// Local memory
				__local double partial_sums[) + to_string(WORKGROUP_SIZE) + R(];
				partial_sums[l_id] = v1 * v1 + v2 * v2;

				barrier(CLK_LOCAL_MEM_FENCE);

				for (int i = 2; i <= group_size; i <<= 1)
				{
					if (l_id % i == 0)
						partial_sums[l_id] += partial_sums[l_id + (i >> 1)];
					barrier(CLK_LOCAL_MEM_FENCE);
				}

				if (l_id == 0)
					chi_sqr_total[get_group_id(0)] = partial_sums[0];
			}

	) + R(	__kernel void multipliaction_arrarr(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] *= B[n];
			}

			__kernel void division_arrarr(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] /= B[n];
			}

			__kernel void subtraction_arrarr(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] -= B[n];
			}

			__kernel void addition_arrarr(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] += B[n];
			}

	) + R(	__kernel void multiplication_arrconst(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] *= B[0];
			}

			__kernel void division_arrconst(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] /= B[0];
			}

			__kernel void subtraction_arrconst(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] -= B[0];
			}

			__kernel void addition_arrconst(__global double* A, __global double* B)
			{
				int n = get_global_id(0);
				A[n] += B[0];
			}

	) + R(	__kernel void add_kernel(__global double* A, __global double* B, __global double* C)
			{
				int n = get_global_id(0);
				C[n] = A[n] + B[n];
			}
			__kernel void sub_kernel(__global double* A, __global double* B, __global double* C)
			{
				int n = get_global_id(0);
				C[n] = A[n] - B[n];
			}

			); // end of OpenCL C code
}

/* replaces all the variable in string expression with necessary gpu buffer pointer calls */
static std::string get_expression_as_string_accounting_variables(std::string && str_x)
{
	std::string str_return = workspace.s_expression;
	std::vector<name_replacement_position> replacements_and_places;

	/* searching cycle */
	for (size_t i = 0, pos = 0; i < workspace.names_of_expression_variables.size(); ++i, pos = 0)
	{
		std::string replacement = i == 0 ? str_x : "vparams[" + std::to_string(i - 1) + "]";

		while ((pos = str_return.find(workspace.names_of_expression_variables[i], pos)) != std::string::npos)
			replacements_and_places.push_back(name_replacement_position(workspace.names_of_expression_variables[i], replacement, pos)), pos += workspace.names_of_expression_variables[i].size();
	}

	/* position sorting */
	quick_sort_for_name_replacement_position_struct(replacements_and_places.data(), replacements_and_places.size());

	/* position correcting and replacement cycle */
	for (size_t i = 0, replacements_length = 0; i < replacements_and_places.size(); ++i)
		str_return.replace(replacements_and_places[i].position += replacements_length, replacements_and_places[i].name.size(), replacements_and_places[i].replacement),
		replacements_length += replacements_and_places[i].replacement.size() - 1;

	return str_return;
}

/* ascending quick sort for name_replacement_position structure, declared above */
inline void quick_sort_for_name_replacement_position_struct(name_replacement_position* mas, size_t size)
{
	//Указатели в начало и в конец массива
	int64_t i = 0;
	int64_t j = size - 1;

	//Центральный элемент массива
	size_t mid = mas[size / 2].position;

	//Делим массив
	do
	{
		//Пробегаем элементы, ищем те, которые нужно перекинуть в другую часть
		//В левой части массива пропускаем(оставляем на месте) элементы, которые меньше центрального
		while (mas[i].position < mid) i++;
		//В правой части пропускаем элементы, которые больше центрального
		while (mas[j].position > mid) j--;

		//Меняем элементы местами
		if (i <= j) std::swap(mas[i++], mas[j--]);
	} while (i <= j);


	//Рекурсивные вызовы, если осталось, что сортировать
	if (j > 0)
		//"Левый кусок"
		quick_sort_for_name_replacement_position_struct(mas, j + 1);
	if (i < size)
		//"Првый кусок"
		quick_sort_for_name_replacement_position_struct(&mas[i], size - i);
}