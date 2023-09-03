#pragma once

#ifndef PROC_RES_H
#define PROC_RES_H

#include <vector>       // for vector

template <typename T>
	requires std::is_convertible_v<T, double>
class Results
{
private:

	std::vector <T> _ramp;							//	Результирующая пила (1 период (или отрезок)). Размер может быть порядка 1000 - 4000. Период пилы например 1 / (50 Гц) сек
	std::vector <myspace::matrix <T>> _results;		/*	Все промежуточные результаты математических преобразований.размер _results может быть например 3 - 4.
														В каждом i-ом элементе _results хранится матрица с промежуточными вычислениями определенных
														математических операций, например:
														0 - оригинальный сигнал каждого отрезка, "порубленный пилой";
														1 - фильтрация каждого отрезка;
														2 - производная фильтрации каждого отрезка;
														3 - аппроксимаия фильтрации каждого отрезка и т.д.
														Размер всех матриц одинаков и обычно сотавляет 90 - 100 (отрезков(строки)) на 1000 - 4000 (точек в отрезке(столбцы)).
														Размеры столбцов матриц и размер _ramp всегда совпадают. */
	myspace::matrix <T> _scalings;					/*	Матрица со всеми скейлингами размера 90 - 100 (отрезков(строки)) на 4 - 5 (параметров(столбцы)).
														Для каждого отрезка в этой матрице хранится одинаковое количество параметров,
														расчитанных для него в результате обработки (первый параметр всегда время начала отрезка). Например:
														1 - максимальный сигнал в этом отрезке;
														2 - средняя энергия частиц в этом отрезке;
														3 - потенциал пика энергораспределения в этом отрезке;
														4 - плотность плазмы в этом отрезке;
														Количество строк матрицы _scalings всегда совпадает с количеством строк матриц -results[i]. */

	size_t _number_of_results;		// Количество сохраненных (для последующего вывода на графики) результатов работы промежуточных математических действий (3 - 4 штуки)					
	size_t _number_of_segments;		// Количество обнаруженных отрезков (периодов прохождения пилобразного напряжения) за время разряда (90 - 100 штук)
	size_t _size_of_segment;		// Количесто точек в одном отрезке (точки, снятые за период пилы (1000 - 4000 штук)). Частота дискретизации АЦП например порядка 150 - 200 кГц.
	size_t _number_of_scalings;		// Количество результатов расчета (они же скейлинги от времени (первый столбец)) (4 - 5 штук)

public:

	Results() :
		_number_of_results(0),
		_number_of_segments(0),
		_size_of_segment(0),
		_number_of_scalings(0)
	{}
	Results(const size_t number_of_results, const size_t number_of_segments, const size_t size_of_segment, const size_t number_of_scalings) :
		_number_of_results(number_of_results),
		_number_of_segments(number_of_segments),
		_size_of_segment(size_of_segment),
		_number_of_scalings(number_of_scalings)
	{
		_ramp.resize(size_of_segment);

		_results.resize(number_of_results);
		for (size_t i = 0; i < number_of_results; ++i)
			_results[i](number_of_segments, size_of_segment);

		_scalings(number_of_scalings, number_of_segments);
	}

	~Results()
	{
		_ramp.clear();
		_results.clear();
	}



	/* RESIZING */
	void set_number_of_results(const size_t number_of_results)
	{
		_number_of_results = number_of_results;

		_results.resize(number_of_results);
		if (_number_of_segments != 0 && _size_of_segment != 0)
			for (size_t i = 0; i < number_of_results; ++i)
				if (_results[i].empty())
					_results[i](_number_of_segments, _size_of_segment);
	}
	void set_number_of_segments(const size_t number_of_segments)
	{
		_number_of_segments = number_of_segments;

		if (_number_of_results != 0 && _size_of_segment != 0)
			for (size_t i = 0; i < _number_of_results; ++i)
					_results[i].set_size(number_of_segments, _size_of_segment);

		if (_number_of_scalings != 0)
			_scalings.set_size(_number_of_scalings, number_of_segments);
	}
	void set_size_of_segment(const size_t size_of_segment)
	{
		_size_of_segment = size_of_segment;

		_ramp.resize(size_of_segment);

		if (_number_of_results != 0 && _number_of_segments != 0)
			for (size_t i = 0; i < _number_of_results; ++i)
				_results[i].set_size(_number_of_segments, size_of_segment);
	}
	void set_number_of_scalings(const size_t number_of_scalings)
	{
		_number_of_scalings = number_of_scalings;

		if (_number_of_segments != 0)
			_scalings.set_size(number_of_scalings, _number_of_segments);
	}



	/* GETTING ACCESS */
	std::vector <T> & const ramp() const
	{
		//static_assert(/*!_ramp.empty()*//*fs != NULL*/, "You must set size_of_segment first");
		return const_cast<std::vector <T> &>(_ramp);
	}
	myspace::matrix <T> & const operator [] (size_t i) const
	{
		return _results[i];
	}
	myspace::matrix <T> & const scalings() const
	{
		return const_cast<std::vector <T> &>(scalings);
	}



	/* DIMENSION GETTING */
	size_t get_number_of_results() const
	{
		return _number_of_results;
	}
	size_t get_number_of_segments() const
	{
		return _number_of_segments;
	}
	size_t get_size_of_segment() const
	{
		return _size_of_segment;
	}
	size_t get_number_of_scalings() const
	{
		return _number_of_scalings;
	}
};

extern Results<double> results;

class asd
{
public:
private:
};

#endif