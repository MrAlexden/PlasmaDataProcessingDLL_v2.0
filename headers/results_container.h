#pragma once

#ifndef PROC_RES_H
#define PROC_RES_H

#include <vector>       // for vector

template <typename T>
	requires std::is_convertible_v<T, double>
class Results
{
private:

	std::vector <T> _ramp;							//	�������������� ���� (1 ������ (��� �������)). ������ ����� ���� ������� 1000 - 4000. ������ ���� �������� 1 / (50 ��) ���
	std::vector <myspace::matrix <T>> _results;		/*	��� ������������� ���������� �������������� ��������������.������ _results ����� ���� �������� 3 - 4.
														� ������ i-�� �������� _results �������� ������� � �������������� ������������ ������������
														�������������� ��������, ��������:
														0 - ������������ ������ ������� �������, "����������� �����";
														1 - ���������� ������� �������;
														2 - ����������� ���������� ������� �������;
														3 - ������������ ���������� ������� ������� � �.�.
														������ ���� ������ �������� � ������ ��������� 90 - 100 (��������(������)) �� 1000 - 4000 (����� � �������(�������)).
														������� �������� ������ � ������ _ramp ������ ���������. */
	myspace::matrix <T> _scalings;					/*	������� �� ����� ����������� ������� 90 - 100 (��������(������)) �� 4 - 5 (����������(�������)).
														��� ������� ������� � ���� ������� �������� ���������� ���������� ����������,
														����������� ��� ���� � ���������� ��������� (������ �������� ������ ����� ������ �������). ��������:
														1 - ������������ ������ � ���� �������;
														2 - ������� ������� ������ � ���� �������;
														3 - ��������� ���� ������������������� � ���� �������;
														4 - ��������� ������ � ���� �������;
														���������� ����� ������� _scalings ������ ��������� � ����������� ����� ������ -results[i]. */

	size_t _number_of_results;		// ���������� ����������� (��� ������������ ������ �� �������) ����������� ������ ������������� �������������� �������� (3 - 4 �����)					
	size_t _number_of_segments;		// ���������� ������������ �������� (�������� ����������� ������������ ����������) �� ����� ������� (90 - 100 ����)
	size_t _size_of_segment;		// ��������� ����� � ����� ������� (�����, ������ �� ������ ���� (1000 - 4000 ����)). ������� ������������� ��� �������� ������� 150 - 200 ���.
	size_t _number_of_scalings;		// ���������� ����������� ������� (��� �� ��������� �� ������� (������ �������)) (4 - 5 ����)

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