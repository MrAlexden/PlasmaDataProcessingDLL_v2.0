#pragma once

#ifndef PROC_RES_H
#define PROC_RES_H

#include <vector>       // for vector

template <typename T>
	requires std::is_convertible_v<T, double>
class Results
{
public:

	Results() :
		NumberOfSegments(0),
		SizeOfSegment(0),
		NumberOfParameters(0)
	{
	}
	~Results()
	{
		ramp.clear();

		for (int i = 0; i < NumberOfSegments; ++i)
		{
			mOriginalData[i].clear();
			mFiltratedData[i].clear();
			mApproximatedData[i].clear();
			mDifferentiatedData[i].clear();
			mParametersData[i].clear();
		}

		mOriginalData.clear();
		mFiltratedData.clear();
		mApproximatedData.clear();
		mDifferentiatedData.clear();
		mParametersData.clear();
	}
	Results(int Number_Of_Segments, int Size_Of_Segment, int Number_Of_Parameters)
	{
		NumberOfSegments = Number_Of_Segments;
		SizeOfSegment = Size_Of_Segment;
		NumberOfParameters = Number_Of_Parameters;

		mOriginalData.resize(Number_Of_Segments);
		mFiltratedData.resize(Number_Of_Segments);
		mApproximatedData.resize(Number_Of_Segments);
		mDifferentiatedData.resize(Number_Of_Segments);
		mParametersData.resize(Number_Of_Segments);
	}
	void set_SegmentsNumber(int Number_Of_Segments)
	{
		NumberOfSegments = Number_Of_Segments;

		mOriginalData.resize(Number_Of_Segments);
		mFiltratedData.resize(Number_Of_Segments);
		mApproximatedData.resize(Number_Of_Segments);
		mDifferentiatedData.resize(Number_Of_Segments);
		mParametersData.resize(Number_Of_Segments);
	}
	void set_SegmentsSize(int Size_Of_Segment)
	{
		SizeOfSegment = Size_Of_Segment;
	}
	void set_ParamsNumber(int Number_Of_Parameters)
	{
		NumberOfParameters = Number_Of_Parameters;
	}


	/* DATA SETTING */
	void set_ramp(std::vector <T> & v)
	{
		ramp = v;
	}
	int set_OriginSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mOriginalData[i] = v;
	}
	int set_FiltedSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mFiltratedData[i] = v;
	}
	int set_ApproxSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mApproximatedData[i] = v;
	}
	int set_DiffedSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mDifferentiatedData[i] = v;
	}
	int set_ParamsSegment(std::vector <T> & v, int i)
	{
		if (i > NumberOfSegments || i < 0) return -1;
		mParametersData[i] = v;
	}


	/* DATA GETTING */
	std::vector <T> & const get_ramp()
	{
		return ramp;
	}
	myspace::matrix <T> & const get_mOriginalData() const 
	{
		return mOriginalData;
	}
	myspace::matrix <T> & const get_mFiltratedData() const
	{
		return mFiltratedData;
	}
	myspace::matrix <T> & const get_mApproximatedData() const
	{
		return mApproximatedData;
	}
	myspace::matrix <T> & const get_mDifferentiatedData() const
	{
		return mDifferentiatedData;
	}
	myspace::matrix <T> & const get_mParametersData() const
	{
		return mParametersData;
	}


	/* DIMENSION GETTING */
	int get_NumberOfSegments() const
	{
		return NumberOfSegments;
	}
	int get_SizeOfSegment() const
	{
		return SizeOfSegment;
	}
	int get_NumberOfParameters() const
	{
		return NumberOfParameters;
	}

private:

	std::vector <T> ramp;

	myspace::matrix <T> mOriginalData;
	myspace::matrix <T> mFiltratedData;
	myspace::matrix <T> mApproximatedData;
	myspace::matrix <T> mDifferentiatedData;
	myspace::matrix <T> mParametersData;

	int NumberOfSegments;
	int SizeOfSegment;
	int NumberOfParameters;
};

extern Results<double> results;

#endif