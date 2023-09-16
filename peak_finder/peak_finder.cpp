#include "peak_finder.h"

static inline void vectorElementsProduct(_In_ const std::vector <double>& a,
										 _In_ const std::vector <double>& b,
										 _Out_ std::vector <double>& out)
{
	out.resize(a.size());

	for (int i = 0; i < a.size(); ++i)
		out[i] = a[i] * b[i];
}

static inline void findIndicesLessThan(_In_ const std::vector <double>& in,
									   _In_ const double threshold,
									   _Out_ std::vector <long long>& indices)
{
	for (int i = 0; i < in.size(); ++i)
		if (in[i] < threshold)
			indices.push_back(i + 1);
}

template <typename T>
	requires std::is_convertible_v<T, double>
static inline void selectElementsFromIndices(_In_ const std::vector <T>& in,
											 _In_ const std::vector <long long>& indices,
											 _Out_ std::vector <T>& out)
{
	out.resize(indices.size());

	for (int i = 0; i < indices.size(); ++i)
		out[i] = in[indices[i]];
}

static inline void signVector(const std::vector <double>& in, std::vector<long long>& out)
{
	out.resize(in.size());

	for (int i = 0; i < in.size(); ++i)
	{
		if (in[i] > 0)
			out[i] = 1;
		else if (in[i] < 0)
			out[i] = -1;
		else
			out[i] = 0;
	}
}

int PeakFinder::findPeaks(_In_ const std::vector <double>& in,
						  _Out_ std::vector <long long>& peakInds,
						  _In_opt_ bool includeEndpoints,
						  _In_opt_ int extrema)
{
	if (includeEndpoints == NULL) includeEndpoints = false;
	if (extrema == NULL) extrema = 1;
	std::vector <double> x0 = in;

	int minIdx = distance(x0.begin(), min_element(x0.begin(), x0.end()));
	int maxIdx = distance(x0.begin(), max_element(x0.begin(), x0.end()));

	double sel = (x0[maxIdx] - x0[minIdx]) / 4.0f;
	int len0 = x0.size(), err = 0;

	if (extrema == -1) x0 *= extrema;

	std::vector <double> dx = diff(x0);
	replace(dx.begin(), dx.end(), (double)0.0, -PeakFinder::EPS);
	std::vector <double> dx0(dx.begin(), dx.end() - 1);
	std::vector <double> dx0_1(dx.begin() + 1, dx.end());
	std::vector <double> dx0_2;

	vectorElementsProduct(dx0, dx0_1, dx0_2);

	std::vector<long long> ind;
	findIndicesLessThan(dx0_2, 0, ind); // Find where the derivative changes sign	
	std::vector <double> x;
	double leftMin;
	int minMagIdx;
	double minMag;

	if (includeEndpoints)
	{
		//x = [x0(1);x0(ind);x0(end)];	
		selectElementsFromIndices(x0, ind, x);
		x.insert(x.begin(), x0[0]);
		x.insert(x.end(), x0.back());
		//ind = [1;ind;len0];
		ind.insert(ind.begin(), 1);
		ind.insert(ind.end(), len0);
		minMagIdx = distance(x.begin(), min_element(x.begin(), x.end()));
		minMag = x[minMagIdx];
		//std::cout<<"Hola"<<std::endl;
		leftMin = minMag;
	}
	else
	{
		selectElementsFromIndices(x0, ind, x);
		if (x.size() > 2)
		{
			minMagIdx = distance(x.begin(), min_element(x.begin(), x.end()));
			minMag = x[minMagIdx];
			leftMin = x[0] < x0[0] ? x[0] : x0[0];
		}
	}

	int len = x.size();

	if (len > 2)
	{
		double tempMag = minMag;
		bool foundPeak = false;
		int ii;

		if (includeEndpoints)
		{
			// Deal with first point a little differently since tacked it on
			// Calculate the sign of the derivative since we tacked the first
			//  point on it does not neccessarily alternate like the rest.
			std::vector <double> xSub0(x.begin(), x.begin() + 3); //tener cuidado subvector
			std::vector <double> xDiff = diff(xSub0); //tener cuidado subvector

			std::vector<long long> signDx;
			signVector(xDiff, signDx);

			if (signDx[0] <= 0) // The first point is larger or equal to the second
			{
				if (signDx[0] == signDx[0]) // Want alternating signs
				{
					x.erase(x.begin() + 1);
					ind.erase(ind.begin() + 1);
					len = len - 1;
				}
			}
			else // First point is smaller than the second
			{
				if (signDx[0] == signDx[0]) // Want alternating signs
				{
					x.erase(x.begin());
					ind.erase(ind.begin());
					len = len - 1;
				}
			}
		}

		//Skip the first point if it is smaller so we always start on maxima
		if (x[0] >= x[1])
			ii = 0;
		else
			ii = 1;

		//Preallocate max number of maxima
		double maxPeaks = ceil((double)len / 2.0f);
		std::vector<int> peakLoc(maxPeaks, 0);
		std::vector <double> peakMag(maxPeaks, 0.0);
		int cInd = 1;
		int tempLoc;

		while (ii < len)
		{
			ii = ii + 1;//This is a peak
			//Reset peak finding if we had a peak and the next peak is bigger
			//than the last or the left min was small enough to reset.
			if (foundPeak)
			{
				tempMag = minMag;
				foundPeak = false;
			}

			//Found new peak that was lager than temp mag and selectivity larger
			//than the minimum to its left.

			if (x[ii - 1] > tempMag && x[ii - 1] > leftMin + sel)
			{
				tempLoc = ii - 1;
				tempMag = x[ii - 1];
			}

			//Make sure we don't iterate past the length of our std::vector
			if (ii == len)
				break; //We assign the last point differently out of the loop

			ii = ii + 1; // Move onto the valley

			//Come down at least sel from peak
			if (!foundPeak && tempMag > sel + x[ii - 1])
			{
				foundPeak = true; //We have found a peak
				leftMin = x[ii - 1];
				peakLoc[cInd - 1] = tempLoc; // Add peak to index
				peakMag[cInd - 1] = tempMag;
				cInd = cInd + 1;
			}
			else if (x[ii - 1] < leftMin) // New left minima
				leftMin = x[ii - 1];

		}

		// Check end point
		if (includeEndpoints)
		{
			if (x.back() > tempMag && x.back() > leftMin + sel)
			{
				peakLoc[cInd - 1] = len - 1;
				peakMag[cInd - 1] = x.back();
				cInd = cInd + 1;
			}
			else if (!foundPeak && tempMag > minMag)// Check if we still need to add the last point
			{
				peakLoc[cInd - 1] = tempLoc;
				peakMag[cInd - 1] = tempMag;
				cInd = cInd + 1;
			}
		}
		else if (!foundPeak)
		{
			double minAux = x0.back() < x.back() ? x0.back() : x.back();
			if (x.back() > tempMag && x.back() > leftMin + sel)
			{
				peakLoc[cInd - 1] = len - 1;
				peakMag[cInd - 1] = x.back();
				cInd = cInd + 1;
			}
			else if (!tempMag > minAux + sel)// Check if we still need to add the last point
			{
				peakLoc[cInd - 1] = tempLoc;
				peakMag[cInd - 1] = tempMag;
				cInd = cInd + 1;
			}
		}

		//Create output
		if (cInd > 0)
		{
			try
			{
				if (cInd < peakLoc.size())
					selectElementsFromIndices(ind, std::vector<long long>(peakLoc.begin(), peakLoc.begin() + cInd - 1), peakInds);
				else
					selectElementsFromIndices(ind, std::vector<long long>(peakLoc.begin(), peakLoc.begin() + cInd - 1), peakInds);
			}
			catch (...)
			{
				err = ERR_Exception;
			}
		}
	}

	return err;
}