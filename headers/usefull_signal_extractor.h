#pragma once

#ifndef USEFULL_SIGNAL_EXTRACTOR_H
#define USEFULL_SIGNAL_EXTRACTOR_H

#include <vector>       // for vector
#include <algorithm>    // for min_element|max_element
#include "../peak_finder/peak_finder.h"

/****************** ������� ������� �� ���� ******************/
/*
����� ������ ������� ������ �������� ���������, ���� ������� ��������� ����� ���� (mnL_mxL_mnR_mxR[5]), ���� ������� �������� ����� ���� (mnL_mxL_mnR_mxR[4])
-- ��� ����� ����� ��������� ��������� �� ����� ��� ����� ��������� ��������� �������(����� ���)  --
&&
	����� ������ ������� && ����������� ����� ������ �������� ����������, �������� �������� ���� ����� (mnL_mxL_mnR_mxR[0]) ������������ �� ����������� (coefficient)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
	||
	����� ������ ������� && ����������� ����� ������ �������� ����������, ������� ��������� ���� ����� (mnL_mxL_mnR_mxR[1]) ������������ �� ����������� (coefficient)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
	||
	����� ������ ������� && ����������� ����� ������ �������� ����������, �������� �������� ���� ������ (mnL_mxL_mnR_mxR[2]) ������������ �� ����������� (coefficient)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
	||
	����� ������ ������� && ����������� ����� ������ �������� ����������, ������� ��������� ���� ������ (mnL_mxL_mnR_mxR[3]) ������������ �� ����������� (coefficient)
	-- ��� ����� ����� ����� ������ �� ���������, �������� ���� ��� ���, � ������ ������� ��������� --
	-- ��� ����� ������ ����� ��� ���������� ���������� �������� --
*/
/****************** ������� ������� �� ���� ******************/

static inline int find_signal_beginning(_In_ const std::vector <double>& signal,
										_In_ const std::vector <double>& mnL_mxL_mnR_mxR,
										/*
										������ mnL_mxL_mnR_mxR:
											-[0] ������� ���� �����
											-[1] �������� ���� ������
											-[2] ������� ���� ������
											-[3] �������� ���� ������
											-[4] ������� ����� ����
											-[5] �������� ����� ����
										*/
										_In_ const unsigned int myboost = 1)
{

#define I workspace.segments_beginning_indices
#define LEN workspace.one_segment_length
#define L workspace.left_cut_off
#define R workspace.right_cut_off

	int ignored_gap_beginning = I[0] + LEN * (1 - R),
		ignored_gap_ending = I[0] + LEN * (1 + L);

	for (int i = I[0] + LEN * L, k = 0; i < I.back(); i += myboost)
	{
		if (i < ignored_gap_beginning || i > ignored_gap_ending)
		{
			if (i > ignored_gap_ending)
			{
				++k;
				ignored_gap_beginning = I[k] + LEN * (1 - R);
				ignored_gap_ending = I[k] + LEN * (1 + L);
			}

			if ((signal[i] < mnL_mxL_mnR_mxR[4] || signal[i] > mnL_mxL_mnR_mxR[5])
				&& ((signal[i] < mnL_mxL_mnR_mxR[0] && signal[i + 1] < mnL_mxL_mnR_mxR[0])
					|| (signal[i] > mnL_mxL_mnR_mxR[1] && signal[i + 1] > mnL_mxL_mnR_mxR[1])
					|| (signal[i] < mnL_mxL_mnR_mxR[2] && signal[i + 1] < mnL_mxL_mnR_mxR[2])
					|| (signal[i] > mnL_mxL_mnR_mxR[3] && signal[i + 1] > mnL_mxL_mnR_mxR[3])))
				return i;
		}
	}

	return 0;

#undef I
#undef LEN
#undef L
#undef R

};

static inline int find_signal_ending(_In_ const std::vector <double>& signal,
									 _In_ const std::vector <double>& mnL_mxL_mnR_mxR,
									 /*
									������ mnL_mxL_mnR_mxR:
										-[0] ������� ���� �����
										-[1] �������� ���� ������
										-[2] ������� ���� ������
										-[3] �������� ���� ������
										-[4] ������� ����� ����
										-[5] �������� ����� ����
									*/
									 _In_ const unsigned int myboost = 1)
{

#define I workspace.segments_beginning_indices
#define LEN workspace.one_segment_length
#define L workspace.left_cut_off
#define R workspace.right_cut_off

	int ignored_gap_beginning = I.back() - LEN * (1 - R),
		ignored_gap_ending = I.back() - LEN * (1 + L);

	for (int i = I.back() - LEN * R, counter = 0, k = 0; i > I[0]; i -= myboost)
	{
		if (i > ignored_gap_beginning || i < ignored_gap_ending)
		{
			if (i < ignored_gap_ending)
			{
				++k;
				ignored_gap_beginning = I[I.size() - 1 - k] - LEN * (1 - R);
				ignored_gap_ending = I[I.size() - 1 - k] - LEN * (1 + L);
			}

			if ((signal[i] < mnL_mxL_mnR_mxR[4] || signal[i] > mnL_mxL_mnR_mxR[5])
				&& ((signal[i] < mnL_mxL_mnR_mxR[0] && signal[i - 1] < mnL_mxL_mnR_mxR[0])
					|| (signal[i] > mnL_mxL_mnR_mxR[1] && signal[i - 1] > mnL_mxL_mnR_mxR[1])
					|| (signal[i] < mnL_mxL_mnR_mxR[2] && signal[i - 1] < mnL_mxL_mnR_mxR[2])
					|| (signal[i] > mnL_mxL_mnR_mxR[3] && signal[i - 1] > mnL_mxL_mnR_mxR[3])))
				return i;
		}
	}

	return 0;

#undef I
#undef LEN
#undef L
#undef R

};

/* ��������� ������� (��������������������), �� ��� ��������. ����� ��� �� ����� �������,
�������� �� ������� ������ ��������(������������ ������� ����). ��� �� ���������� ������ �������
������ ������ �������� +-�����, ������� ���������� ������������(�������������� ��������� � workspace).
����� ����� ������������ � ����������� �������� ���� ����� � ������ � ���������� ��� ��������� ������ */
static inline std::vector<double> get_noise_amp(_In_ const std::vector <double>& signal,
												_In_opt_ const unsigned int k = 0,
												_In_opt_ const double coefficient = 1)
{

#define I workspace.segments_beginning_indices
#define LEN workspace.one_segment_length
#define L workspace.left_cut_off
#define R workspace.right_cut_off

	std::vector <double> left_noise(LEN),
		right_noise(LEN),
		mnL_mxL_mnR_mxR(6);	/*
							������ mnL_mxL_mnR_mxR:
								-[0] ������� ���� �����
								-[1] �������� ���� ������
								-[2] ������� ���� ������
								-[3] �������� ���� ������
								-[4] ������� ����� ����
								-[5] �������� ����� ����
							*/

		/* ������ ��� ���������� left_noise � right_noise */
	unsigned int i = 0,
		/* ������ ������ ��������, ������������ �����, � ����� ������� */
		l_ind = k,
		/* ������ ������ ��������, ������������ �����, � ����� ������� */
		r_ind = I.size() - 1 - k,
		/* "��������", ����������� ������ ����� �������� �� signal (��� left_noise) */
		l_it = I[l_ind] > LEN * (L + R) ?
		LEN * L : I[l_ind++] + LEN * L,
		/* "��������", ����������� ������ ����� �������� �� signal (��� right_noise) */
		r_it = right_noise.size() - I[r_ind] > LEN * (L + R) ?
		I[r_ind] + LEN * L : I[--r_ind] + LEN * L;

	/* ��������� ����� ������ ���� */
	for (; i < LEN && l_it < signal.size() * 0.5;
		l_it = I[l_ind++] + LEN * L)
	{
		unsigned int step_length = I[l_ind] - LEN * R - l_it;
		if (i + step_length > LEN)
			memcpy(&left_noise[i], &signal[l_it], LEN - i), i = LEN;
		else memcpy(&left_noise[i], &signal[l_it], step_length), i += step_length;
	}

	/* ��������� ������ ������ ���� */
	for (i = 0; i < LEN && r_it > signal.size() * 0.5;
		r_it = I[--r_ind] + LEN * L)
	{
		unsigned int step_length = I[r_ind] + LEN * (1 - R) - r_it;
		if (i + step_length > LEN)
			memcpy(&right_noise[i], &signal[r_it], LEN - i), i = LEN;
		else memcpy(&right_noise[i], &signal[r_it], step_length), i += step_length;
	}

	/* ��������� ������ mnL_mxL_mnR_mxR */
	double avg = (myspace::sum(left_noise) + myspace::sum(right_noise)) / (LEN * 2),
		mnL = *min_element(left_noise.begin(), left_noise.end()),
		mxL = *max_element(left_noise.begin(), left_noise.end()),
		mnR = *min_element(right_noise.begin(), right_noise.end()),
		mxR = *max_element(right_noise.begin(), right_noise.end());

	/* coefficient - ����������� ��� ����������� ������ ����,
	�� ����������� �� ����� ������ ������� �� ������� ����� (������� ���� �������) */
	if (avg != 0)
	{
		mnL_mxL_mnR_mxR[0] = avg - abs(avg - mnL) * coefficient;
		mnL_mxL_mnR_mxR[1] = avg + abs(mxL - avg) * coefficient;
		mnL_mxL_mnR_mxR[2] = avg - abs(avg - mnR) * coefficient;
		mnL_mxL_mnR_mxR[3] = avg + abs(mxR - avg) * coefficient;
	}
	else
	{
		if (mnL < 0) mnL_mxL_mnR_mxR[0] = mnL * coefficient;
		else mnL_mxL_mnR_mxR[0] = mnL / coefficient;
		if (mxL > 0) mnL_mxL_mnR_mxR[1] = mxL * coefficient;
		else mnL_mxL_mnR_mxR[1] = mxL / coefficient;
		if (mnR < 0) mnL_mxL_mnR_mxR[2] = mnR * coefficient;
		else mnL_mxL_mnR_mxR[2] = mnR / coefficient;
		if (mxR > 0) mnL_mxL_mnR_mxR[3] = mxR * coefficient;
		else mnL_mxL_mnR_mxR[3] = mxR / coefficient;
	}

	/* ���������� ������� ������� ���� */
	if (mxL > mxR) mnL_mxL_mnR_mxR[5] = avg + abs(mxL - avg);
	else mnL_mxL_mnR_mxR[5] = avg + abs(mxR - avg);
	if (mnL < mnR) mnL_mxL_mnR_mxR[4] = avg - abs(avg - mnL);
	else mnL_mxL_mnR_mxR[4] = avg - abs(avg - mnR);

	return mnL_mxL_mnR_mxR;

#undef I
#undef LEN
#undef L
#undef R

};

/* �������� �������� �������� ����� ���� � ������� �������(����������� ��� �������, �������� ������������� ���������� ��������) */
static inline int match_ramp_and_signal_periods(_In_ const unsigned int voltage_size,
												_In_ const unsigned int signal_size)
{

#define ISIGN workspace.segments_beginning_indices

	/* ���� ������� ������ � �������� � ����������� ��������� - �� ����� ����������� ������� ���� �� ������� */
	if (voltage_size == signal_size) return 0;

	/* ����� ��������� ����� ������� ������� int � ��������� ������, ��� ������� ���������� ���� */
	std::vector <double> holder(ISIGN.begin(), ISIGN.end());
	
	voltage_size > signal_size ? holder /= (double)voltage_size / signal_size : holder *= (double)signal_size / voltage_size;

	/* ���������� ������ � �������� ������ ���� int */
	ISIGN.assign(holder.begin(), holder.end());

	/*��������� ���������� �������� �������*/
	while (ISIGN.back() > signal_size)
		ISIGN.pop_back();

	return 0;

#undef ISIGN

};

/* �������������� ������ (signal � voltage) ����� �� workspace. �������� ������ ���� � ������ � workspace.ramp_peaks_indices.
��� �� ������� � ������� ��� �� ����� � ������� ���, ������ ������ ��������(� �������� ����) � �������� �������� � workspace.segments_beginning_indices */
inline int erase_noise()
{

#define IVOLT workspace.ramp_peaks_indices
#define ISIGN workspace.segments_beginning_indices
#define LEN workspace.one_segment_length
#define RFR workspace.ramp_freq		// ��� ���� � workspace ������ ���� ��������� �������
#define BTIME workspace.begin_time	// ��� ���� � workspace ������ ���� ��������� �������
#define ETIME workspace.end_time	// ��� ���� � workspace ������ ���� ��������� �������

	if (workspace.is_empty()) return ERR_StorageIsEmpty;

	int err = 0, begin, end;
	std::vector <double> voltage(workspace.voltage), signal(workspace.get_signal_ptr(), workspace.get_signal_ptr() + workspace.get_signal_length());

	/* ������� ��� ������� �� ���� ��� ������ ����� */
	err = PeakFinder::findPeaks(voltage, IVOLT), ISIGN = IVOLT;
	if (IVOLT.size() <= 3) // 3 ������ ��� ������ � ������� � ������� ������ ������� � �� ����� �������, ������ �� ������� �� �������� �������
		return ERR_TooFewSegs;

	/* ���������� ����� ������� ����������� ���� */
	int todel = (int)abs(IVOLT[2] - IVOLT[1]) % 100; double length_between_voltage_peaks;
	if ((double)todel / 100 <= 0.5) length_between_voltage_peaks = abs(IVOLT[2] - IVOLT[1]) - todel;
	else length_between_voltage_peaks = abs(IVOLT[2] - IVOLT[1]) + 100 - todel;

	/* ���������� ����� ������� ������� */
	LEN = (voltage.size() != signal.size()) ? length_between_voltage_peaks * ((double)signal.size() / voltage.size()) : length_between_voltage_peaks;

	/* �������� ������ */
	if (LEN == 0
		|| length_between_voltage_peaks == 0
		|| myspace::is_invalid(LEN)
		|| myspace::is_invalid(length_between_voltage_peaks)
		|| LEN > signal.size()
		|| length_between_voltage_peaks > signal.size())
		return ERR_BadSegsLength;

	/* ��������� ������� ����� �������� ������� �� ��������� ����� �������/���������� ����� �������� �������������� ������
	(�� ����� ��� ��������� ������� �������� ����) */
	match_ramp_and_signal_periods(voltage.size(), signal.size());

	/* ���� ������������ ��� ������ �������� ��������� ->
	���������� ������� �� ���� � ������ ���������� ������� �������(� ������ ������� ����), ��������� ������������� */
	if (BTIME != ETIME)
	{
		double total_time = (double)signal.size() / LEN / RFR;

		if (BTIME < 0
			|| BTIME > total_time
			|| BTIME > ETIME
			|| ETIME < 0
			|| ETIME > total_time)
			return ERR_BadStEndTime;

		begin = signal.size() * BTIME / total_time;
		end = signal.size() * BTIME / total_time;

		goto SkipWhile;
	}/*   ^_______________________   */
	 /*                          |   */
	{/* <-- ������ - ������� �� goto */
		double coefficient = 1.5; int segments = 0, k = 0;

		/* ���� ��������� ����, ����� ����� ����� �������� ������(�� ������ ���� ������� ���������, ��� ���) */
		std::vector <double> mnL_mxL_mnR_mxR;

		/* �������� ���� ������ ������� */
		do
		{
			/* ���� 5 ������� ����� ������ �� ���� */
			if (k > 5) return ERR_TooManyAttempts;

			/* ���� ��������� ����, ����� ����� ����� �������� ������(�� ������ ���� ������� ���������, ��� ���) */
			mnL_mxL_mnR_mxR = get_noise_amp(signal, k++, coefficient);

			/* ���� ������ � ����� ���� ��� ������ ������� � ���������� ����(mnL_mxL_mnR_mxR) */
			begin = find_signal_beginning(signal, mnL_mxL_mnR_mxR, (double)signal.size() / 1E5);
			end = find_signal_ending(signal, mnL_mxL_mnR_mxR, (double)signal.size() / 1E5);

			/* �������� ������ */
			if (begin < end || (begin != 0 && end != 0))
				segments = abs(end - begin) / LEN;
			else return ERR_BadStartEnd;

			/* ���� ���������� ���� �������� ������ 0.1 �� ������������ �������������, �� �� ������� ���� -> ������������� (����� ��������� ����������) */
			if (segments <= (int)(0.1f * signal.size() / LEN))
			{
				/* ��������� ������������ ����� �������� ������ ������ �������� */
				coefficient -= 0.3f;

				/* �������� ���-�� �������� ����� ����� ����� � ���� while */
				segments = 0;
			}
			/* ���� ���������� ���� �������� ������ 0.9 �� ������������ �������������, �� �� ������� ����� -> ������������� (����� ��������� ����������) */
			if (segments >= (int)(0.9f * signal.size() / LEN))
			{
				/* ����������� ������������ ����� �������� ������ ������ �������� */
				coefficient += 0.3f;

				/* �������� ���-�� �������� ����� ����� ����� � ���� while */
				segments = 0;
			}

		} while (segments == 0);
	}

SkipWhile:

	/* ���������� ������� ����� �������� � ������ ������ � ����� ������� */
	for (int i = 0; i < ISIGN.size(); ++i)
		if (ISIGN[i] > begin)
		{
			ISIGN.assign(ISIGN.begin() + i, ISIGN.begin() + i + (abs(end - begin) / LEN + i < ISIGN.size() ? abs(end - begin) / LEN : ISIGN.size() - i));
			break;
		}

	/* ������� ������ � ����� */
	for (int i = ISIGN.size() - 1; i > 1; --i)
		if (ISIGN[i] != 0)
		{
			ISIGN.resize(i + 1);
			break;
		}

	if (ISIGN.size() <= 0)
		return ERR_NoSegs;

	return err;

#undef IVOLT
#undef ISIGN
#undef LEN
#undef RFR
#undef BTIME
#undef ETIME

};

#endif
