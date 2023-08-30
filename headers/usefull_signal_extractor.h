#pragma once

#ifndef USEFULL_SIGNAL_EXTRACTOR_H
#define USEFULL_SIGNAL_EXTRACTOR_H

#include <vector>       // for vector
#include <algorithm>    // for min_element|max_element
#include "../peak_finder/peak_finder.h"

/****************** УСЛОВИЕ ОЧИСТКИ ОТ ШУМА ******************/
/*
Точка начала сигнала должна обладать значением, либо большим максимума всего шума (mnL_mxL_mnR_mxR[5]), либо меньшим минимума всего шума (mnL_mxL_mnR_mxR[4])
-- это нужно чтобы исключить вхождение по массе при очень маленькой амплитуде сигнала(почти шум)  --
&&
	Точка начала сигнала && последующая точка должны обладать значениями, меньшими минимума шума слева (mnL_mxL_mnR_mxR[0]) домноженного на коэффициент (coefficient)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
	||
	Точка начала сигнала && последующая точка должны обладать значениями, большим максимума шума слева (mnL_mxL_mnR_mxR[1]) домноженного на коэффициент (coefficient)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
	||
	Точка начала сигнала && последующая точка должны обладать значениями, меньшими минимума шума справа (mnL_mxL_mnR_mxR[2]) домноженного на коэффициент (coefficient)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
	||
	Точка начала сигнала && последующая точка должны обладать значениями, большим максимума шума справа (mnL_mxL_mnR_mxR[3]) домноженного на коэффициент (coefficient)
	-- это может сразу найти сигнал по амплитуде, например если шум мал, а сигнал высокой амплитуды --
	-- две точки подряд взяты для исключения случайного всплеска --
*/
/****************** УСЛОВИЕ ОЧИСТКИ ОТ ШУМА ******************/

static inline int find_signal_beginning(_In_ const std::vector <double>& signal,
										_In_ const std::vector <double>& mnL_mxL_mnR_mxR,
										/*
										вектор mnL_mxL_mnR_mxR:
											-[0] минимум шума слева
											-[1] максимум шума справа
											-[2] минимум шума справа
											-[3] максимум шума справа
											-[4] минимум всего тока
											-[5] максимум всего тока
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
									вектор mnL_mxL_mnR_mxR:
										-[0] минимум шума слева
										-[1] максимум шума справа
										-[2] минимум шума справа
										-[3] максимум шума справа
										-[4] минимум всего тока
										-[5] максимум всего тока
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

/* дебильная функция (суперНЕуниверсальная), но что поделать. Берет шум по краям сигнала,
опираясь на индексы начала отрезков(образованных попилом пилы). Так же игнорирует учетки сигнала
вокруг границ отрезков +-части, которые отсекаются пользователе(предварительно загнанные в workspace).
Далее берет максимальные и минимальные значения шума слева и справа и возвращяет как результат работы */
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
							вектор mnL_mxL_mnR_mxR:
								-[0] минимум шума слева
								-[1] максимум шума справа
								-[2] минимум шума справа
								-[3] максимум шума справа
								-[4] минимум всего тока
								-[5] максимум всего тока
							*/

		/* индекс для заполнения left_noise и right_noise */
	unsigned int i = 0,
		/* индекс начала сегмента, распиленного пилой, с левой стороны */
		l_ind = k,
		/* индекс начала сегмента, распиленного пилой, с левой стороны */
		r_ind = I.size() - 1 - k,
		/* "итератор", указывающий откуда берем значения из signal (для left_noise) */
		l_it = I[l_ind] > LEN * (L + R) ?
		LEN * L : I[l_ind++] + LEN * L,
		/* "итератор", указывающий откуда берем значения из signal (для right_noise) */
		r_it = right_noise.size() - I[r_ind] > LEN * (L + R) ?
		I[r_ind] + LEN * L : I[--r_ind] + LEN * L;

	/* заполняем левый вектор шума */
	for (; i < LEN && l_it < signal.size() * 0.5;
		l_it = I[l_ind++] + LEN * L)
	{
		unsigned int step_length = I[l_ind] - LEN * R - l_it;
		if (i + step_length > LEN)
			memcpy(&left_noise[i], &signal[l_it], LEN - i), i = LEN;
		else memcpy(&left_noise[i], &signal[l_it], step_length), i += step_length;
	}

	/* заполняем правый вектор шума */
	for (i = 0; i < LEN && r_it > signal.size() * 0.5;
		r_it = I[--r_ind] + LEN * L)
	{
		unsigned int step_length = I[r_ind] + LEN * (1 - R) - r_it;
		if (i + step_length > LEN)
			memcpy(&right_noise[i], &signal[r_it], LEN - i), i = LEN;
		else memcpy(&right_noise[i], &signal[r_it], step_length), i += step_length;
	}

	/* заполняем вектор mnL_mxL_mnR_mxR */
	double avg = (myspace::sum(left_noise) + myspace::sum(right_noise)) / (LEN * 2),
		mnL = *min_element(left_noise.begin(), left_noise.end()),
		mxL = *max_element(left_noise.begin(), left_noise.end()),
		mnR = *min_element(right_noise.begin(), right_noise.end()),
		mxR = *max_element(right_noise.begin(), right_noise.end());

	/* coefficient - коэффициент для регулировки высоты шума,
	он подгоняется во время поиска сигнала во внешнем цикле (снаружи этой функции) */
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

	/* определяем средний уровень шума */
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

/* приводит значения индексов пиков пилы к размеру сигнала(растягивает или сжимает, оставляя относительное расстояния отрезков) */
static inline int match_ramp_and_signal_periods(_In_ const unsigned int voltage_size,
												_In_ const unsigned int signal_size)
{

#define ISIGN workspace.segments_beginning_indices

	/* если размеры данных с сигналом и напряжением совпадают - не нужно выравнивать индексы пилы по сигналу */
	if (voltage_size == signal_size) return 0;

	/* нужен держатель чтобы сделать индексы int с плавающей точкой, для точного приведения фазы */
	std::vector <double> holder(ISIGN.begin(), ISIGN.end());
	
	voltage_size > signal_size ? holder /= (double)voltage_size / signal_size : holder *= (double)signal_size / voltage_size;

	/* возвращаем данные в основной вектор типа int */
	ISIGN.assign(holder.begin(), holder.end());

	/*проверяем превышение крайнего индекса*/
	while (ISIGN.back() > signal_size)
		ISIGN.pop_back();

	return 0;

#undef ISIGN

};

/* первоначальные данные (signal и voltage) берет из workspace. Выделяет период пилы и кладет в workspace.ramp_peaks_indices.
Так же находит в сигнале шум по краям и убирает его, кладет начала отрезков(с периодом пилы) с полезным сигналом в workspace.segments_beginning_indices */
inline int erase_noise()
{

#define IVOLT workspace.ramp_peaks_indices
#define ISIGN workspace.segments_beginning_indices
#define LEN workspace.one_segment_length
#define RFR workspace.ramp_freq		// это поле в workspace должно быть заполнено заранее
#define BTIME workspace.begin_time	// это поле в workspace должно быть заполнено заранее
#define ETIME workspace.end_time	// это поле в workspace должно быть заполнено заранее

	if (workspace.is_empty()) return ERR_StorageIsEmpty;

	int err = 0, begin, end;
	std::vector <double> voltage(workspace.voltage), signal(workspace.get_signal_ptr(), workspace.get_signal_ptr() + workspace.get_signal_length());

	/* находим все участки на пиле при помощи пиков */
	err = PeakFinder::findPeaks(voltage, IVOLT), ISIGN = IVOLT;
	if (IVOLT.size() <= 3) // 3 потому что дальше в функцию я передаю начало второго и он нужен целиком, тоесть от второго до третьего индекса
		return ERR_TooFewSegs;

	/* определяем длину участка изначальной пилы */
	int todel = (int)abs(IVOLT[2] - IVOLT[1]) % 100; double length_between_voltage_peaks;
	if ((double)todel / 100 <= 0.5) length_between_voltage_peaks = abs(IVOLT[2] - IVOLT[1]) - todel;
	else length_between_voltage_peaks = abs(IVOLT[2] - IVOLT[1]) + 100 - todel;

	/* определяем длину участка сигнала */
	LEN = (voltage.size() != signal.size()) ? length_between_voltage_peaks * ((double)signal.size() / voltage.size()) : length_between_voltage_peaks;

	/* проверка ошибок */
	if (LEN == 0
		|| length_between_voltage_peaks == 0
		|| myspace::is_invalid(LEN)
		|| myspace::is_invalid(length_between_voltage_peaks)
		|| LEN > signal.size()
		|| length_between_voltage_peaks > signal.size())
		return ERR_BadSegsLength;

	/* домножаем индексы начал отрезков сигнала на отношение длины сигнала/напряжения чтобы получить действительные начала
	(до этого они равняются началам отрезков пилы) */
	match_ramp_and_signal_periods(voltage.size(), signal.size());

	/* если пользователь сам выбрал диапазон обработки ->
	пропускаем очистку от шума и просто возвращаем участок сигнала(с учетом периода пилы), выбранный пользователем */
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
	{/* <-- скобки - протект от goto */
		double coefficient = 1.5; int segments = 0, k = 0;

		/* ищем амплитуду шума, чтобы затем найти полезный сигнал(он должен быть большей амплитуды, чем шум) */
		std::vector <double> mnL_mxL_mnR_mxR;

		/* основной цикл поиска сигнала */
		do
		{
			/* есть 5 попыток найти сигнал из шума */
			if (k > 5) return ERR_TooManyAttempts;

			/* ищем амплитуду шума, чтобы затем найти полезный сигнал(он должен быть большей амплитуды, чем шум) */
			mnL_mxL_mnR_mxR = get_noise_amp(signal, k++, coefficient);

			/* ищем начало и конец шума при помощи вектора с амплитудой шума(mnL_mxL_mnR_mxR) */
			begin = find_signal_beginning(signal, mnL_mxL_mnR_mxR, (double)signal.size() / 1E5);
			end = find_signal_ending(signal, mnL_mxL_mnR_mxR, (double)signal.size() / 1E5);

			/* проверка ошибок */
			if (begin < end || (begin != 0 && end != 0))
				segments = abs(end - begin) / LEN;
			else return ERR_BadStartEnd;

			/* если количество всех отрезков меньше 0.1 от потенциально максимального, то их слишком мало -> пересчитываем (нужно увеличить количество) */
			if (segments <= (int)(0.1f * signal.size() / LEN))
			{
				/* уменьшаем коэффициенты чтобы понизить пороги отбора отрезков */
				coefficient -= 0.3f;

				/* зануляем кол-во отрезков чтобы снова зайти в цикл while */
				segments = 0;
			}
			/* если количество всех отрезков больше 0.9 от потенциально максимального, то их слишком много -> пересчитываем (нужно уменьшить количество) */
			if (segments >= (int)(0.9f * signal.size() / LEN))
			{
				/* увеличиваем коэффициенты чтобы повысить пороги отбора отрезков */
				coefficient += 0.3f;

				/* зануляем кол-во отрезков чтобы снова зайти в цикл while */
				segments = 0;
			}

		} while (segments == 0);
	}

SkipWhile:

	/* записываем индексы начал отрезков с учетом начала и конца сигнала */
	for (int i = 0; i < ISIGN.size(); ++i)
		if (ISIGN[i] > begin)
		{
			ISIGN.assign(ISIGN.begin() + i, ISIGN.begin() + i + (abs(end - begin) / LEN + i < ISIGN.size() ? abs(end - begin) / LEN : ISIGN.size() - i));
			break;
		}

	/* убираем лишнее с конца */
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
