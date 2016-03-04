/*
 * Copyright 2015 Paul Chote
 * This file is part of tsspectrogram, which is free software. It is made available
 * to you under the terms of version 3 or later of the GNU General Public License,
 * as published by the Free Software Foundation. For more information, see LICENSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cpgplot.h>
#include <stdbool.h>
#include <stdarg.h>

// Prints an vararg error to stderr then returns 1
int error(const char * format, ...)
{
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n");
    return 1;
}

#define error_jump(label, ret, ...) do { ret = error(__VA_ARGS__); goto label; } while(0)

struct ts_data
{
    double *time;
    double *mmi;
    double *err;
    size_t obs_count;
};

// ".ts" data format is three space delimited columns:
// time (bjd) intensity (mmi) sigma (mmi)
static int load_tsfile(const char *ts_path, struct ts_data *data)
{
    int ret = 0;
    char linebuf[1024];

    FILE *file = fopen(ts_path, "r+");
    if (!file)
        error_jump(file_error, ret, "Unable to load %s", ts_path);

    size_t total = 0;
    while (fgets(linebuf, sizeof(linebuf) - 1, file))
        if (linebuf[0] != '#' && linebuf[0] != '\n')
            total++;
    rewind(file);

    data->time = calloc(total, sizeof(double));
    data->mmi = calloc(total, sizeof(double));
    data->err = calloc(total, sizeof(double));
    if (!data->time || !data->mmi || !data->err)
        error_jump(allocation_error, ret, "Allocation error");

    size_t count = 0;
    while (fgets(linebuf, sizeof(linebuf) - 1, file) && count < total)
    {
        // Skip comment / empty lines
        if (linebuf[0] == '#' || linebuf[0] == '\n')
            continue;

        int read = sscanf(linebuf, "%lf %lf %lf\n", &data->time[count], &data->mmi[count], &data->err[count]);

        // No error defined; set to unity
        if (read == 2)
            data->err[count] = 1;

        // Convert to seconds
        data->time[count] *= 86400;
        count++;
    }
    data->obs_count = count;
    fclose(file);

    return ret;

allocation_error:
    free(data->time); data->time = NULL;
    free(data->mmi); data->mmi = NULL;
    free(data->err); data->err = NULL;

file_error:
    data->obs_count = 0;
    return ret;
}

static void ts_data_free(struct ts_data *data)
{
    free(data->time);
    free(data->mmi);
    free(data->err);
    free(data);
}

static void set_color_table()
{
    float l[9] = {0.0, 0.33, 0.66, 1.0};
    float r[9] = {1.0, 0.0, 0.0, 1.0};
    float g[9] = {1.0, 0.0, 1.0, 0.0};
    float b[9] = {1.0, 1.0, 0.0, 0.0};
    cpgctab(l, r, g, b, 4, 1.0, 0.5);

    //cpgsitf(2); // Use a sqrt mapping between value and colour
}

int generate_spectrogram_image(float *image, double window_width,
    double *time, double *mmi, size_t count,
    double time_min, double time_max, size_t time_steps,
    double freq_min, double freq_max, size_t freq_steps)
{
    int ret = 0;
    double dt = (time_max - time_min) / time_steps;
    double df = (freq_max - freq_min) / freq_steps;

    double *temp_time = calloc(count, sizeof(double));
    if (!temp_time)
        error_jump(temp_time_alloc_error, ret, "Error allocating temp_time");

    double *temp_mmi = calloc(count, sizeof(double));
    if (!temp_mmi)
        error_jump(temp_mmi_alloc_error, ret, "Error allocating temp_mmi");

    for (size_t i = 0; i < time_steps; i++)
    {
        // Sample the windowed data in this slice
        double mid_time = time_min + i*dt;
        size_t n = 0;
        for (size_t j = 0; j < count; j++)
        {
            if ((time[j] >= mid_time - window_width / 2) && (time[j] <= mid_time + window_width / 2))
            {
                temp_time[n] = time[j];
                temp_mmi[n] = mmi[j];
                n++;
            }
        }

        // TODO: Filter condition

        // Calculate the DFT amplitudes for this time slice
        for (size_t j = 0; j < freq_steps; j++)
        {
            double real = 0;
            double imag = 0;
            double freq = freq_min + j*df;

            for (size_t k = 0; k < n; k++)
            {
                double phase = -freq*2*M_PI*(temp_time[k] - temp_time[0]);
                real += temp_mmi[k] * cos(phase) / n;
                imag += temp_mmi[k] * sin(phase) / n;
            }

            image[j*time_steps + i] = 2*sqrt(real*real + imag*imag);
        }
    }

    free(temp_mmi);
temp_mmi_alloc_error:
    free(temp_time);
temp_time_alloc_error:
    return ret;
}

int colorplot(const char *ts_path, double time_window_width, double min_amplitude, double max_amplitude, double freq_min, double freq_max, size_t freq_steps, size_t time_steps, const char *output_ps)
{
    int ret = 0;

    struct ts_data *data = calloc(1, sizeof(struct ts_data));
    if (!data)
        return 1;

    if (load_tsfile(ts_path, data))
        error_jump(load_failed_error, ret, "Error loading timeseries data");

    double time_min = data->time[0];
    double time_max = data->time[data->obs_count - 1];

    printf("Generating spectrogram of %s:\n", ts_path);
    printf("   Using %g second spectrogram window.\n", time_window_width);
    printf("   Using %zu frequency steps between %g and %g uHz.\n", freq_steps, freq_min, freq_max);
    printf("   Using %zu time steps between %g and %g seconds.\n", time_steps, time_min, time_max);
    if (output_ps == NULL)
        printf("   Displaying plot in X Window 3.\n\n");
    else
        printf("   Saving plot as %s.\n\n", output_ps);

    printf("Loaded %zu data points...\n", data->obs_count);

    // Plot layout
    const float plot_left = 0.15;
    const float plot_right = 0.85;
    const float plot_ts_top = 0.225;
    const float plot_ts_bottom = 0.075;
    const float plot_data_bottom = 0.075+0.18;
    const float plot_data_top = 0.5925+0.18;
    const float plot_window_bottom = 0.6+0.18;
    const float plot_window_top = 0.7125+0.18;
    const float plot_scale_top = 0.92;
    const float plot_scale_bottom = 0.90;
    const size_t plot_scale_steps = 50;
    const float plot_label_margin = 4;

    float plot_window_data_ratio = (plot_window_top - plot_window_bottom) / (plot_data_top - plot_data_bottom);

    char device[128];
    if (output_ps != NULL)
        snprintf(device, 128, "%s/cps", output_ps);
    else
        strcpy(device, "3/xs");

    if (cpgopen(device) <= 0)
        error_jump(pgplot_open_error, ret, "Unable to open PGPLOT window");

    cpgpap(9.5, 0.8);
    cpgask(0);
    cpgslw(1);
    cpgsfs(2);
    cpgscf(2);
    cpgsch(1.0);

    // The same buffer used for the data, window, and scale images.
    // The data image is always the largest
    float *image = calloc(time_steps * freq_steps, sizeof(float));
    if (!image)
        error_jump(image_alloc_error, ret, "Error allocating image buffer");

    set_color_table();

    // Time series panel
    {
        float *temp_time = calloc(data->obs_count, sizeof(float));
        float *temp_mmi = calloc(data->obs_count, sizeof(float));
        double min_mmi = data->mmi[0];
        double max_mmi = data->mmi[0];

        for (size_t i = 0; i < data->obs_count; i++)
        {
            temp_time[i] = data->time[i];
            temp_mmi[i] = data->mmi[i];
            min_mmi = fmin(min_mmi, data->mmi[i]);
            max_mmi = fmax(max_mmi, data->mmi[i]);
        }

        double mid_mmi = (min_mmi + max_mmi) / 2;
        double range_mmi = (max_mmi - min_mmi) / 2;
        printf("%f %f\n", mid_mmi, range_mmi);
        min_mmi = mid_mmi - 1.2 * range_mmi;
        max_mmi = mid_mmi + 1.2 * range_mmi;

        cpgsvp(plot_left, plot_right, plot_ts_bottom, plot_ts_top);
        cpgswin(time_min, time_max, min_mmi, max_mmi);
        cpgsci(4);
        cpgpt(data->obs_count, temp_time, temp_mmi, 229);
        cpgsci(1);
        cpgswin(time_min / 3600, time_max / 3600, min_mmi, max_mmi);
        cpgbox("bcstn", 0, 0, "bcstnv", 0, 0);
        cpgmtxt("l", plot_label_margin, 0.5, 0.5, "(mmi)");
        cpgmtxt("B", 2.3, 0.5, 0.5, "Time (hours)");
    }

    // Data panel
    {
        printf("Generating spectrogram...\n");
        generate_spectrogram_image(image, time_window_width,
            data->time, data->mmi, data->obs_count,
            time_min, time_max, time_steps,
            freq_min * 1e-6, freq_max * 1e-6, freq_steps);

        cpgsvp(plot_left, plot_right, plot_data_bottom, plot_data_top);
        cpgmtxt("l", plot_label_margin, 0.5, 0.5, "Frequency (\\gmHz)");
        cpgswin(time_min, time_max, freq_min, freq_max);

        float x_scale = (time_max - time_min) / (time_steps - 1);
        float y_scale = (freq_max - freq_min) / (freq_steps - 1);
        float tr[] = { time_min - x_scale, x_scale, 0, freq_min - y_scale, 0, y_scale };

        cpgimag(image, time_steps, freq_steps, 1, time_steps, 1, freq_steps, min_amplitude, max_amplitude, tr);
        cpgswin(time_min / 86400, time_max / 86400, freq_min, freq_max);
        cpgbox("bcst", 0, 0, "bstnv", 0, 0);
        cpgbox("0", 0, 0, "c", 0, 0);

        double period_ticks[] = {110, 120, 135, 150, 170, 200, 250, 350, 500, 1000, 10000};
        size_t period_tick_count = sizeof(period_ticks) / sizeof(double);
        cpgswin(0, 1, freq_min, freq_max);
        for (size_t i = 0; i < period_tick_count; i++)
        {
            double period = period_ticks[i];
            double freq = 1e6/period;
            if (freq < freq_min || freq > freq_max)
                continue;

            char label[128];
            snprintf(label, 127, "%.0f", period);
            cpgmove(0.9925, freq);
            cpgdraw(1, freq);
            cpgmtxt("RV", 0.5, (freq - freq_min) / (freq_max - freq_min), 0, label);
        }

        cpgmtxt("r", 4, 0.5, 0.5, "Period (s)");
    }

    // Window
    {
        printf("Generating window spectrogram...\n");
        double dft_window_half_width = (freq_max - freq_min) * plot_window_data_ratio / 2;
        double dft_window_freq = (freq_max + freq_min) / 2;
        double window_freq_min = dft_window_freq - dft_window_half_width;
        double window_freq_max = dft_window_freq + dft_window_half_width;
        size_t window_freq_steps = (size_t)(freq_steps * plot_window_data_ratio);

        // Reuse the data intensity array for calculating the window
        for (size_t j = 0; j < data->obs_count; j++)
            data->mmi[j] = max_amplitude * sin(2 * M_PI * dft_window_freq * 1e-6 * data->time[j]);

        generate_spectrogram_image(image, time_window_width,
            data->time, data->mmi, data->obs_count,
            time_min, time_max, time_steps,
            window_freq_min * 1e-6, window_freq_max * 1e-6, window_freq_steps);

        cpgsvp(plot_left, plot_right, plot_window_bottom, plot_window_top);
        cpgmtxt("l", plot_label_margin, 0.5, 0.5, "Window");
        cpgswin(time_min, time_max, window_freq_min, window_freq_max);

        float x_scale = (time_max - time_min) / (time_steps - 1);
        float y_scale = 2 * dft_window_half_width / (window_freq_steps - 1);
        float tr[] = { time_min - x_scale, x_scale, 0, window_freq_min - y_scale, 0, y_scale };

        cpgimag(image, time_steps, window_freq_steps, 1, time_steps, 1, window_freq_steps, min_amplitude, max_amplitude, tr);

        cpgswin(time_min / 86400, time_max / 86400, -dft_window_half_width, dft_window_half_width);
        cpgbox("bcst", 0, 0, "bcstnv", 0, 0);
    }

    // Scale
    {
        cpgsvp(plot_left, plot_right, plot_scale_bottom, plot_scale_top);

        // Reuse t
        float ampl_scale = (max_amplitude - min_amplitude) / (plot_scale_steps - 1);
        for (size_t i = 0; i < plot_scale_steps; i++)
            image[i] = min_amplitude + i * ampl_scale;

        float tr[] = { min_amplitude - ampl_scale, 0, ampl_scale, -0.5, 1, 0, 0, 1 };
        cpgswin(min_amplitude, max_amplitude, 0, 1);
        cpgimag(image, 1, plot_scale_steps, 1, 1, 1, plot_scale_steps, min_amplitude, max_amplitude, tr);

        cpgbox("bcsm", 0, 0, "bc", 0, 0);
        cpgmtxt("t", 2, 0.5, 0.5, "Amplitude (mma)");
    }

    // Window size indicators
    {
        double half_dft = 1.0 / (time_steps - 1) / 2;
        double half_mmi = time_window_width / (time_max - time_min) / 2;
        cpgsvp(plot_left, plot_right, 0, 1);
        cpgswin(0, 1, 0, 1);
        cpgsci(2);
        double indicator_offset = 0.01;
        cpgmove(0.5 - half_mmi, plot_ts_top - indicator_offset);
        cpgdraw(0.5 - half_mmi, plot_ts_top + indicator_offset);
        cpgdraw(0.5 - half_dft, plot_data_bottom - indicator_offset);
        cpgdraw(0.5 - half_dft, plot_data_bottom + indicator_offset);
        cpgmove(0.5 + half_mmi, plot_ts_top - indicator_offset);
        cpgdraw(0.5 + half_mmi, plot_ts_top + indicator_offset);
        cpgdraw(0.5 + half_dft, plot_data_bottom - indicator_offset);
        cpgdraw(0.5 + half_dft, plot_data_bottom + indicator_offset);
    }

    printf("Done.\n");

    free(image);
image_alloc_error:
    cpgend();
pgplot_open_error:
    ts_data_free(data);
load_failed_error:
    return ret;
}

int main(int argc, char *argv[])
{
    if (argc >= 9)
    {
        const char *tsfile = argv[1];
        double window_width = atof(argv[2]);
        double ampl_min = atof(argv[3]);
        double ampl_max = atof(argv[4]);
        double freq_min = atof(argv[5]);
        double freq_max = atof(argv[6]);
        double freq_steps = atoi(argv[7]);
        double time_steps = atoi(argv[8]);
        const char *psfile = argc == 9 ? NULL : argv[9];
        return colorplot(tsfile, window_width, ampl_min, ampl_max,
            freq_min, freq_max, freq_steps, time_steps, psfile);
    }
    else
    {
        printf("tsspectrogram tsfile window_width min_amplitude max_amplitude freq_min freq_max freq_steps time_steps [output_ps]\n\n");
        printf("Arguments:\n");
        printf("    tsfile: timeseries data file.  Space delimited with time in days and amplitude in mmi.\n");
        printf("    window_width: Width of the spectrogram window in seconds.\n");
        printf("    ampl_min: Minimum DFT amplitude to display in mma.\n");
        printf("    ampl_max: Maximum DFT amplitude to display in mma.\n");
        printf("    freq_min: Minimum frequency to display in uHz.\n");
        printf("    freq_max: Maximum frequency to display in uHz.\n");
        printf("    freq_steps: Number of frequency samples.\n");
        printf("    time_steps: Number of time samples.\n");
        printf("    output_ps (optional): filename for output postscript plot.\n");
        return 1;
    }

    return 0;
}
