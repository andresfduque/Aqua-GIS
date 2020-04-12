#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 08:46:04 2016

PROCESS TIMESERIES:
    + Resample timeseries (on progress)
    + Convert from matrix to vector
    + Convert from vector to matrix

REQUIREMENTS:
    + Python 3.5

@author:    Andrés Felipe Duque Pérez
email:      andresfduque@gmail.com
"""


# %% Convert timeseries to matrix
def vectorToMatrix(timeseries):
    """Convert a vectorial timeseries into a matricial timeseries

    Timeseries: Vector[date, value]

    Matrix structure daily: [monthly date, day1, day2,......, day31]
    """
    import pandas as pd
    import numpy as np

    if timeseries.index.is_all_dates:
        if timeseries.index.freqstr == 'D':
            values = timeseries.values
            dates = timeseries.index
            nyears = dates[-1].year - dates[0].year + 1
            nmonths = nyears * 12
            matrix = np.zeros((nmonths, 33))    # matrix form of the ts
            matrix[:, :] = np.nan

            # locate each timeseries value in the matrix
            for i in range(0, len(values)):
                row = 12 * (dates[i].year - dates[0].year) + dates[i].month - 1
                col = dates[i].day + 1
                matrix[row, col] = values[i]

            # fill year and month
            month = 1
            year = dates[0].year
            for i in range(0, nmonths):
                matrix[i, 0] = year
                matrix[i, 1] = month
                month += 1
                if month == 13:
                    month = 1
                    year += 1

            df_columns = ['Year', 'Month', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6',
                          'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13', 'D14',
                          'D15', 'D16', 'D17', 'D18', 'D19', 'D20', 'D21',
                          'D22', 'D23', 'D24', 'D25', 'D26', 'D27', 'D28',
                          'D29', 'D30', 'D31']

            matrix = pd.DataFrame(matrix, columns=df_columns)
            return matrix

        elif timeseries.index.freqstr == 'M':

            values = timeseries.values
            dates = timeseries.index
            nyears = dates[-1].year - dates[0].year + 1
            nmonths = nyears * 12
            matrix = np.zeros((nyears, 13))     # matrix form of the ts
            matrix[:, :] = np.nan

            # locate each timeseries value in the matrix
            year = dates[0].year
            for i in range(0, nmonths):
                row = dates[i].year - dates[0].year
                col = dates[i].month
                matrix[row, col] = values[i]

            # years
            year = dates[0].year
            for i in range(0, nyears):
                matrix[i, 0] = year
                year += 1

            df_columns = ['Year', 'Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun',
                          'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic']

            matrix = pd.DataFrame(matrix, columns=df_columns)
            return matrix

    else:
        print('!!!Error: not a Pandas timeseries!!!')


# %% Convert matrix to timeseries
def matrixToVector(dataframe, freq='D'):
    """Convert a vectorial timeseries into a matricial timeseries

    dataframe: pandas dataframe (SIRH form)

    Matrix structure daily: [year, month, day1, day2,......, day31]
    """
    import pandas as pd
    import numpy as np

    if freq == 'D':
        # Create empty timeseries
        index = dataframe.index
        start_year = dataframe['Año'][index[0]]
        end_year = dataframe['Año'][index[-1]]
        start_date = pd.datetime(start_year, 1, 1)
        end_date = pd.datetime(end_year, 12, 31)
        date_range = pd.date_range(start_date, end_date, freq='D')
        series_values = np.zeros((len(date_range)))
        series_values[:] = np.nan
        ts = pd.Series(series_values, date_range)

        # Fill timeseries from dataframe
        nmonths = len(index)
        row = 0
        for i in range(nmonths):
            for day in range(1, 32):
                year = dataframe['Año'][index[row]]
                month = dataframe['Mes'][index[row]]
                value = dataframe['Dia ' + str(day)][index[row]]
                try:
                    date = pd.datetime(year, month, day)
                    if not np.isnan(value):
                        ts[date] = value
                except ValueError:
                    pass

            row += 1

        return ts


# %% Basic timeseries analysis
def basicTimeseriesProcessing(timeseries, variable=int, missingThreshold=float):
    """Basic timeseries processing\n
        - Timeseries: pandas timeseries [daily, monthly or annual]\n
        - Variable: 1 for precipitation, 2 for daily mean flows, 3 for monthly max instantly flows, 4 for monthly min mean flows\n
        - MissingThreshold: Threshold for missing data in resample process\n
        \n
        Functions:\n
        + Resample from daily to mothly timeseries (if daily timeseries)\n
        + Annual cycle\n
        + Missing analysis
    """

    import pandas as pd
    import numpy as np

    month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if isinstance(timeseries, pd.Series):
        # for daily or monthly data
        if timeseries.index.freqstr == 'D' or timeseries.index.freqstr == 'M':
            if timeseries.index.freqstr == 'D':

                daily_ts = pd.DataFrame(timeseries)

                # resample to monthly timeseries
                monthly_ts = timeseries.resample('M').mean()

                # monthly missing values
                timeseries_miss = timeseries.isnull()
                monthly_miss_ts_days = timeseries_miss.resample('M').sum()
                monthly_miss_ts_perc = monthly_miss_ts_days/monthly_miss_ts_days.index.day

                # correct monthly mean [depending on missing trheshold]
                monthly_ts[monthly_miss_ts_perc >= missingThreshold] = np.nan

                # put data in matrix form [pandas dataframe]
                matrix_daily_ts = vectorToMatrix(timeseries)
                matrix_missing_monthly_ts = vectorToMatrix(monthly_miss_ts_perc * 100)

            elif timeseries.index.freqstr == 'M':
                daily_ts = None
                monthly_ts = timeseries
                matrix_daily_ts = None
                matrix_missing_monthly_ts = None

            # annual cycle
            annual_cycle_mean = monthly_ts.groupby(lambda x: x.month).mean()

            if variable == 1:   # if data is cummulative [i.e. precipitation]
                annual_cycle_mean = pd.Series(annual_cycle_mean.values * month_days)
            else:
                annual_cycle_mean = monthly_ts.groupby(lambda x: x.month).mean()
                annual_cycle_std = monthly_ts.groupby(lambda x: x.month).std()
                annual_cycle_months = monthly_ts.groupby(lambda x: x.month).count()
                annual_cycle_error = annual_cycle_std/annual_cycle_months

            # annual timeseries (mean, max, min)
            annual_mean_ts = monthly_ts.resample('A').mean()
            annual_max_ts = monthly_ts.resample('A').max()
            annual_min_ts = monthly_ts.resample('A').min()

            annual_miss = monthly_ts.isnull()
            annual_miss_ts_months = annual_miss.resample('A').sum()
            annual_miss_ts_perc = annual_miss_ts_months / 12

            annual_mean_ts[annual_miss_ts_perc > missingThreshold] = np.nan
            annual_max_ts[annual_miss_ts_perc > missingThreshold] = np.nan
            annual_min_ts[annual_miss_ts_perc > missingThreshold] = np.nan

            annual_ts = pd.DataFrame({'Mean': annual_mean_ts.values, 'Max': annual_max_ts.values,
                                      'Min': annual_min_ts.values}, index=annual_mean_ts.index)

            # interanual mean
            if variable == 1:   # if data is cummulative [i.e. precipitation]
                monthly_ts = monthly_ts*monthly_miss_ts_days.index.day  # accumulate pcp
                inter_mean = annual_cycle_mean.sum()
                annual_cycle_std = monthly_ts.groupby(lambda x: x.month).std()
                annual_cycle_months = monthly_ts.groupby(lambda x: x.month).count()
                annual_cycle_error = annual_cycle_std/annual_cycle_months
            else:
                inter_mean = annual_cycle_mean.mean()

            # monthly timeseries to dataframe
            monthly_ts = pd.DataFrame(monthly_ts)

            # put data in matrix form [pandas dataframe]
            matrix_monthly_ts = vectorToMatrix(monthly_ts)
            annual_cycle = annual_cycle_mean.values.reshape((1, 12))
            annual_cycle_error = annual_cycle_error.values.reshape((1, 12))
            annual_cycle_columns = ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul',
                                    'Ago', 'Sep', 'Oct', 'Nov', 'Dic']
            matrix_annual_cycle = pd.DataFrame(annual_cycle, columns=annual_cycle_columns)
            matrix_annual_cycle_err = pd.DataFrame(annual_cycle_error, columns=annual_cycle_columns)

        elif timeseries.index.freqstr == 'A':
            daily_ts = None
            monthly_ts = None
            annual_ts = pd.DataFrame([timeseries, timeseries, timeseries],
                                     index=annual_mean_ts.index, columns=['Mean', 'Max', 'Min'])
            annual_cycle_mean = None
            matrix_daily_ts = None
            matrix_monthly_ts = None
            matrix_missing_monthly_ts = None
            matrix_annual_cycle = None
            matrix_annual_cycle_err = None

        return {'daily_ts': daily_ts,
                'monthly_ts': monthly_ts,
                'annual_ts': annual_ts,
                'matrix_daily_ts': matrix_daily_ts,
                'matrix_monthly_ts': matrix_monthly_ts,
                'matrix_missing_monthly_ts': matrix_missing_monthly_ts,
                'matrix_annual_cycle': matrix_annual_cycle,
                'matrix_annual_cycle_error': matrix_annual_cycle_err,
                'interanual_mean': inter_mean}


# %% Export daily data
def exportToExcel(dataDictionary, output_path, filename):
    """Export daily, monthly or annual timeseries basic processing to excel
    """

    import pandas as pd

    # dictionary variables
    inter_mean = dataDictionary['interanual_mean']
    daily_ts = dataDictionary['daily_ts']
    monthly_ts = dataDictionary['monthly_ts']
    annual_ts = dataDictionary['annual_ts']
    matrix_daily_ts = dataDictionary['matrix_daily_ts']
    matrix_monthly_ts = dataDictionary['matrix_monthly_ts']
    matrix_annual_cycle = dataDictionary['matrix_annual_cycle']
    matrix_annual_cycle_err = dataDictionary['matrix_annual_cycle_error']
    matrix_missing_monthly_ts = dataDictionary['matrix_missing_monthly_ts']

    # excel writer
    writer_path = (output_path + filename)
    writer = pd.ExcelWriter(writer_path, engine='xlsxwriter')
    workbook = writer.book

    # excel formats
    format1 = workbook.add_format({'num_format': 'dd/mm/yyyy', 'align': 'right',
                                   'valign': 'vcenter'})
    format2 = workbook.add_format({'num_format': '0', 'align': 'right',
                                   'valign': 'vcenter'})
    format3 = workbook.add_format({'num_format': '#,##0.00', 'align': 'right',
                                   'valign': 'vcenter'})
    format4 = workbook.add_format({'num_format': '#,##0.00', 'align': 'right',
                                   'valign': 'vcenter', 'bold': True})

    # write excel report [series]
    if daily_ts is not None:
        daily_ts.to_excel(writer, sheet_name='Series_Tiempo')
        monthly_ts.to_excel(writer, sheet_name='Series_Tiempo', startcol=2)
    else:
        if monthly_ts is not None:
            monthly_ts.to_excel(writer, sheet_name='Series_Tiempo')
        else:
            annual_ts.to_excel(writer, sheet_name='Series_Tiempo')

    worksheet = writer.sheets['Series_Tiempo']
    worksheet.set_column('A:A', 20, format1)
    worksheet.set_column('B:B', 18, format3)
    worksheet.set_column('C:C', 20, format1)
    worksheet.set_column('D:D', 18, format3)

    # write excel report [daily matrix]
    if daily_ts is not None:
        matrix_daily_ts.to_excel(writer, sheet_name='Matriz_Diaria',
                                 startcol=0, startrow=0, index=False)

        workbook = writer.book
        worksheet = writer.sheets['Matriz_Diaria']

        worksheet.set_column('A:B', 18, format2)
        worksheet.set_column('C:AH', 18, format3)

    # write excel report [monthly matrix]
    matrix_monthly_ts.to_excel(writer, sheet_name='Matriz_Mensual',
                               startcol=0, startrow=0, index=False)
    workbook = writer.book
    worksheet = writer.sheets['Matriz_Mensual']
    worksheet.set_column('A:A', 18, format2)
    worksheet.set_column('B:S', 18, format3)

    # write excel report [annual timeseries]
    annual_ts.to_excel(writer, sheet_name='Matriz_Mensual',
                       startcol=16, startrow=0, index=False)

    workbook = writer.book
    worksheet = writer.sheets['Matriz_Mensual']

    # Write excel report [annual cycle]
    matrix_annual_cycle.to_excel(writer, sheet_name='Matriz_Mensual', header=False, index=False,
                                 startcol=1, startrow=len(matrix_monthly_ts) + 1)
    matrix_annual_cycle_err.to_excel(writer, sheet_name='Matriz_Mensual', header=False, index=False,
                                     startcol=1, startrow=len(matrix_monthly_ts) + 2)

    worksheet.set_row(len(matrix_monthly_ts) + 1, None, format4)
    worksheet.set_row(len(matrix_monthly_ts) + 2, None, format4)

    # write interanual mean
    worksheet.write_string(len(matrix_monthly_ts), 13, 'mean')
    worksheet.write_number(len(matrix_monthly_ts) + 1, 13, inter_mean)

    # Write excel report [mothly missing matrix]
    matrix_missing_monthly_ts.to_excel(writer, sheet_name='Matriz_Mensual', index=False, startcol=0,
                                       startrow=len(matrix_monthly_ts) + 4)

    # Save excel file
    writer.save()
