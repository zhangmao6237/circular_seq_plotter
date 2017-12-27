# coding=utf-8

import xlrd


def read_excel(path, title_row_num=0, col_range=None):
    wb = xlrd.open_workbook(path)
    ws = wb.sheets()[0]

    row_range = ws.nrows
    row_range -= title_row_num
    _get_cell = lambda row, col: ws.cell(row, col).value

    if col_range is None:
        col_range = ws.ncols

    # -------- row & col --------
    info = [[_get_cell(row + title_row_num, col) for col in xrange(col_range)] for row in xrange(row_range)]

    return info
