#!/usr/bin/gnuplot
#
#  Vector Network Analyzer Calibration Library
#  Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Plot script for example.c
#
set key left top
set title 'S Parameter Calibration'
set xlabel 'frequency (Hz)'
unset ylabel
set logscale x
plot 'vnacal-example.out' \
    index 0 using 1:2 title 'expected   S11_r' lt 1 with points, \
    '' index 1 using 1:2 title 'measured   S11_r' lt 1 dt 2 with lines, \
    '' index 2 using 1:2 title 'corrected  S11_r' lt 1 dt solid with lines, \
    '' index 0 using 1:3 title 'expected   S11_i' lt 2 with points, \
    '' index 1 using 1:3 title 'measured   S11_i' lt 2 dt 2 with lines, \
    '' index 2 using 1:3 title 'corrected  S11_i' lt 2 dt solid with lines, \
    '' index 0 using 1:4 title 'expected   S12_r' lt 3 with points, \
    '' index 1 using 1:4 title 'measured   S12_r' lt 3 dt 2 with lines, \
    '' index 2 using 1:4 title 'corrected  S12_r' lt 3 dt solid with lines, \
    '' index 0 using 1:5 title 'expected   S12_i' lt 4 with points, \
    '' index 1 using 1:5 title 'measured   S12_i' lt 4 dt 2 with lines, \
    '' index 2 using 1:5 title 'corrected  S12_i' lt 4 dt solid with lines, \
    '' index 0 using 1:6 title 'expected   S21_r' lt 5 with points, \
    '' index 1 using 1:6 title 'measured   S21_r' lt 5 dt 2 with lines, \
    '' index 2 using 1:6 title 'corrected  S21_r' lt 5 dt solid with lines, \
    '' index 0 using 1:7 title 'expected   S21_i' lt 6 with points, \
    '' index 1 using 1:7 title 'measured   S21_i' lt 6 dt 2 with lines, \
    '' index 2 using 1:7 title 'corrected  S21_i' lt 6 dt solid with lines, \
    '' index 0 using 1:8 title 'expected   S22_r' lt 7 with points, \
    '' index 1 using 1:8 title 'measured   S22_r' lt 7 dt 2 with lines, \
    '' index 2 using 1:8 title 'corrected  S22_r' lt 7 dt solid with lines, \
    '' index 0 using 1:9 title 'expected   S22_i' lt 8 with points, \
    '' index 1 using 1:9 title 'measured   S22_i' lt 8 dt 2 with lines, \
    '' index 2 using 1:9 title 'corrected  S22_i' lt 8 dt solid with lines
pause -1
