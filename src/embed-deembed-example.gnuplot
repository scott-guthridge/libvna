#!/usr/bin/gnuplot
#
#  Vector Network Analyzer Library
#  Copyright © 2020-2026 D Scott Guthridge <scott_guthridge@rompromity.net>
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
# Plot script for embed-deembed-example
#
set title 'Embed Example'
set xlabel 'Frequency (Hz)'
set ylabel 'Reflection Coefficient'
plot 'embed-deembed-example.out' \
       index 0 using 1:2 title 'original S11_r' lt 1 dt solid with lines, \
    '' index 0 using 1:3 title 'original S11_i' lt 1 dt 2 with lines, \
    '' index 0 using 1:4 title 'embedded S11_r' lt 2 dt solid with lines, \
    '' index 0 using 1:5 title 'embedded S11_i' lt 2 dt 2 with lines
pause -1

set title 'De-Embed Example'
set xlabel 'Frequency (Hz)'
set ylabel 'Reflection Coefficient'
plot  [] [-1:1] 'embed-deembed-example.out' \
       index 1 using 1:2 title 'measured S11_r' lt 3 dt solid with lines, \
    '' index 1 using 1:3 title 'measured S11_i' lt 3 dt 2 with lines, \
    '' index 1 using 1:4 title 'de-embedded S11_r' lt 4 dt solid with lines, \
    '' index 1 using 1:5 title 'de-embedded S11_i' lt 4 dt 2 with lines
pause -1
