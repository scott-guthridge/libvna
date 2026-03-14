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
