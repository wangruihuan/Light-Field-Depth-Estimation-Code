
rem 3x3 views with color pixels
rprf.exe -ndm_jndm_auto_fit 1 clean -i stillLife_lf_3x3 -dstep 0.5 -datap 4.000000 6.000000 -smoothp 300.000000 -0.01 9.0 24.0 0.05 1.5  

rem five crosshair views with color pixels
rprf.exe -+ -ndm_jndm_auto_fit 1 clean -i stillLife_lf_3x3 -dstep 0.5 -datap 4.000000 6.000000 -smoothp 300.000000 -0.01 9.0 24.0 0.05 1.5  

rem five crosshair views with grey-scale pixels
rprf.exe -+ -luma -ndm_jndm_auto_fit 1 clean -i stillLife_lf_3x3 -dstep 0.5 -datap 4.000000 6.000000 -smoothp 300.000000 -0.01 9.0 24.0 0.05 1.5  

