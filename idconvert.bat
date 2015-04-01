F:
cd F:\MSGFplus\
::  run idconvert without commands to ask it for help with usage
idconvert.exe
pause
:: push any key after a pause to run remaining commands
:: run idconvert
idconvert.exe output\exactive\2015\Mg132\mg132pl_4_light1.mzid -v -e .mzid
idconvert.exe output\exactive\2015\Mg132\mg132pl_4_heavy1.mzid -v -e .mzid
pause
