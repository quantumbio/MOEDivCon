set MOE=C:\moe2011
set PLTFM_VER=windows-x86_64
set WD=C:\testing-grounds\quantumbio
cd %CD%
set PATH=%WD%\bin;%MOE%\bin-i4w9;%WD%\%PLTFM_VER%\bin;%WD%\%PLTFM_VER%\lib;C:\Windows\SysWOW64;%PATH%

set QBHOME=%WD%
set OPAL_PREFIX=%QBHOME%

set MOE_SVL_RUNPATH=%QBHOME%\svl
set MOE_SVL_LOAD=%QBHOME%\svl

set QBEXEC_INFO=Unknown
moe

