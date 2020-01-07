@echo off

set BTYPE=Release

:forargs
if [%1] == [] goto :endforargs
REM echo %1

if [%1] == [Release] (
	set BTYPE=Release
) else if [%1] == [Debug] (
	set BTYPE=Debug
) else (
	echo Warning: unknown cmd argument %1
	echo.
)

shift
goto :forargs
:endforargs

echo Using build type %BTYPE%

REM goto :endbuild

REM call "C:\Program Files (x86)\Microsoft Visual Studio 12.0\vc\vcvarsall.bat" x86

set CMAKE="C:\Program Files\cmake\bin\cmake.exe"

set TARGET=target
mkdir %TARGET%
pushd %TARGET%

REM If cl is not in the path, source the environment for the most recent VS version available.
where cl > NUL 2>&1 || call "%VS190COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS180COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS170COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS160COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS150COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS130COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1 || call "%VS120COMNTOOLS%\..\..\VC\vcvarsall.bat" > NUL 2>&1

%CMAKE% .. -DCMAKE_BUILD_TYPE=%BTYPE%
%CMAKE% --build . --config %BTYPE%

REM from TARGET
popd

:endbuild

