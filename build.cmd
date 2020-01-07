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

%CMAKE% .. -DCMAKE_BUILD_TYPE=%BTYPE%
%CMAKE% --build . --config %BTYPE%

REM from TARGET
popd

:endbuild

