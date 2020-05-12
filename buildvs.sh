#!/bin/bash

# For VS, this needs to be run from a git bash shell which has also sourced the
# environment variables for Visual Studio.  To do so, start "VS2015 x64 Native
# Tools Command Prompt", and then run "C:\Program Files\Git\git-bash.exe".

# When compiled by VS, maph will run from cmd or git bash.  When compiled by
# gcc, maph will only run from bash.

use_defaultgen="true"

source ./submodules/bat/build.sh

