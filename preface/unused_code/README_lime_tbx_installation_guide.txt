# Unfortunately, you cannot install lime_tbx just by typing pip install lime_tbx into the command prompt.
# Here's how you can install lime_tbx as a Python package:

1. Download and setup Git: https://git-scm.com/downloads/
2. Type the following commands into your command prompt

git clone https://github.com/LIME-ESA/lime_tbx.git
cd lime_tbx
pip install -e ".[pyside6]"


3. lime_tbx should be successfully installed afterwards.
4. Clone the needed kernels and eocfi data folders to %appdata% from the following commands

cd %APPDATA%
git clone https://github.com/lime-esa/lime_tbx.git LimeTBX


5. If your Python interpreter still doesn't detect lime_tbx imports, try:

pip install git+https://github.com/lime-esa/lime_tbx.git


# For Mac and Linux users:
Try running the pipeline once with Moon_Noise = 'Y_Moon'.
InputCheck should alert you where the kernels and eocfi folders should be at for your respective device.
Then, replace %APPDATA% in step 4 with the appropriate path (From %APPDATA%\LimeTBX for Windows).