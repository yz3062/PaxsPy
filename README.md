# PaxsPy
This software calculates protactinium excess, a proxy for ocean sediments usually used in conjunction with thorium excess, from ICP mass spectrometer results of uranium, thorium, and protactinium isotopes. It has a sister project, [ThxsPy](https://github.com/yz3062/ThxsPy), that calculates thorium excess.

## Getting started

### Prerequisites

If you're new to python, I recommend you simply download [Anaconda](https://www.anaconda.com/download/), which includes all the prerequisites. Alternatively, download the followings individually:

Python 3 (the newest version of this program ended support for Python 2, since it has reached its end of life)

[scipy](https://www.scipy.org/)

[numpy](http://www.numpy.org/)

[pandas](https://pandas.pydata.org/)

### Installing

Click "Clone or download" on top right -> "Download ZIP". Unzip all files.

## Running the tests

Open Command Prompt in Windows or Terminal in MacOS, change to the directory where ThxsPy.py exsits, and type
```
python PaxsPy.py
```
You'll be prompted to confirm the version of Pa and Th spikes you're using

![alt text](/README_screenshots/Spike_prompt_PaxsPy.png)

Hit "Enter", and you'll be asked whether you want to inspect the data in figures

![alt text](/README_screenshots/inspect_figures_prompt.png)

Hit "Enter". Figures of all the ICPMS counts will be saved in the same folder as the input files. You can check the figures to see if there's anything abnormal, e.g. a spike in counts or trailing in counts midway. Notice that the all isotopes are plotted on the same y-axis, meaning you'll mostly see the variations in major isotopes like 238U and 232Th. You'll then select data files as well as a sample info file. Notice that the file selector window sometimes doesn't pop up and is open in the background.

![alt text](/README_screenshots/data_selection_PaxsPy.JPG)

Where you should double click "data" folder and select all the files in that folder

![alt text](/README_screenshots/data_select_all_PaxsPy.JPG)

Notice that alongside the data files, there's also a "sample_info.xlsx" file that looks like this

![alt text](/README_screenshots/sample_info_screenshot.JPG)

Next, in the command window you'll be prompted to enter the days between Pa spike preparation and Pa clean up column

![alt text](/README_screenshots/decay_days_PaxsPy.JPG)

And Voila! Calculation is done and you're asked to save the output file

![alt text](/README_screenshots/save_output_PaxyPy.JPG)

## License

[BSD](https://opensource.org/licenses/BSD-2-Clause)
