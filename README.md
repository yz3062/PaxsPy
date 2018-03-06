# PaxsPy
This software calculates protactinium excess, a proxy for ocean sediments usually used in conjunction with thorium excess, from ICP mass spectrometer results of uranium, thorium, and protactinium isotopes. It has a sister project, [ThxsPy](https://github.com/yz3062/ThxsPy), that calculates thorium excess.

## Getting started

### Prerequisites

If you're new to python, I recommend you simply download [Anaconda](https://www.anaconda.com/download/), which includes all the prerequisites. Alternatively, download the followings individually:

Python 2

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

![alt text](/README_screenshots/Spike_prompt_PaxsPy.JPG)

Click "Yes", and you'll then be prompted to select data files as well as a sample info file

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
