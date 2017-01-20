# Requirements
* Python 2.7
* Scikit Image 0.11.x (for CHT and basic image processing (gaussian filter))
* Matplotlib (for I/O)
* Numpy (for numerical applications)
* Pillow (Python Imaging Library) (for saving)

# Development Environment
For my own use, I developed/tested on Windows, with the latest Anaconda Python distribution as of Dec 2, 2015. For some reason, there an issue with their version of matplotlib and I ran "conda update matplotlib" to update to the latest matplotlib and the error was resolved.

# How to Run
To run the example code, change into the directory of interest and run the dlrse.py script.

For example
    $ cd cell_seq
    $ python ../dlrse.py
will process all the cell images in the cell_seq directory. There should be both a real-time display of the processing and output images saved to the <location of dlrse.py>/out/ folder.