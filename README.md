## Table of content
- [Installation](#installation)
    * [Prerequisites](#prerequisites)
    * [Packages](#packages)
    * [Scripts](#scripts)
- [Usage](#usage)
    * [Examples](#examples)
- [Contact](#contact)

## Installation
The installation of the Chimera Detector is rather straight forward, enter this command in a command line: ```git clone https://github.com/IJsbarnd35/ChimeraDetector``` Now the tool is ready for use if the correct packages discussed below are installed and present in the correct version.
### Prerequisites
For the Chimera Detector Python 3.8 is required, it can be a higher version, but not lower. Furthermore, MiniMap2 is a requirement for the mapping method. MiniMap2 can be installed [here](https://github.com/lh3/minimap2), a full installation guide is present.
For the Chimera Detector to run smoothly MiniMap must be present in the same folder, otherwise the subprocess statement will not find it.
### Packages
| Python Library | Version | Usage within the script          |
|----------------|---------|----------------------------------|
| bioPython      | 1.81    | Parsing the FastQ data set.      |
| docopt         | 0.6.2   | Creating command line interface. |
| os             | 3.10    | Finding paths to the data set.   |
| pandas         | 2.0.1   | Parsing PAF data set.            |
| random         | 3.8     | To create atrificial data.       |
| re             | 3.8     | Creation of regular expressions. |
| subprocess     | 3.10    | Calling on Minimap2.             |
| unittest       | 3.10    | Building tests for script.       |
### Scripts
The Chimera detector runs with the chimera_detector.py, this script calls upon three separate scripts to align the reads, map these and detect the chimeras and remove these. These functions can be run separately. 
## Usage
The script can be called by using the following statements: -f for the fastq file, -m for the method, -s for the slice size and -l for the minimum length of the reads that should be checked.
The parameters used are visible in the table below:

| Parameter | Required | Default                                                               |
|-----------|----------|-----------------------------------------------------------------------|
|    -f     |   Yes    |   Is required                                                         |
|    -m     |   No     |   Is required, either , "mapping", "self_aligned" or "both"           |
|    -s     |   No     |   75                                                                  |
|    -l     |   No     |   2000                                                                |
Some info on parameters here.
 
### Examples
Some examples for usage of the different methods are:<br />
1: ``python3 chimera_detector.py -f "/path/to/fastq/file.fq" -m "mapping"``
2: ``python3 chimera_detector.py -f "/path/to/fastq/file.fq" -m "self_aligned" -s 100 -l 1000``


The first statement does a regular run with the mapping method. The second statement does a run with the self_aligned method using slice sizes of 100 and with a minimum length of 1000.
## Contact
For any questions or issues please contact the creator as specified below.
* IJsbrand Pool
  * i.j.a.pool@st.hanze.nl