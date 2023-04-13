# **OpenIBIS: Open Implementation of BIS Algorithm**

OpenIBIS is an open-source implementation of the BIS (Bispectral Index) algorithm used in anesthesia monitoring. BIS scores are commonly used to measure depth of anesthesia, but the proprietary nature of these algorithms limits our understanding of their physiological basis and clinical significance. This repository aims to provide an open and clear implementation of the BIS algorithm, enabling researchers and clinicians to better understand, develop, and improve upon existing anesthesia monitoring techniques.

## **Background**
Based on the work by Connor (2022), this project focuses on the reimplementation of the BIS algorithm to provide a clear understanding of its physiological basis and clinical significance. BIS monitors process electroencephalogram (EEG) data to provide a depth-of-anesthesia score. However, the proprietary nature of these algorithms restricts our understanding of their physiological basis, limiting advances in closed-loop control and the avoidance of postoperative cognitive dysfunction or juvenile neurological injury.

This project is based on a previous forensic disassembly of an A-2000 BIS monitor, which allowed BIS scores to be calculated from arbitrary EEG data for the first time. The current project focuses on understanding and reimplementing the BIS Engine to provide a clear-box depth-of-anesthesia algorithm.

## **Citation**
Connor, C. W. (2022). Open Reimplementation of the BIS Algorithms for Depth of Anesthesia. Anesthesia & Analgesia, 135(4), 855-864.


``bibtex
@article{connor_open_2022,
	title = {Open {Reimplementation} of the {BIS} {Algorithms} for {Depth} of {Anesthesia}},
	volume = {135},
	issn = {0003-2999},
	url = {https://journals.lww.com/10.1213/ANE.0000000000006119},
	doi = {10.1213/ANE.0000000000006119},
	language = {en},
	number = {4},
	urldate = {2023-03-13},
	journal = {Anesthesia \& Analgesia},
	author = {Connor, Christopher W.},
	month = oct,
	year = {2022},
	pages = {855--864},,
}
``

## **Installation**
To install the OpenIBIS library, please follow these instructions:

``bash
Copy code
git clone https://github.com/user/OpenIBIS.git
cd OpenIBIS
pip install -r requirements.txt
``
Usage
Please refer to the provided example scripts for usage and instructions on integrating OpenIBIS into your projects.

Contributing
We welcome contributions to improve OpenIBIS. Please submit a pull request or open an issue to contribute to the project or report bugs.

License
This project is released under the MIT License. Please see the LICENSE file for more information.