# Compositional Drift

This program tool allows polymer chemists easily to simulate the results of controlled polymerizations given all required reaction parameters (monomer amounts, monomer reactivity ratio, average polydispersity, ect.). The simulations use the Mayo-Lewis equation along Monte Carlo methods to predict polymer chain growth. The tool allows for rapid pre-experimental screening of reaction conditions along with critical post-synthesis analysis of the chain structure. Please refer to the [manuscript](https://pubs.acs.org/doi/full/10.1021/acsmacrolett.8b00813) published in Macro Letters for more details on theory, implementation, and applications.

![polymer image](https://i.imgur.com/mElYP4x.png)

## Getting Started
The program can be accessed in several ways.

### Running via Executable (Windows only)

* Download and save the .exe file in the most [current release](https://github.com/vince-wu/CompositionalDrift/releases)

* Run the program as an executable by double clicking (certification is on its way)

### Running via Web App (Demo Version)

* Go to https://vince-wu.github.io/CompositionalDrift/ 

* Note that the web app does not have the full capability of the executable program

### Running Locally with Python 

* Download and install Python 3.5+

* (Optional but recommend) Set up a python virtual environment using `virtualenv`

* Clone into your local repository:

`git clone https://github.com/vince-wu/CompositionalDrift.git`

* Navigate into the CompositionDrift directory and install all dependencies:

`pip install -r requirements.txt`

* Run the program:

`python polymerApp.py`

## Built With

* Python
* Javascript

## Contributing/ Development

The code for this application is open source and available to everyone. Feel free to clone or fork the repository if
you want to alter or add onto the codebase. You can refer to the wiki for a quick rundown of the code.

Please report any bugs to vincent.wu@berkeley.edu

## Versioning

For the versions available, see the [tags on this repository](https://github.com/vince-wu/CompositionalDrift/tags). 

## Authors

* **Vincent Wu** 

## License

This project is licensed under the MIT License - see the [LICENSE.txt](https://github.com/vince-wu/CompositionalDrift/blob/master/LISCENCE.txt) file for details

## Acknowledgments

* Anton A. A. Smith, PhD
