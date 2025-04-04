# atomipy Web Interface

A web application for processing molecular structure files using the atomipy package. This application allows users to upload GRO, PDB, or XYZ files, process them with either MINFF or CLAYFF forcefield atom typing, and generate various topology files for molecular dynamics simulations.

## Features

- Upload GRO, PDB, or XYZ structure files
- Automatically assign chemical elements
- Calculate bonds and angles with periodic boundary conditions
- Support for multiple forcefields:
  - MINFF forcefield for mineral systems
  - CLAYFF forcefield for clay minerals, hydroxides, and oxyhydroxides
- Multiple charge assignment methods:
  - Formal charges
  - MINFF partial charges
  - CLAYFF partial charges
  - Charge balancing
- Parallel processing of both MINFF and CLAYFF forcefields
- Generate processed structure files with the new atom types
- Create topology files for multiple simulation packages:
  - GROMACS (.itp)
  - NAMD/CHARMM (.psf)
  - LAMMPS (.data)
- Organized results with clear categorization by forcefield type

## Force Field References

### CLAYFF
Cygan, R.T.; Liang, J.J.; Kalinichev, A.G. Molecular Models of Hydroxide, Oxyhydroxide, and Clay Phases and the Development of a General Force Field. *J. Phys. Chem. B* **2004**, *108*, 1255-1266.

### MINFF
[Include MINFF reference here if available]

## Installation

1. Make sure you have Python 3.7 or newer installed
2. Clone this repository or download the files
3. Install the required dependencies:

```bash
pip install -r requirements.txt
```

4. Make sure the atomipy package is accessible (it should be in the parent directory)

## Usage

1. Start the web server:

```bash
python app.py
```

2. Open a web browser and navigate to:

```
http://localhost:5000
```

3. Upload a GRO or PDB file through the web interface
4. The application will process the file and provide links to download the results

## Project Structure

```
atomipy-web/
├── app.py                 # Main Flask application
├── requirements.txt       # Python dependencies
├── static/                # Static files (CSS, JS)
│   └── css/
│       └── style.css      # Custom styling
├── templates/             # HTML templates
│   ├── about.html         # About page
│   ├── base.html          # Base template with layout
│   ├── index.html         # Homepage with upload form
│   └── results.html       # Results page
├── uploads/               # Directory for uploaded files
└── results/               # Directory for processed files
```

## Technical Details

This web application is built with:
- Flask: A lightweight Python web framework
- Bootstrap 5: For responsive UI design
- atomipy: The molecular structure processing package

The application follows this processing workflow:
1. Import the structure file (GRO or PDB)
2. Assign elements to atoms based on their names/residues
3. Calculate bonds and angles with periodic boundary awareness
4. Assign specialized MINFF atom types based on chemical environment
5. Generate structure files with the new atom types
6. Create topology files for various simulation packages

## License

This project is released under the MIT License, the same as the atomipy package.
